"""Meaningful leakage/unit/risk tests, final training fit, and reporting diagnostics."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json, hashlib, time, warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test, logrank_test
from scipy.stats import spearmanr
import validate_pipeline as vp

ROOT=vp.ROOT; OUT=ROOT/'03_model_development'

def main():
    data=pd.read_csv(OUT/'analysis_matrix.csv',dtype={'patient_id':str}); tests=[]
    tr,va=next(vp.grouped_splits(data,5,20261014)); train=data.iloc[tr].copy(); valid=data.iloc[va].copy()
    assert not set(train.patient_id)&set(valid.patient_id)
    tests.append(dict(test='patient_disjoint_outer_split',status='PASS'))
    first=vp.develop(train,20261014); transformed=vp.apply(valid,first['transform']); altered=valid.copy(); altered['death_event']=1-altered.death_event; altered['os_months']=altered.os_months.iloc[::-1].to_numpy()
    second=vp.develop(train,20261014)
    frozen=lambda f:json.dumps({k:v for k,v in f.items() if k!='model'},sort_keys=True,allow_nan=True)
    assert frozen(first)==frozen(second); assert np.allclose(vp.rad_predict(valid,first),vp.rad_predict(altered,second))
    tests.append(dict(test='validation_outcome_perturbation_does_not_change_training_or_score',status='PASS'))
    assert np.array_equal(transformed,vp.apply(valid,second['transform']))
    tests.append(dict(test='training_transform_and_seed_reproducible',status='PASS'))
    toy=pd.DataFrame(dict(patient_id=['a','b','c','d'],death_event=[1,1,1,1],os_months=[1.,2.,3.,4.]))
    assert vp.cstat(toy,np.array([4.,3.,2.,1.]))['C_index']==1.; assert vp.cstat(toy,np.array([1.,2.,3.,4.]))['C_index']==0.
    tests.append(dict(test='Harrell_risk_direction_and_event_coding',status='PASS'))
    a=pd.read_csv(ROOT/'01_data_audit/phase_pairing_audit.csv',dtype={'patient_id':str}); selected=a.set_index('sample_id').loc[data.sample_id]
    expected=selected.physical_volume_mm3.to_numpy()/1000
    assert np.allclose(data.selected_volume_cm3,expected)
    tests.append(dict(test='original_grid_mm3_to_cm3_volume',status='PASS'))
    pd.DataFrame(tests).to_csv(ROOT/'04_validation/critical_integrity_tests.csv',index=False)
    fit=vp.develop(data,20261014); data['radiomics_score']=vp.rad_predict(data,fit)
    freeze=[]; keep=np.array(fit['transform']['keep']); names=np.array(fit['transform']['columns'])[keep]
    for col,beta in zip(names,fit['coefficients']):
        pos=fit['transform']['columns'].index(col)
        freeze.append(dict(feature=col,beta=beta,mean=fit['transform']['mean'][pos],SD_population=fit['transform']['sd'][pos],median_imputation=fit['transform']['median'][pos],kind=fit['kind'],penalty=fit['penalty']))
    pd.DataFrame(freeze).to_csv(OUT/'FINAL_REVISION_MODEL_FREEZE.csv',index=False)
    (OUT/'final_development_tuning.json').write_text(json.dumps({k:v for k,v in fit.items() if k!='model'},indent=2,allow_nan=True),encoding='utf-8')
    data.to_csv(OUT/'FINAL_REVISION_PATIENT_SCORES.csv',index=False)
    reconstructed=sum(r['beta']*(data[r['feature']]-r['mean'])/r['SD_population'] for r in freeze)
    offset=getattr(fit['model'],'offset_',np.array([0.]))[0]
    assert np.allclose(data.radiomics_score,reconstructed-offset,atol=1e-8)
    k=KaplanMeierFitter().fit(data.os_months,event_observed=data.death_event)
    reverse=KaplanMeierFitter().fit(data.os_months,event_observed=1-data.death_event)
    median=float(k.median_survival_time_); follow=float(reverse.median_survival_time_)
    summary=dict(n=len(data),events=int(data.death_event.sum()),censored=int((1-data.death_event).sum()),phase='arterial',lesions=int(data.n_segmented_arterial_lesions.sum()),median_KM_OS=median,median_reverse_KM_followup=follow,observed_followup_range=[float(data.os_months.min()),float(data.os_months.max())],final_kind=fit['kind'],final_penalty=fit['penalty'],features=names.tolist(),score_offset=float(offset),apparent_radiomics_C=vp.cstat(data,data.radiomics_score)['C_index'])
    cutoff=float(data.radiomics_score.median()); hi=data.radiomics_score>cutoff
    km=logrank_test(data.loc[hi,'os_months'],data.loc[~hi,'os_months'],data.loc[hi,'death_event'],data.loc[~hi,'death_event'])
    kmdata=data[['os_months','death_event']].copy(); kmdata['high_risk']=hi.astype(int)
    kmcox=CoxPHFitter().fit(kmdata,'os_months','death_event')
    summary['apparent_KM']=dict(cutoff=cutoff,high_n=int(hi.sum()),low_n=int((~hi).sum()),HR=float(np.exp(kmcox.params_.iloc[0])),CI_lower=float(np.exp(kmcox.confidence_intervals_.iloc[0,0])),CI_upper=float(np.exp(kmcox.confidence_intervals_.iloc[0,1])),logrank_P=float(km.p_value))
    ph=[]; coef=[]
    for label,cols in [('clinical',vp.CLINICAL),('primary_combined',vp.SPECS['clinical_diameter_mean_hu_radiomics'])]:
        tr=vp.transforms(data,cols); z=vp.apply(data,tr); z=pd.DataFrame(z,columns=np.array(cols)[np.array(tr['keep'])]); z['os_months']=data.os_months; z['death_event']=data.death_event
        cox=CoxPHFitter(penalizer=.1).fit(z,'os_months','death_event')
        for transform in ['rank','km']:
            test=proportional_hazard_test(cox,z,time_transform=transform)
            for variable,row in test.summary.iterrows(): ph.append(dict(model=label,variable=variable,time_transform=transform,test_statistic=row.test_statistic,P=row.p,scope='supplementary penalised Cox assumption diagnostic; not the predictive performance estimator'))
        for variable,row in cox.summary.iterrows(): coef.append(dict(model=label,variable=variable,HR_per_training_SD=row['exp(coef)'],CI_lower=row['exp(coef) lower 95%'],CI_upper=row['exp(coef) upper 95%'],P=row.p,conditional_on_selected_score=True))
    pd.DataFrame(ph).to_csv(OUT/'PH_DIAGNOSTICS.csv',index=False); pd.DataFrame(coef).to_csv(OUT/'CONDITIONAL_ASSOCIATION_DIAGNOSTICS.csv',index=False)
    corr=[]
    for x in ['original_shape_Maximum2DDiameterRow','selected_diameter_mm','selected_volume_cm3','total_segmented_volume_cm3','clinical_max_axial_diameter_mm','p10_hu','mean_hu']:
        for y in ['selected_diameter_mm','selected_volume_cm3','mean_hu']:
            rho,p=spearmanr(data[x],data[y]); corr.append(dict(feature=x,baseline=y,Spearman_rho=rho,P=p))
    pd.DataFrame(corr).to_csv(OUT/'SIZE_INTENSITY_COLLINEARITY.csv',index=False)
    (OUT/'final_analysis_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=True),encoding='utf-8'); print(summary,flush=True)

if __name__=='__main__': main()
