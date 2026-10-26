"""Training-only radiomics development and paired nested patient bootstrap validation."""
import os
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_key]='1'
from pathlib import Path
import argparse, hashlib, json, gzip, time, warnings, concurrent.futures, sys, traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
from validation_integrity import (ensure_validation_plan, checkpoint_issues,
                                  audit_validation, verify_paired_summary, ValidationIntegrityError)
from validation_execution import execute_validation_tasks, ValidationTasksFailed
from survival_evidence import fitted_baseline, monthly_survival, source_fingerprint, SCHEMA_VERSION, PRIMARY_MODELS
from run_new_cohort import require_frozen_model_inputs
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import SplineTransformer
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc, brier_score, integrated_brier_score
from sksurv.util import Surv
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.functions import StepFunction

ROOT=Path(os.environ.get('TACE_PROJECT_ROOT','.'))/'revision_R1_20261014'
CFG=json.loads((ROOT/'00_protocol_and_provenance/revision_config.json').read_text())
OUT=ROOT/'04_validation'; OUT.mkdir(exist_ok=True)
CLINICAL=['age','sex_female','child_b','log_afp','log_tumour_count','bclc_b']
SPECS={
 'clinical':CLINICAL,
 'clinical_diameter':CLINICAL+['selected_diameter_mm'],
 'clinical_selected_volume':CLINICAL+['log_selected_volume'],
 'clinical_total_volume':CLINICAL+['log_total_volume'],
 'clinical_mean_hu':CLINICAL+['mean_hu'],
 'clinical_median_hu':CLINICAL+['median_hu'],
 'clinical_p10_hu':CLINICAL+['p10_hu'],
 'clinical_diameter_mean_hu':CLINICAL+['selected_diameter_mm','mean_hu'],
 'clinical_diameter_p10_hu':CLINICAL+['selected_diameter_mm','p10_hu'],
 'clinical_axial_diameter':CLINICAL+['clinical_max_axial_diameter_mm'],
 'clinical_legacy_seven':['age','sex_female','child_b','log_afp','afp_high_400','tumour_multiple','bclc_b'],
}
for _name,_cols in list(SPECS.items()):
 if _name in ['clinical','clinical_diameter','clinical_selected_volume','clinical_total_volume','clinical_diameter_mean_hu','clinical_diameter_p10_hu']: SPECS[_name+'_radiomics']=_cols+['radiomics_score']
SPECS['radiomics_only']=['radiomics_score']
SPECS['clinical_age_spline']=CLINICAL

def labels(d): return Surv.from_arrays(d.death_event.astype(bool),d.os_months.astype(float))

def cstat(d,risk):
    val=concordance_index_censored(d.death_event.astype(bool).to_numpy(),d.os_months.to_numpy(),np.asarray(risk))
    return dict(C_index=float(val[0]),concordant=int(val[1]),discordant=int(val[2]),tied_risk=int(val[3]),comparable_pairs=int(val[1]+val[2]+val[3]))

def transforms(d,cols):
    a=d[cols].to_numpy(dtype=float); med=np.nanmedian(a,axis=0); med=np.nan_to_num(med,nan=0.0); a=np.where(np.isnan(a),med,a)
    mean=a.mean(0); sd=a.std(0); keep=sd>1e-10
    return dict(columns=list(cols),median=med.tolist(),mean=mean.tolist(),sd=np.where(keep,sd,1).tolist(),keep=keep.tolist())

def apply(d,tr):
    a=d[tr['columns']].to_numpy(dtype=float); a=np.where(np.isnan(a),tr['median'],a)
    return ((a-np.array(tr['mean']))/np.array(tr['sd']))[:,np.array(tr['keep'])]

def features(d):
    cols=[c for c in d if c.startswith('original_')]
    ranked=[]
    for col in cols:
        a=d[col].to_numpy(dtype=float); a=np.where(np.isnan(a),np.nanmedian(a),a)
        if np.std(a)<1e-10: continue
        try: ranked.append((abs(cstat(d,a)['C_index']-.5),col))
        except Exception: continue
    ranked.sort(key=lambda x:(-x[0],x[1])); chosen=[]
    corr=d[[x[1] for x in ranked]].corr(method='spearman').abs()
    for _,col in ranked:
        if all(corr.loc[col,c]<=CFG['correlation_threshold'] for c in chosen): chosen.append(col)
        if len(chosen)>=CFG['max_filtered_features']: break
    if not chosen: raise ValueError('No train-variable radiomics feature')
    return chosen

def fit_radiomics(d,kind,penalty):
    chosen=features(d)
    if kind!='ridge': chosen=chosen[:2]
    tr=transforms(d,chosen); x=apply(d,tr)
    if kind=='top_two_elastic_net':
        model=CoxnetSurvivalAnalysis(l1_ratio=.5,alphas=[penalty],max_iter=100000,fit_baseline_model=True).fit(x,labels(d))
        coef=model.coef_[:,0]
    else:
        model=CoxPHSurvivalAnalysis(alpha=penalty,ties='breslow',n_iter=100).fit(x,labels(d)); coef=model.coef_
    return dict(model=model,transform=tr,kind=kind,penalty=penalty,features=chosen,coefficients=np.asarray(coef).tolist())

def rad_predict(d,fit): return fit['model'].predict(apply(d,fit['transform']))

def grouped_splits(d,folds,seed):
    ids=np.array(sorted(d.patient_id.unique())); kf=KFold(folds,shuffle=True,random_state=seed)
    for train,val in kf.split(ids):
        tr=np.flatnonzero(d.patient_id.isin(ids[train])); va=np.flatnonzero(d.patient_id.isin(ids[val])); yield tr,va

def develop(d,seed):
    settings=[(kind,float(p)) for kind in CFG['candidates'] for p in (CFG['ridge_penalties'] if kind=='ridge' else CFG['two_feature_penalties'])]
    score={s:[] for s in settings}; failure={s:[] for s in settings}
    for fold,(i,j) in enumerate(grouped_splits(d,CFG['inner_folds'],seed)):
        train=d.iloc[i]; val=d.iloc[j]
        for spec in settings:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore'); fit=fit_radiomics(train,*spec)
                stat=cstat(val,rad_predict(val,fit)); score[spec].append(stat)
            except Exception as exc: failure[spec].append({'fold':fold,'seed':seed,'exception_type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
    tuning=[]
    for spec in settings:
        st=score[spec]; total=sum(x['comparable_pairs'] for x in st)
        value=sum(x['C_index']*x['comparable_pairs'] for x in st)/total if total else None
        tuning.append(dict(kind=spec[0],penalty=spec[1],C_index=value,success_folds=len(st),failed_folds=len(failure[spec]),failure=failure[spec]))
    valid=[v for v in tuning if v['success_folds']==CFG['inner_folds'] and v['C_index'] is not None and np.isfinite(v['C_index'])]
    if not valid: raise ValueError('All inner development candidates failed')
    best=max(v['C_index'] for v in valid); shortlist=[v for v in valid if v['C_index']>=best-.01]
    chosen=sorted(shortlist,key=lambda v:(v['kind']=='ridge',-v['penalty'],-v['C_index'],v['kind']))[0]
    with warnings.catch_warnings(): warnings.simplefilter('ignore'); fit=fit_radiomics(d,chosen['kind'],chosen['penalty'])
    fit['tuning']=tuning
    return fit

def fit_baseline(train,val,cols,spline=False):
    tr=transforms(train,cols); x=apply(train,tr); z=apply(val,tr); basis=None
    if spline:
        basis=SplineTransformer(n_knots=2,degree=3,include_bias=False,extrapolation='linear')
        x=np.column_stack((basis.fit_transform(train[['age']]),x[:,1:])); z=np.column_stack((basis.transform(val[['age']]),z[:,1:]))
    if x.shape[1]==0:
        times,probability=kaplan_meier_estimator(train.death_event.astype(bool),train.os_months)
        class NullModel:
            coef_=np.array([])
            def predict(self,design): return np.zeros(len(design))
            def predict_survival_function(self,design): return [StepFunction(times,probability) for _ in range(len(design))]
        return NullModel(),z,tr,basis
    with warnings.catch_warnings():
        warnings.simplefilter('ignore'); model=CoxPHSurvivalAnalysis(alpha=CFG['baseline_penalizer'],ties='breslow',n_iter=100).fit(x,labels(train))
    return model,z,tr,basis

def time_metrics(train,val,risk,model,z):
    result=[]; survival=model.predict_survival_function(z); months=np.array(CFG['prediction_months'],dtype=float)
    max_time=min(max(labels(train)['time']),max(labels(val)['time']))
    for t in months:
        if t>=max_time or t<=min(val.os_months): result.append(dict(time=t,status='unsupported_followup')); continue
        try:
            st=labels(val).copy(); st['event'][st['time']>t]=False; st['time']=np.minimum(st['time'],t+1e-7)
            pp=np.array([f(t) for f in survival]); auc,_=cumulative_dynamic_auc(labels(train),st,risk,[t]); _,bs=brier_score(labels(train),st,pp[:,None],[t])
            result.append(dict(time=float(t),AUC=float(auc[0]),Brier=float(bs[0]),n=len(val),status='evaluated',predicted_risk=(1-pp).tolist()))
        except Exception as exc: result.append(dict(time=float(t),status='failed',error=str(exc)))
    try:
        grid=np.arange(max(1,int(np.ceil(val.os_months.min()+1e-7))),min(36,int(np.floor(max_time-1e-7)))+1,dtype=float)
        if len(grid)>=12:
            st=labels(val).copy(); st['event'][st['time']>grid[-1]]=False; st['time']=np.minimum(st['time'],grid[-1]+1e-7)
            pp=np.array([[f(t) for t in grid] for f in survival]); ibs=integrated_brier_score(labels(train),st,pp,grid)
            result.append(dict(time='IBS',IBS=float(ibs),grid=grid.tolist(),n=len(val),status='evaluated'))
    except Exception as exc: result.append(dict(time='IBS',status='failed',error=str(exc)))
    return result

def one_run(run_id,mode):
    require_frozen_model_inputs(ROOT.parent)
    plan=ensure_validation_plan(ROOT)
    path=OUT/'checkpoints'/f'{mode}_{run_id:04d}.json.gz'; path.parent.mkdir(exist_ok=True)
    if path.exists():
        with gzip.open(path,'rt',encoding='utf-8') as handle: saved=json.load(handle)
        problems=checkpoint_issues(saved,plan,path)
        if problems:raise RuntimeError('Checkpoint reuse refused: '+json.dumps(problems,ensure_ascii=False))
        return {'path':str(path),'reused':True}
    base=pd.read_csv(ROOT/'03_model_development/analysis_matrix.csv',dtype={'patient_id':str})
    seed=CFG['seed_point']+run_id if mode=='point' else CFG['seed_bootstrap']+run_id
    if mode=='bootstrap':
        rng=np.random.default_rng(seed); draw=rng.integers(0,len(base),len(base)); d=base.iloc[draw].reset_index(drop=True)
    else: draw=np.arange(len(base)); d=base.copy()
    run=dict(run_id=run_id,mode=mode,seed=seed,n=len(d),unique_patients=d.patient_id.nunique(),folds=[],failures=[],analysis_matrix_sha256=hashlib.sha256((ROOT/'03_model_development/analysis_matrix.csv').read_bytes()).hexdigest(),configuration_sha256=hashlib.sha256(json.dumps(CFG,sort_keys=True).encode()).hexdigest())
    run['configuration_file_sha256']=hashlib.sha256((ROOT/'00_protocol_and_provenance/revision_config.json').read_bytes()).hexdigest()
    run.update(prediction_evidence_schema_version=SCHEMA_VERSION,
               execution_source_sha256=source_fingerprint(Path(__file__).resolve().parent),
               draw_source_row_indices=draw.tolist(),draw_patient_ids=d.patient_id.tolist())
    t0=time.time()
    for fold,(i,j) in enumerate(grouped_splits(d,CFG['outer_folds'],seed)):
        train=d.iloc[i].copy(); val=d.iloc[j].copy(); assert not set(train.patient_id)&set(val.patient_id)
        record=dict(fold=fold,train_ids=train.patient_id.tolist(),validation_ids=val.patient_id.tolist(),models={})
        record.update(train_draw_row_indices=i.tolist(),validation_draw_row_indices=j.tolist(),
                      training_events=train.death_event.tolist(),training_times=train.os_months.tolist(),
                      inner_splits=[dict(fold=k,seed=seed+10000+fold,
                         train_outer_training_row_indices=ii.tolist(),validation_outer_training_row_indices=jj.tolist(),
                         train_ids=train.iloc[ii].patient_id.tolist(),validation_ids=train.iloc[jj].patient_id.tolist())
                         for k,(ii,jj) in enumerate(grouped_splits(train,CFG['inner_folds'],seed+10000+fold))])
        try:
            fit=develop(train,seed+10000+fold); train['radiomics_score']=rad_predict(train,fit); val['radiomics_score']=rad_predict(val,fit)
            record['radiomics_development']={k:v for k,v in fit.items() if k!='model'}
            record['radiomics_development']['fitted_baseline']=fitted_baseline(fit['model'],train.os_months)
            for name,cols in SPECS.items():
                try:
                    model,z,tr,basis=fit_baseline(train,val,cols,name=='clinical_age_spline'); risk=model.predict(z)
                    st=cstat(val,risk); record['models'][name]=dict(**st,transform=tr,coefficients=model.coef_.tolist(),risk=risk.tolist(),validation_events=val.death_event.tolist(),validation_times=val.os_months.tolist())
                    record['models'][name]['fitted_baseline']=fitted_baseline(model,train.os_months)
                    if name in PRIMARY_MODELS:
                        record['models'][name]['monthly_survival']=monthly_survival(record['models'][name]['fitted_baseline'],risk)
                    if basis is not None:
                        record['models'][name]['spline_basis']={'n_features_in':int(basis.n_features_in_),
                            'n_features_out':int(basis.n_features_out_),
                            'degree':int(basis.degree),'extrapolation':basis.extrapolation,
                            'include_bias':bool(basis.include_bias),
                            'bsplines':[dict(t=b.t.tolist(),c=b.c.tolist(),k=int(b.k)) for b in basis.bsplines_]}
                    if mode=='point' and name in ['clinical','clinical_diameter_mean_hu','clinical_diameter_mean_hu_radiomics','radiomics_only']:
                        record['models'][name]['time_metrics']=time_metrics(train,val,risk,model,z)
                except Exception as exc: record['models'][name]=dict(error=str(exc),exception_type=type(exc).__name__,traceback=traceback.format_exc())
            if len(base)>=100 and base.death_event.sum()>=70:
                for suffix,cols in [('RSF',CLINICAL+['selected_diameter_mm','mean_hu']),('radiomics_RSF',CLINICAL+['selected_diameter_mm','mean_hu','radiomics_score'])]:
                    tr=transforms(train,cols); x=apply(train,tr); z=apply(val,tr)
                    rsf=RandomSurvivalForest(n_estimators=100,min_samples_leaf=10,max_depth=3,max_features='sqrt',random_state=seed+fold,n_jobs=1).fit(x,labels(train)); risk=rsf.predict(z)
                    record['models']['clinical_size_intensity_'+suffix]=dict(**cstat(val,risk),risk=risk.tolist(),transform=tr,validation_events=val.death_event.tolist(),validation_times=val.os_months.tolist())
        except Exception as exc: run['failures'].append(dict(fold=fold,seed=seed,error=str(exc),exception_type=type(exc).__name__,traceback=traceback.format_exc()))
        run['folds'].append(record)
    summaries={}
    for name in list(SPECS)+['clinical_size_intensity_RSF','clinical_size_intensity_radiomics_RSF']:
        stats=[f['models'].get(name,{}) for f in run['folds']]; valid=[s for s in stats if 'C_index' in s]
        total=sum(s['comparable_pairs'] for s in valid)
        summaries[name]=dict(C_index=sum(s['C_index']*s['comparable_pairs'] for s in valid)/total if total else None,successful_folds=len(valid),attempted_folds=CFG['outer_folds'])
    run['summary']=summaries; run['elapsed_seconds']=time.time()-t0
    with gzip.open(path,'wt',encoding='utf-8') as f: json.dump(run,f,allow_nan=True)
    problems=checkpoint_issues(run,plan,path)
    if problems:raise RuntimeError('Incomplete validation checkpoint retained for failure audit: '+json.dumps(problems,ensure_ascii=False))
    return {'path':str(path),'reused':False}

def summarize():
    # Never manufacture a formal summary from a subset of successful runs.
    integrity=audit_validation(ROOT)
    records=[]; failures=[]; selections=[]
    for p in sorted((OUT/'checkpoints').glob('*.json.gz')):
        with gzip.open(p,'rt',encoding='utf-8') as f: r=json.load(f)
        for model,s in r['summary'].items(): records.append(dict(run=r['run_id'],mode=r['mode'],model=model,**s))
        failures+= [dict(run=r['run_id'],mode=r['mode'],**x) for x in r['failures']]
        for f in r['folds']:
            if 'radiomics_development' in f:
                for name,beta in zip(f['radiomics_development']['features'],f['radiomics_development']['coefficients']): selections.append(dict(mode=r['mode'],run=r['run_id'],fold=f['fold'],feature=name,coefficient=beta,selected=abs(beta)>1e-8))
    df=pd.DataFrame(records); df.to_csv(OUT/'validation_run_metrics.csv',index=False)
    pd.DataFrame(failures,columns=['run','mode','fold','error']).to_csv(OUT/'validation_failures.csv',index=False)
    pd.DataFrame(selections).to_csv(OUT/'feature_recurrence_records.csv',index=False)
    summary=[]
    for name,g in df.groupby('model'):
        pt=g[(g['mode']=='point')&(g.successful_folds==5)].C_index.dropna(); bs=g[(g['mode']=='bootstrap')&(g.successful_folds==5)].C_index.dropna()
        primary=g[(g['mode']=='point')&(g.run==0)&(g.successful_folds==5)].C_index
        if name not in integrity['plan']['model_names']:continue
        counts=integrity['counts']['bootstrap']
        summary.append(dict(model=name,cross_validated_C_index=primary.iloc[0] if len(primary) else np.nan,repeat_mean_C_index=pt.mean(),repeat_SD_C_index=pt.std(),point_repeats_completed=len(pt),bootstrap_planned=counts['planned'],bootstrap_attempted=counts['attempted'],bootstrap_failed=counts['failed'],bootstrap_complete=len(bs),CI_lower=bs.quantile(.025),CI_upper=bs.quantile(.975),estimator='predefined run 0 pair-weighted grouped 5-fold CV; percentile patient bootstrap of the same randomized full CV procedure; ten-repeat mean is split sensitivity only'))
    pd.DataFrame(summary).to_csv(OUT/'full_pipeline_performance.csv',index=False)
    paired=[]
    for base in ['clinical','clinical_diameter','clinical_selected_volume','clinical_total_volume','clinical_diameter_mean_hu','clinical_diameter_p10_hu']:
        combo=base+'_radiomics'; sub=df[df.model.isin([base,combo])]; wide=sub.pivot(index=['mode','run'],columns='model',values='C_index')
        ok=sub[sub.successful_folds==5].groupby(['mode','run']).size(); wide=wide.loc[ok[ok==2].index]; delta=wide[combo]-wide[base]; pt=delta.loc['point']; bs=delta.loc['bootstrap']
        paired.append(dict(baseline=base,combined=combo,paired_delta_C_index=pt.loc[0],repeat_mean_delta=pt.mean(),repeat_SD_delta=pt.std(),CI_lower=bs.quantile(.025),CI_upper=bs.quantile(.975),bootstrap_joint_complete=len(bs),point_joint_complete=len(pt),ci_crosses_zero=bs.quantile(.025)<=0<=bs.quantile(.975)))
    pd.DataFrame(paired).to_csv(OUT/'paired_incremental_performance.csv',index=False)
    verify_paired_summary(ROOT,integrity)
    print(pd.DataFrame(summary).to_string(index=False),flush=True); print(pd.DataFrame(paired).to_string(index=False),flush=True)

def main():
    require_frozen_model_inputs(ROOT.parent)
    parser=argparse.ArgumentParser(); parser.add_argument('--mode',choices=['point','bootstrap','summary'],required=True); parser.add_argument('--workers',type=int,default=4); parser.add_argument('--start',type=int,default=0); parser.add_argument('--stop',type=int); args=parser.parse_args()
    if args.mode=='summary': summarize(); return
    limit=args.stop if args.stop is not None else CFG['point_repeats'] if args.mode=='point' else CFG['bootstrap_replicates']
    result=execute_validation_tasks(ROOT,args.mode,range(args.start,limit),one_run,workers=args.workers)
    print(json.dumps(result,ensure_ascii=False),flush=True)
    # Only the explicit summary entry creates formal tables after checking all
    # point and bootstrap tasks, model specifications, folds and paired results.

if __name__=='__main__':
    try:main()
    except Exception as exc:
        mode=sys.argv[sys.argv.index('--mode')+1] if '--mode' in sys.argv and sys.argv.index('--mode')+1<len(sys.argv) else None
        from validation_integrity import _fingerprints
        failure={'status':'FAILED','mode':mode,'run_id':None,'seed':None,
                 'scope':'validation_entry_exception; per-task seeds and run IDs are retained in VALIDATION_TASK_FAILURES.jsonl',
                 'exception':{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()},
                 'input_fingerprints':_fingerprints(ROOT),'arguments':sys.argv[1:]}
        with (OUT/'VALIDATION_ENTRY_FAILURES.jsonl').open('a',encoding='utf-8') as handle:
            handle.write(json.dumps(failure,ensure_ascii=False)+'\n')
        print('VALIDATION_FAILED',str(exc),file=sys.stderr,flush=True)
        sys.exit(1)
