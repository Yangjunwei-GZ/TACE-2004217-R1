"""Draw manuscript and supplementary figures exclusively from completed revision outputs."""
from pathlib import Path
import os,json,gzip,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from PIL import Image
import SimpleITK as sitk
from publication_data import LABELS
from validation_integrity import audit_validation, verify_paired_summary
ROOT=Path(os.environ.get('TACE_PROJECT_ROOT','.'))/'revision_R1_20261014'; OUT=ROOT/'05_figures_and_tables/Figures'
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':2,'lines.linewidth':2,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
BLUE='#477B96'; GREEN='#4C997A'; RED='#B75D65'; GRAY='#ADB7BC'; GOLD='#B4945B'

def save(fig,name):
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ['png','tiff','pdf','svg']: fig.savefig(OUT/f'{name}.{ext}',dpi=600,bbox_inches='tight',facecolor='white')
    for ext in ['png','tiff']:
        p=OUT/f'{name}.{ext}'; im=Image.open(p).convert('RGB'); im.save(p,dpi=(600,600),**({'compression':'tiff_lzw'} if ext=='tiff' else {}))
    plt.close(fig)

def panel(ax,letter): ax.text(-.08,1.05,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,va='bottom')

def draw_interval_point(ax,estimate,lower,upper,y,color):
    """Draw a percentile interval and point independently; neither is clipped."""
    if not np.isfinite([estimate,lower,upper]).all() or lower>upper:
        raise ValueError('A finite ordered interval is required for publication')
    ax.plot([lower,upper],[y,y],color=color,lw=1.5)
    ax.plot([lower,lower],[y-.06,y+.06],color=color,lw=1.5)
    ax.plot([upper,upper],[y-.06,y+.06],color=color,lw=1.5)
    ax.plot([estimate],[y],marker='o',color=color,linestyle='none')

def main():
    validation=audit_validation(ROOT)
    verify_paired_summary(ROOT,validation)
    data=pd.read_csv(ROOT/'03_model_development/FINAL_REVISION_PATIENT_SCORES.csv'); s=json.loads((ROOT/'03_model_development/final_analysis_summary.json').read_text()); audit=json.loads((ROOT/'01_data_audit/audit_summary.json').read_text())
    perf=pd.read_csv(ROOT/'04_validation/full_pipeline_performance.csv').set_index('model'); delta=pd.read_csv(ROOT/'04_validation/paired_incremental_performance.csv')
    eligibility=pd.read_csv(ROOT/'01_data_audit/phase_analysis_eligibility.csv'); pairs=pd.read_csv(ROOT/'01_data_audit/phase_pairing_audit.csv'); time=pd.read_csv(ROOT/'05_figures_and_tables/cross_validated_time_performance.csv'); cal=pd.read_csv(ROOT/'05_figures_and_tables/cross_validated_calibration.csv')
    legends={}
    provenance_path=ROOT.parent/'ACTUAL_RUN_PROVENANCE.json'
    provenance=json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
    qc_scope=('Previously confirmed physician QC and17-mask work retain their physician attribution. The eight original masks from P120/P125 were retained without contour changes following a physician decision communicated by the author. This attribution correction supplies no new signed or independent review record and no confirmation that L03/L04 are two independent clinical lesions. The physician decision date, decision-specific independence and blinding remain unprovided.')
    source_clinical=pd.read_csv(ROOT/'01_data_audit/source_clinical_233.csv')
    arterial=pairs[pairs.ct_phase==1];candidate_n=arterial.patient_id.nunique();candidate_masks=len(arterial)
    fig,ax=plt.subplots(figsize=(7.09,6.1)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    boxes=[(.02,.82,.58,.14,f'Reconstructed source: {len(source_clinical)} patients / {len(pairs)} masks\nPrevious subset: 195 patients / 318 pairs\n38 patients / 59 pairs recovered from incomplete archive'),(.02,.59,.58,.14,f'Source-labelled arterial masks\n{candidate_n} patients / {candidate_masks} source records'),(.69,.59,.30,.14,f'Other source phases\n{len(source_clinical)-candidate_n} patients / {len(pairs)-candidate_masks} masks\nNot in primary analysis'),(.02,.36,.58,.14,f'Current eligible arterial input\n{s["n"]} patients / {s["lesions"]} active records\n107 Original features per lesion'),(.69,.36,.30,.14,f'Outside current primary input\n{candidate_n-s["n"]} patients\n{candidate_masks-s["lesions"]} inactive source records\nReasons and aliases in audits'),(.02,.08,.58,.18,f'Primary arterial analysis: {s["n"]} patients\n{s["events"]} deaths; {s["censored"]} censored\nLargest lesion by physical volume\nTraining-only development and paired validation')]
    for x,y,w,h,text in boxes:
        ax.add_patch(Rectangle((x,y),w,h,facecolor='#F2F5F6',edgecolor='#677C85',linewidth=2)); ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=8)
    for ya,yb in [(.82,.73),(.59,.50),(.36,.26)]:ax.annotate('',(.31,yb),(.31,ya),arrowprops={'arrowstyle':'->','lw':2,'color':'#677C85'})
    for start,end in [((.31,.775),(.84,.73)),((.31,.545),(.84,.50))]:ax.annotate('',end,start,arrowprops={'arrowstyle':'->','connectionstyle':'angle,angleA=0,angleB=90,rad=0','lw':2,'color':'#677C85'})
    save(fig,'Figure_1'); legends['Figure 1']='Source-cohort reconstruction and arterial-phase analysis. Patient and lesion counts derive from source audits and the actual frozen input. Archive omission is distinguished from analysis exclusion and computational aliases. '+qc_scope

    names=['clinical','radiomics_only','clinical_diameter','clinical_diameter_mean_hu','clinical_diameter_mean_hu_radiomics']
    fig,(a,b)=plt.subplots(1,2,figsize=(7.09,4.0),gridspec_kw={'width_ratios':[1,1]})
    short=['Clinical','Radiomics','Clinical + diameter','Clinical + diameter\n+ mean intensity','Clinical + diameter\n+ mean intensity + radiomics']
    for i,name in enumerate(names):
        r=perf.loc[name]; draw_interval_point(a,r.cross_validated_C_index,r.CI_lower,r.CI_upper,i,BLUE)
    a.set_yticks(range(len(names)),short); a.invert_yaxis(); a.set_xlabel('Full-pipeline cross-validated\nC-index'); panel(a,'A')
    labels=['Clinical','+ Diameter','+ Selected volume','+ Total volume','+ Diameter and mean','+ Diameter and 10th percentile']
    for i,(_,r) in enumerate(delta.iterrows()): draw_interval_point(b,r.paired_delta_C_index,r.CI_lower,r.CI_upper,i,GREEN)
    b.set_yticks(range(len(labels)),labels); b.invert_yaxis(); b.axvline(0,color='#777777',ls='--',lw=2); b.set_xlabel('Paired radiomics increment\nin C-index'); panel(b,'B')
    fig.subplots_adjust(left=.26,right=.98,top=.93,bottom=.18,wspace=1.45); save(fig,'Figure_2')
    legends['Figure 2']='Complete-development internal validation. (A) Pair-weighted grouped five-fold cross-validated C-index from predefined run 0. (B) Paired C-index increments from the same run when adding a training-developed radiomics score to each matched clinical/measurement baseline. Bars show percentile intervals from 500 patient bootstrap repetitions of the same randomized full cross-validation estimator; they are not intervals for the ten-repeat mean, fixed-score optimism corrections or OOB distribution ranges.'

    hi=data.radiomics_score>data.radiomics_score.median(); fig,ax=plt.subplots(figsize=(7.09,5.0)); fits=[]
    for group,color,name in [(~hi,BLUE,'Lower apparent score'),(hi,RED,'Higher apparent score')]:
        km=KaplanMeierFitter().fit(data.loc[group,'os_months'],data.loc[group,'death_event'],label=name); km.plot_survival_function(ax=ax,ci_show=True,color=color,show_censors=True); fits.append(km)
    ax.set_xlabel('Months from first TACE'); ax.set_ylabel('Overall survival probability'); ax.set_xlim(0,72); ax.set_ylim(0,1.02); ax.legend(loc='upper right',fontsize=8)
    k=s['apparent_KM']; ax.text(.03,.07,f'Apparent median split\nHR {k["HR"]:.2f} (95% CI {k["CI_lower"]:.2f}-{k["CI_upper"]:.2f})\nLog-rank P = {k["logrank_P"]:.3g}',transform=ax.transAxes,fontsize=8.5)
    add_at_risk_counts(*fits,ax=ax,rows_to_show=['At risk'],xticks=[0,12,24,36,48,60,72]); fig.subplots_adjust(bottom=.30); save(fig,'Figure_3')
    legends['Figure 3']='Descriptive apparent survival stratification using the full-cohort fitted revision score and its median. Numbers at risk and conditional HR with 95% CI are shown. This same-cohort display is not independent validation, and no alternative cut-point was searched.'

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.09,3.5))
    for name,color,label in [('clinical',BLUE,'Clinical'),('radiomics_only',GOLD,'Radiomics'),('clinical_diameter_mean_hu',GRAY,'Clinical + size/intensity'),('clinical_diameter_mean_hu_radiomics',GREEN,'Clinical + size/intensity + radiomics')]:
        g=time[time.model==name].sort_values('time'); ax1.plot(g.time,g.AUC,marker='o',color=color,label=label); ax2.plot(g.time,g.Brier,marker='o',color=color,label=label)
    ax1.set_ylabel('Cross-validated IPCW AUC'); ax2.set_ylabel('Cross-validated IPCW Brier score')
    for letter,ax in [('A',ax1),('B',ax2)]: panel(ax,letter); ax.set_xticks([12,24,36]); ax.set_xlabel('Months from first TACE')
    fig.legend(*ax1.get_legend_handles_labels(),loc='lower center',ncol=2,fontsize=8,bbox_to_anchor=(.5,-.13)); fig.subplots_adjust(left=.10,right=.97,bottom=.21,wspace=.5,top=.91); save(fig,'Figure_4')
    legends['Figure 4']='Internal time-dependent discrimination and prediction error from outer-fold predictions. Each model and Breslow baseline survival were fitted using training patients, and censoring weights were estimated from the same training set. Means are across evaluable folds/repeats; unsupported estimates are retained in the audit. These are internally cross-validated, not external validation.'

    chosen=data.sort_values('patient_id').iloc[[0,len(data)//2]]; fig,axes=plt.subplots(2,2,figsize=(7.09,6.5))
    for i,(_,r) in enumerate(chosen.iterrows()):
        im=sitk.DICOMOrient(sitk.ReadImage(r.ct_path),'LPS'); ma=sitk.ReadImage(r.mask_path); mapped=sitk.Resample(ma,im,sitk.Transform(),sitk.sitkNearestNeighbor,0,sitk.sitkUInt8); ct=sitk.GetArrayFromImage(im); m=sitk.GetArrayFromImage(mapped)>0; z=np.argmax(m.sum((1,2)))
        for j in range(2):
            ax=axes[i,j]; ax.imshow(ct[z],cmap='gray',vmin=-125,vmax=225); ax.axis('off'); ax.text(.01,.99,'ABCD'[i*2+j],transform=ax.transAxes,color='white',weight='bold',va='top',fontsize=11)
            ax.text(.02,.5,'R',transform=ax.transAxes,color='white',va='center'); ax.text(.98,.5,'L',transform=ax.transAxes,color='white',va='center',ha='right')
            if j: ax.contour(m[z],levels=[.5],colors=[RED],linewidths=1.3)
    fig.subplots_adjust(wspace=.02,hspace=.02); save(fig,'Figure_5')
    legends['Figure 5']='Representative arterial-analysis CT inputs and corresponding final analysis outlines, chosen by source identifier without selecting survival outcomes or scores. Original CT intensities and physical-coordinate overlays are shown in paired panels. '+qc_scope

    fig,(a,b)=plt.subplots(1,2,figsize=(7.09,3.8)); a.bar(['Successful','Failed'],[len(pd.read_csv(ROOT/'02_extraction/arterial_lesion_features.csv')),len(pd.read_csv(ROOT/'02_extraction/extraction_failures.csv'))],color=[BLUE,RED]); a.set_ylabel('Eligible arterial lesion pairs'); panel(a,'A')
    vals=[107,int((~pairs.geometry_tiny_tolerance.fillna(False)).sum()),int((~eligibility.phase_analysis_eligible).sum())]; b.barh(['Original features','Physical geometry\nrequires review','Patients outside\nprimary phase analysis'],vals,color=[BLUE,GRAY,GOLD]); b.set_xlabel('Count'); panel(b,'B'); fig.subplots_adjust(left=.1,right=.97,wspace=.8,bottom=.2,top=.9); save(fig,'Supplementary_Figure_S1')
    legends['Supplementary Figure S1']='Radiomics extraction and eligibility quality control. Counts include the final eligible arterial lesion set; geometry and phase exclusions are documented separately and are not silently discarded.'

    fig,axes=plt.subplots(1,3,figsize=(7.09,2.9))
    for ax,t,letter in zip(axes,[12,24,36],'ABC'):
        ax.plot([0,1],[0,1],color='#777777',ls='--',lw=1.5)
        for model,col,label in [('clinical',BLUE,'Clinical'),('radiomics_only',GOLD,'Radiomics'),('clinical_diameter_mean_hu',GRAY,'Clinical + size/intensity'),('clinical_diameter_mean_hu_radiomics',GREEN,'Clinical + size/intensity + radiomics')]:
            g=cal[(cal.model==model)&(cal.time_months==t)]; ax.plot(g.mean_predicted_risk,g.observed_km_risk,'o-',color=col,label=label)
        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xlabel(f'Mean predicted risk at {t} months'); ax.set_ylabel('Observed KM event risk'); panel(ax,letter)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,fontsize=8,bbox_to_anchor=(.5,-.23)); fig.subplots_adjust(wspace=.55,bottom=.24,left=.09,right=.98,top=.9); save(fig,'Supplementary_Figure_S2')
    legends['Supplementary Figure S2']='Cross-validated internal grouped calibration. Patient-level mean validation event probabilities over ten repeats are grouped into five quantiles and compared with group Kaplan-Meier event risks. This is an internal descriptive calibration assessment, with group uncertainty and small-sample limitations; it is not external calibration.'

    fig,ax=plt.subplots(figsize=(7.09,4.2)); corr=pd.read_csv(ROOT/'03_model_development/SIZE_INTENSITY_COLLINEARITY.csv'); selected=corr[(corr.feature=='original_shape_Maximum2DDiameterRow')&(corr.baseline=='selected_diameter_mm')].iloc[0]
    ax.scatter(data.selected_diameter_mm,data.original_shape_Maximum2DDiameterRow,color=BLUE,alpha=.75,s=15); ax.set_xlabel('Selected-lesion maximum 3D diameter (mm)'); ax.set_ylabel('Maximum 2D diameter, row (mm)'); ax.text(.03,.93,f'Spearman rho = {selected.Spearman_rho:.3f}',transform=ax.transAxes); save(fig,'Supplementary_Figure_S3')
    legends['Supplementary Figure S3']='Overlap between a historical size-related radiomics feature and maximum tumour diameter. This new sensitivity figure replaces the retired exploratory decision curve, which lacked a justified clinical action and thresholds.'

    records=pd.read_csv(ROOT/'04_validation/feature_recurrence_records.csv'); group=records[records['mode']=='point'].groupby('feature').selected.sum().sort_values(ascending=True).tail(12)/50
    fig,ax=plt.subplots(figsize=(7.09,4.2)); ax.barh(group.index.str.replace('original_',''),group,color=BLUE); ax.set_xlabel('Recurrence in the 50 outer training fits'); ax.set_xlim(0,1); fig.subplots_adjust(left=.47,right=.97,bottom=.15,top=.97); save(fig,'Supplementary_Figure_S4')
    legends['Supplementary Figure S4']='Feature recurrence across the 50 training fits used for repeated outer-fold point estimation. Correlated alternatives, shrinkage and variable selection are reflected in these frequencies; no feature-level biological stability claim is made.'

    boot=pd.read_csv(ROOT/'04_validation/validation_run_metrics.csv'); g=boot[(boot['mode']=='bootstrap')&(boot.model=='clinical_diameter_mean_hu_radiomics')&(boot.successful_folds==5)]
    fig,ax=plt.subplots(figsize=(7.09,3.5)); ax.hist(g.C_index,bins=25,color=GREEN,edgecolor='white'); ax.set_xlabel('Bootstrap full-pipeline cross-validation C-index'); ax.set_ylabel('Patient bootstrap repetitions'); save(fig,'Supplementary_Figure_S5')
    legends['Supplementary Figure S5']='Distribution of the full cross-validation performance statistic under patient bootstrap sampling. This is the sampling distribution used for the stated estimator intervals and is distinct from the historical per-resample OOB performance distribution.'

    selected=['clinical_diameter','clinical_selected_volume','clinical_total_volume']; fig,(a,b)=plt.subplots(1,2,figsize=(7.09,3.7))
    for i,name in enumerate(selected):
        r=perf.loc[name]; rr=perf.loc[name+'_radiomics']; a.plot([r.cross_validated_C_index,rr.cross_validated_C_index],[i,i],color=GRAY,lw=2); a.scatter(r.cross_validated_C_index,i,color=BLUE,label='Baseline' if i==0 else None); a.scatter(rr.cross_validated_C_index,i,color=GREEN,label='+ radiomics' if i==0 else None)
    a.set_yticks(range(3),['Maximum diameter','Selected physical volume','Total segmented volume']); a.set_xlabel('Full-pipeline cross-validated C-index'); a.legend(loc='lower right',fontsize=8); panel(a,'A')
    for i,name in enumerate(selected):
        r=delta[delta.baseline==name].iloc[0]; b.errorbar(r.paired_delta_C_index,i,xerr=[[r.paired_delta_C_index-r.CI_lower],[r.CI_upper-r.paired_delta_C_index]],fmt='o',color=GREEN,capsize=3)
    b.axvline(0,color='#777777',ls='--',lw=2); b.set_yticks(range(3),['Diameter','Selected volume','Total volume']); b.set_xlabel('Paired radiomics increment\nin C-index (95% CI)'); panel(b,'B'); fig.subplots_adjust(left=.25,right=.98,bottom=.22,top=.92,wspace=.85); save(fig,'Supplementary_Figure_S6')
    legends['Supplementary Figure S6']='Paired tumour-burden comparisons under the revised development protocol. (A) Same-run cross-validated performance for matched clinical plus size/physical-volume baselines, with and without a newly training-developed radiomics score. (B) Paired radiomics increments with full-pipeline bootstrap intervals. All estimates correspond to the same arterial cohort and validation design; historical independently resampled corrected values are not combined here.'

    fig,(a,b)=plt.subplots(1,2,figsize=(7.09,3.7)); counts=pairs.groupby('phase_name').patient_id.nunique(); a.barh(counts.index.str.replace('_',' ').str.capitalize(),counts,color=BLUE); a.set_xlabel('Source patients with a tumour mask'); panel(a,'A')
    included=eligibility.phase_analysis_eligible.astype(bool); no_source=eligibility.source_arterial_masks.eq(0)
    groups=['Included','No source arterial mask','Other documented exclusion']; counts=[int(included.sum()),int((~included&no_source).sum()),int((~included&~no_source).sum())]
    b.barh(groups,counts,color=[GREEN,GRAY,GRAY]); b.set_xlabel('Patients'); panel(b,'B'); fig.subplots_adjust(left=.15,right=.98,bottom=.20,top=.9,wspace=1.4); save(fig,'Supplementary_Figure_S7')
    legends['Supplementary Figure S7']='Official source phase distribution and actual patient-level primary-analysis eligibility. Source labels and recorded final phase entries remain separate from research eligibility decisions; no phase-invariant performance claim is made. '+qc_scope
    (ROOT/'05_figures_and_tables/FIGURE_LEGENDS.json').write_text(json.dumps(legends,indent=2),encoding='utf-8')
    print('Figures 1-5 and supplementary figures S1-S7 exported.',flush=True)

if __name__=='__main__':main()
