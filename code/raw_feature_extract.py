"""Fresh image extraction using the frozen extractor and unchanged matrix equations."""
from pathlib import Path
import os,sys,json,argparse,concurrent.futures,logging,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import SimpleITK as sitk
from extract_matrix import extractor,execute
from frozen_input_pipeline import PACKAGE,sha,rows,write,dump,verify_bindings
ROOT=Path(os.environ['TACE_PROJECT_ROOT']);REV=ROOT/'revision_R1_20261014';OUT=REV/'02_extraction'
def worker(row):
    from types import SimpleNamespace
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)
    logging.getLogger('radiomics').setLevel(logging.ERROR)
    before={k:sha(row[k]) for k in ['ct_path','mask_path']}
    if before['ct_path']!=row['ct_sha256'] or before['mask_path']!=row['mask_sha256']:raise ValueError('Input hash differs before extraction: '+row['sample_id'])
    features,diag=execute(SimpleNamespace(**row),extractor())
    if len(features)!=107 or not all(np.isfinite(v) for v in features.values()):raise ValueError('Failed/nonfinite/incomplete 107-feature extraction')
    if any(sha(row[k])!=value for k,value in before.items()):raise ValueError('Input changed while extracting: '+row['sample_id'])
    return features,diag
def build_matrix(a,f,c):
    fcols=[x for x in f if x.startswith("original_")]
    counts=a.groupby('patient_id').size(); successful=f.groupby('patient_id').size(); eligible=[pid for pid,n in counts.items() if successful.get(pid,0)==n]
    selected=f[f.patient_id.isin(eligible)].sort_values(['patient_id','physical_volume_mm3','sample_id'],ascending=[True,False,True]).drop_duplicates('patient_id')
    sums=f.groupby('patient_id').physical_volume_mm3.sum()
    totals={}
    for pid,group in f.groupby('patient_id'):
        union=None; reference=None
        for _,lesion in group.iterrows():
            mask=sitk.ReadImage(lesion.mask_path)
            grid=(mask.GetSize(),mask.GetSpacing(),mask.GetOrigin(),mask.GetDirection())
            if reference is None: reference=grid
            assert grid[0]==reference[0] and all(np.allclose(grid[k],reference[k],atol=1e-4,rtol=0) for k in [1,2,3]), 'Mask union requires a common verified physical grid'
            binary=sitk.GetArrayFromImage(mask)>0
            union=binary.copy() if union is None else np.logical_or(union,binary)
        totals[pid]=float(union.sum()*np.prod(reference[1]))
    voxel_selected=f.sort_values(['patient_id','mask_voxels','sample_id'],ascending=[True,False,True]).drop_duplicates('patient_id').set_index('patient_id').sample_id
    mat=c.merge(selected,on='patient_id',how='inner',validate='one_to_one')
    mat['os_months']=mat.survival_time/30.4375; mat['death_event']=mat.death.astype(int)
    mat['selected_volume_cm3']=mat.physical_volume_mm3/1000; mat['total_segmented_volume_cm3']=mat.patient_id.map(totals)/1000
    mat['summed_mask_volume_cm3_before_overlap_correction']=mat.patient_id.map(sums)/1000
    mat['selected_diameter_mm']=mat.original_shape_Maximum3DDiameter
    mat['clinical_max_axial_diameter_mm']=mat[['lesion1_diameter','lesion2_diameter','lesion3_diameter']].max(axis=1)
    mat['mean_hu']=mat.original_firstorder_Mean; mat['median_hu']=mat.original_firstorder_Median; mat['p10_hu']=mat.original_firstorder_10Percentile
    mat['sex_female']=mat.gender_woman; mat['child_b']= (mat.cps==2).astype(int); mat['bclc_b']=(mat.bclc==1).astype(int)
    mat['log_afp']=np.log1p(mat.lab_afp); mat['log_tumour_count']=np.log1p(mat.lesions_number)
    mat['log_selected_volume']=np.log1p(mat.selected_volume_cm3); mat['log_total_volume']=np.log1p(mat.total_segmented_volume_cm3)
    mat['afp_high_400']=(mat.lab_afp>=400).astype(int); mat['tumour_multiple']=(mat.lesions_number>=2).astype(int)
    mat['lesion_selection_changed_from_voxels']=mat.sample_id!=mat.patient_id.map(voxel_selected)
    mat['n_segmented_arterial_lesions']=mat.patient_id.map(successful)
    assert not mat.patient_id.duplicated().any(); assert (mat.os_months>0).all(); assert np.isfinite(mat[fcols]).all().all()
    expected_columns=list(rows(PACKAGE/'code/data/analysis_matrix.csv')[0])
    if set(mat.columns)!=set(expected_columns):raise ValueError('Rebuilt matrix column schema differs')
    expected_ids=[r['patient_id'] for r in rows(PACKAGE/'code/data/analysis_matrix.csv')]
    if set(mat.patient_id)!=set(expected_ids):raise ValueError('Rebuilt patient set differs from frozen inputs')
    mat=mat.set_index('patient_id',drop=False).loc[expected_ids].reset_index(drop=True)
    mat[expected_columns].to_csv(REV/'03_model_development/analysis_matrix.csv',index=False)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    if any(OUT.iterdir()):raise ValueError('Fresh extraction directory must be empty; no feature checkpoints are resumed')
    binding=json.loads((ROOT/'FROZEN_RUN_INPUT_BINDINGS.json').read_text(encoding='utf-8'));verify_bindings(binding)
    a=pd.read_csv(REV/'01_data_audit/phase_pairing_audit.csv',dtype={'patient_id':str})
    active={r['sample_id'] for r in rows(ROOT/'FINAL_IMAGE_AND_ALIAS_ROSTER.csv') if r['active_analysis_record']=='True'}
    a=a[a.sample_id.isin(active)].sort_values('sample_id').copy()
    if len(a)!=len(active) or a.sample_id.duplicated().any():raise ValueError('Frozen active ROI set differs from pairing rows')
    c=pd.read_csv(REV/'01_data_audit/source_clinical_233.csv',dtype={'patient_id':str})
    results={};failures=[];source=a.to_dict('records')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(worker,r):r for r in source}
        for future in concurrent.futures.as_completed(futures):
            r=futures[future]
            try:
                feat,diag=future.result();results[r['sample_id']]=dict(sample_id=r['sample_id'],patient_id=r['patient_id'],ct_phase=1,lesion_index=int(r['lesion_index']),physical_volume_mm3=float(r['physical_volume_mm3']),mask_voxels=int(r['mask_voxels']),ct_path=r['ct_path'],mask_path=r['mask_path'],feature_provenance='newly extracted with locked corrected settings',**feat)
                dump(OUT/('diagnostics_'+r['sample_id']+'.json'),diag)
            except Exception as exc:failures.append({'sample_id':r['sample_id'],'reason':str(exc),'exception':type(exc).__name__,'ct_sha256':r['ct_sha256'],'mask_sha256':r['mask_sha256']})
            if (len(results)+len(failures))%10==0:print('EXTRACT',len(results)+len(failures),'/',len(source),'failed',len(failures),flush=True)
    write(OUT/'extraction_failures.csv',[{'sample_id':r['sample_id'],'reason':r['reason']} for r in failures],['sample_id','reason'])
    dump(OUT/'EXTRACTION_TASK_AUDIT.json',{'planned':len(source),'attempted':len(source),'success':len(results),'failed':len(failures),'failures':failures,'old_features_reused':False,'model_fits':0})
    if failures:raise RuntimeError('Feature extraction failed; no matrix or training is authorized')
    columns=list(rows(PACKAGE/'code/reporting_inputs/02_extraction/arterial_lesion_features.csv')[0]);f=pd.DataFrame([results[r['sample_id']] for r in source],columns=columns)
    fcols=[x for x in columns if x.startswith('original_')]
    if len(fcols)!=107 or not np.isfinite(f[fcols]).all().all():raise ValueError('Incomplete or nonfinite feature matrix')
    f.to_csv(OUT/'arterial_lesion_features.csv',index=False);build_matrix(a,f,c)
    dump(OUT/'extraction_effective_settings.json',{'settings':extractor().settings,'imageType':extractor().enabledImagetypes,'featureClass':extractor().enabledFeatures,'source':'same execute()/extractor() as the frozen scientific run','all_features_newly_extracted':True,'SimpleITK_threads_per_worker':1})
    verify_bindings(binding)
    print(json.dumps({'status':'FRESH_FEATURES_WRITTEN_PENDING_INPUT_COMPARISON','ROIs':len(f),'features':107,'new_model_fits':0}),flush=True)
if __name__=='__main__':main()
