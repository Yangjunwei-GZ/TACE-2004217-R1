"""Provenance-checked extraction, original-grid physical volume, and arterial patient matrix."""
from pathlib import Path
import json, logging, platform, importlib.metadata, hashlib, warnings, os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import radiomics
from radiomics import featureextractor
import nibabel as nib

ROOT=Path(os.environ.get('TACE_PROJECT_ROOT',str(Path.cwd()))); REV=ROOT/'revision_R1_20261014'; OUT=REV/'02_extraction'

def extractor():
    ex=featureextractor.RadiomicsFeatureExtractor(binWidth=25,resampledPixelSpacing=[1,1,1],interpolator='sitkBSpline',normalize=False,preCrop=False,padDistance=5,force2D=False,force2Ddimension=0,label=1,additionalInfo=True)
    ex.disableAllImageTypes(); ex.enableImageTypeByName('Original'); ex.disableAllFeatures()
    for cls in ['firstorder','shape','glcm','glrlm','glszm','gldm','ngtdm']: ex.enableFeatureClassByName(cls)
    return ex

def execute(row,ex):
    im=sitk.ReadImage(row.ct_path); ma=sitk.ReadImage(row.mask_path); a=sitk.GetArrayFromImage(ma); binary=sitk.GetImageFromArray((a>0).astype(np.uint8)); binary.CopyInformation(ma)
    result=ex.execute(im,binary,label=1)
    feat={k:float(v) for k,v in result.items() if k.startswith('original_')}
    diag={k:str(v) for k,v in result.items() if k.startswith('diagnostics_')}
    roi=sitk.GetArrayFromImage(im)[a>0]
    n=nib.load(row.ct_path); proxy=n.dataobj
    diag.update(nifti_proxy_slope=str(proxy.slope),nifti_proxy_intercept=str(proxy.inter),raw_roi_mean=str(float(np.mean(roi))),raw_roi_median=str(float(np.median(roi))),raw_roi_10percentile=str(float(np.percentile(roi,10))),raw_roi_min=str(float(np.min(roi))),raw_roi_max=str(float(np.max(roi))))
    return feat,diag

def main():
    OUT.mkdir(exist_ok=True)
    logging.getLogger('radiomics').setLevel(logging.ERROR)
    a=pd.read_csv(REV/'01_data_audit/phase_pairing_audit.csv',dtype={'patient_id':str})
    p=pd.read_csv(REV/'01_data_audit/phase_analysis_eligibility.csv',dtype={'patient_id':str})
    c=pd.read_csv(REV/'01_data_audit/source_clinical_233.csv',dtype={'patient_id':str})
    a=a[(a.ct_phase==1)&a.eligible&(a.patient_id.isin(p.loc[p.phase_analysis_eligible,'patient_id']))].sort_values('sample_id').copy()
    legacy=ROOT/'results_corrected_radiomics/20260816_212013/corrected_radiomics_features_20260816_212013.csv'
    if legacy.exists():
        old=pd.read_csv(legacy).set_index('sample_id'); fcols=[x for x in old if x.startswith('original_')]
    else:
        fcols=sorted('original_'+cls+'_'+name for cls,klass in radiomics.getFeatureClasses().items() if cls in ['firstorder','shape','glcm','glrlm','glszm','gldm','ngtdm'] for name,deprecated in klass.getFeatureNames().items() if not deprecated)
        old=pd.DataFrame(columns=fcols)
    assert len(fcols)==107
    ex=extractor(); versions={n:importlib.metadata.version(n) for n in ['numpy','pandas','SimpleITK','pyradiomics-fix','scikit-survival','scikit-learn','scipy']}
    config=dict(python=platform.python_version(),radiomics_module=radiomics.__version__,distribution=versions,settings=ex.settings,imageType=ex.enabledImagetypes,featureClass=ex.enabledFeatures,mask_interpolation='sitkNearestNeighbor used by PyRadiomics resampleImage',geometryTolerance='not overridden; SimpleITK defaults',correctMask=False,source_mask_modification='positive labels binarised in memory only')
    config.update(SimpleITK_coordinate_tolerance=sitk.ProcessObject.GetGlobalDefaultCoordinateTolerance(),SimpleITK_direction_tolerance=sitk.ProcessObject.GetGlobalDefaultDirectionTolerance(),voxelArrayShift=0,voxelArrayShift_evidence='RadiomicsFirstOrder.__init__: kwargs.get(voxelArrayShift, 0), inspected in installed 3.0.1 module')
    (OUT/'extraction_effective_settings.json').write_text(json.dumps(config,indent=2,default=str),encoding='utf-8')
    yaml='imageType:\n  Original: {}\nfeatureClass:\n'+''.join('  '+cl+': []\n' for cl in ex.enabledFeatures)+'setting:\n'+''.join('  '+str(k)+': '+json.dumps(v)+'\n' for k,v in ex.settings.items())
    yaml+='  voxelArrayShift: 0\n  correctMask: false\n  geometryTolerance: null\n# Mask interpolation: nearest neighbour in physical coordinates.\n# SimpleITK coordinate and direction tolerances: 1e-6, not overridden.\n'
    (OUT/'extraction_parameters.yaml').write_text(yaml,encoding='utf-8')
    completed={}
    cp=OUT/'arterial_lesion_features_checkpoint.csv'
    if cp.exists(): completed=pd.read_csv(cp,dtype={'patient_id':str}).set_index('sample_id').to_dict('index')
    verification=[]; reused=a[a.sample_id.isin(old.index)]
    candidates=reused.iloc[np.unique(np.array([0,len(reused)//2,len(reused)-1]))] if len(reused) else reused
    for _,row in candidates.iterrows():
        features,diag=execute(row,ex)
        expected=old.loc[row.sample_id,fcols].astype(float)
        rel=np.abs(expected.to_numpy()-np.array([features[f] for f in fcols]))/np.maximum(1,np.abs(expected.to_numpy()))
        verification.append(dict(sample_id=row.sample_id,max_relative_difference=float(rel.max()),pass_comparison=bool(np.all(rel<1e-8))))
        (OUT/f'diagnostics_verification_{row.sample_id}.json').write_text(json.dumps(diag,indent=2),encoding='utf-8')
    pd.DataFrame(verification).to_csv(OUT/'corrected_feature_reproduction_check.csv',index=False)
    reuse=bool(verification) and all(r['pass_comparison'] for r in verification)
    rows=[]; failures=[]
    for i,(_,row) in enumerate(a.iterrows()):
        if row.sample_id in completed: rows.append(dict(sample_id=row.sample_id,**completed[row.sample_id])); continue
        try:
            if reuse and row.sample_id in old.index:
                feat=old.loc[row.sample_id,fcols].astype(float).to_dict(); provenance='verified corrected extraction reused'
            else:
                feat,diag=execute(row,ex); provenance='newly extracted with locked corrected settings'
                (OUT/f'diagnostics_{row.sample_id}.json').write_text(json.dumps(diag,indent=2),encoding='utf-8')
            rows.append(dict(sample_id=row.sample_id,patient_id=row.patient_id,ct_phase=1,lesion_index=int(row.lesion_index),physical_volume_mm3=float(row.physical_volume_mm3),mask_voxels=int(row.mask_voxels),ct_path=row.ct_path,mask_path=row.mask_path,feature_provenance=provenance,**feat))
        except Exception as exc: failures.append(dict(sample_id=row.sample_id,reason=str(exc)))
        if (i+1)%10==0:
            pd.DataFrame(rows).to_csv(cp,index=False); print(f'EXTRACT {i+1}/{len(a)}',flush=True)
    f=pd.DataFrame(rows); f.to_csv(OUT/'arterial_lesion_features.csv',index=False); pd.DataFrame(failures,columns=['sample_id','reason']).to_csv(OUT/'extraction_failures.csv',index=False)
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
    mat.to_csv(REV/'03_model_development/analysis_matrix.csv',index=False)
    summary=dict(n_patients=len(mat),events=int(mat.death_event.sum()),censored=int((1-mat.death_event).sum()),n_lesions=len(f),features=len(fcols),reproduction_pass=reuse,selection_changed=int(mat.lesion_selection_changed_from_voxels.sum()),extraction_failed=len(failures),phase=1)
    (OUT/'extraction_matrix_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(summary,flush=True)

if __name__=='__main__': main()
