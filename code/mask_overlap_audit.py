"""Audit whether sum of all provided lesion volumes double-counts overlapping mask voxels."""
from pathlib import Path
import os,numpy as np,pandas as pd,SimpleITK as sitk
ROOT=Path(os.environ.get('TACE_PROJECT_ROOT',str(Path.cwd())))/'revision_R1_20261014'
def main():
    f=pd.read_csv(ROOT/'02_extraction/arterial_lesion_features.csv',dtype={'patient_id':str}); rows=[]
    for pid,g in f.groupby('patient_id'):
        union=None; total=0; voxel_volume=None
        for _,r in g.iterrows():
            mask=sitk.ReadImage(r.mask_path); a=sitk.GetArrayFromImage(mask)>0; total+=int(a.sum()); voxel_volume=float(np.prod(mask.GetSpacing()))
            union=a.copy() if union is None else np.logical_or(union,a)
        rows.append(dict(patient_id=pid,n_masks=len(g),summed_voxels=total,union_voxels=int(union.sum()),overlap_voxels=total-int(union.sum()),union_volume_cm3=float(union.sum()*voxel_volume/1000),sum_volume_cm3=float(g.physical_volume_mm3.sum()/1000)))
    out=pd.DataFrame(rows); out.to_csv(ROOT/'01_data_audit/ARTERIAL_MASK_OVERLAP_AUDIT.csv',index=False); print('Patients with overlap',int((out.overlap_voxels>0).sum()),'total overlap voxels',int(out.overlap_voxels.sum()),flush=True)
if __name__=='__main__':main()
