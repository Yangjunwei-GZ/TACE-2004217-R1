"""Audit true WAW inputs and recover missing same-phase images into the revision only."""
from pathlib import Path
import argparse, csv, hashlib, json, struct, zipfile, zlib, datetime, urllib.request, os
import numpy as np
import pandas as pd
import SimpleITK as sitk

ROOT = Path(os.environ.get('TACE_PROJECT_ROOT',str(Path.cwd())))
REV = ROOT / 'revision_R1_20261014'
OUT = REV / '01_data_audit'
OLD = ROOT / 'results_corrected_radiomics/20260816_212013'

def digest(p, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def dump(p,v):
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,default=str),encoding='utf-8')

def salvage_index(path):
    """Index complete local ZIP entries without modifying a truncated original archive."""
    records=[]
    size=path.stat().st_size
    with path.open('rb') as f:
        while f.tell()+30<=size:
            offset=f.tell(); b=f.read(30)
            if b[:4]!=b'PK\x03\x04': break
            _,ver,flag,method,tm,dt,crc,comp,uncomp,nl,el=struct.unpack('<4s5H3I2H',b)
            name=f.read(nl).decode('utf-8' if flag&2048 else 'cp437'); extra=f.read(el)
            if flag&8: break
            if comp==0xffffffff or uncomp==0xffffffff:
                pos=0
                while pos+4<=len(extra):
                    ident,length=struct.unpack_from('<HH',extra,pos); value=extra[pos+4:pos+4+length]; pos+=length+4
                    if ident==1:
                        k=0
                        if uncomp==0xffffffff: uncomp=struct.unpack_from('<Q',value,k)[0]; k+=8
                        if comp==0xffffffff: comp=struct.unpack_from('<Q',value,k)[0]
            start=f.tell()
            if start+comp>size: break
            records.append(dict(name=name,offset=start,size=comp,uncompressed=uncomp,crc=crc,method=method))
            f.seek(start+comp)
    return records

def headers(p):
    r=sitk.ImageFileReader(); r.SetFileName(str(p)); r.ReadImageInformation()
    return dict(size=r.GetSize(),spacing=r.GetSpacing(),origin=r.GetOrigin(),direction=r.GetDirection())

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=[]
    paths=[ROOT/'data_raw/WAW_TACE/clinical_data_wawtace_v2_15_07_2024.xlsx',ROOT/'data_raw/WAW_TACE/ct_hcc_metadata_v2.csv',ROOT/'data_raw/WAW_TACE/supplementary_table_s1_definitions_v2.xlsx',OLD/'corrected_radiomics_features_20260816_212013.csv',OLD/'corrected_radiomics_extraction_qc_20260816_212013.csv']
    paths+=list((ROOT/'code/corrected_radiomics').glob('*.py'))
    paths+=list((REV/'00_protocol_and_provenance/handoff_inputs').glob('submitted_inputs/*'))
    for p in paths:
        if not p.exists(): continue
        item=dict(path=str(p),sha256=digest(p),bytes=p.stat().st_size,mtime=p.stat().st_mtime,role='old_result' if 'results_corrected' in str(p) else 'source_input')
        if p.suffix.lower()=='.csv':
            df=pd.read_csv(p); item.update(rows=len(df),columns=list(df),identifiers=[c for c in df if c in ['patient_id','sample_id','PATPRI']])
        manifest.append(item)
    dump(REV/'00_protocol_and_provenance/input_manifest.json',manifest)
    clinical=pd.read_excel(paths[0]); clinical['patient_id']=clinical.PATPRI.astype(int).astype(str)
    meta=pd.read_csv(paths[1]); meta.patient_id=meta.patient_id.astype(str)
    old=pd.read_csv(paths[3],dtype={'patient_id':str}) if paths[3].exists() else pd.DataFrame(columns=['patient_id','sample_id','ct_path'])
    definitions=pd.read_excel(paths[2]); definitions.to_csv(OUT/'source_clinical_definitions.csv',index=False)
    cv=[]
    for c in clinical:
        if c=='patient_id': continue
        definition=definitions.loc[definitions.Variable==c,'Description']
        cv.append(dict(variable=c,definition=definition.iloc[0] if len(definition) else 'Not provided',n=233,missing=int(clinical[c].isna().sum()),unique=clinical[c].nunique(),range_or_levels=str(clinical[c].dropna().unique()[:15].tolist()),post_treatment=c in ['initial_tace_answer','tace_number','initial_LR_TR','if_nonV','progression','progression_time','nonv','time_to_nonv','death','survival_time']))
    pd.DataFrame(cv).to_csv(OUT/'clinical_variable_audit.csv',index=False)
    archives=[]; valid={}; salvage={}
    for p in sorted((ROOT/'data_raw/WAW_TACE').glob('ct_scans*.zip')):
        ok=zipfile.is_zipfile(p)
        if ok:
            z=zipfile.ZipFile(p); valid[p]=z
            names=z.namelist()
        else:
            records=salvage_index(p); salvage[p]={x['name']:x for x in records}; names=list(salvage[p]); dump(OUT/(p.stem+'_complete_local_entries.json'),records)
        archives.append(dict(path=str(p),bytes=p.stat().st_size,central_directory_readable=ok,complete_local_entries=len(names),contains_patients=sorted(set(x.split('/')[0] for x in names if '/' in x))))
    dump(OUT/'archive_audit.json',archives)
    maskzip=zipfile.ZipFile(ROOT/'data_raw/WAW_TACE/tumor_masks_wawtace_v1_08_05_2024.zip')
    all_masknames=[x for x in maskzip.namelist() if x.endswith('.nrrd') and Path(x).name[0].isdigit()]
    masknames=[x for x in all_masknames if x.endswith('_tumor_seg.nrrd')]
    dump(OUT/'noncanonical_mask_entries.json',[dict(name=x,reason='Noncanonical additional archive entry; not counted as an official tumour mask') for x in all_masknames if x not in masknames])
    cache={Path(p).name:Path(p) for p in old.ct_path.unique()}
    cache.update({p.name:p for p in (ROOT/'data_processed').glob('waw*ct*/**/*.nii.gz')})
    recovered=REV/'02_extraction/source_recovered'; recovered.mkdir(parents=True,exist_ok=True)
    checkpoint=OUT/'phase_pairing_audit_checkpoint.csv'
    completed=pd.read_csv(checkpoint,dtype={'patient_id':str}).set_index('sample_id').to_dict('index') if checkpoint.exists() else {}
    rows=[]
    for i,name in enumerate(sorted(masknames)):
        stem=Path(name).name.removesuffix('_tumor_seg.nrrd'); pid,phase,li=stem.split('_'); phase=int(phase)
        if stem in completed and 'CT unavailable' not in str(completed[stem].get('exclusion_reason','')):
            rows.append(dict(sample_id=stem,**completed[stem])); continue
        mask=ROOT/'data_processed/waw_tumor_masks'/name
        if not mask.exists():
            mask=recovered/'masks'/Path(name).name; mask.parent.mkdir(exist_ok=True)
            if not mask.exists(): mask.write_bytes(maskzip.read(name))
        scanname=f'{pid}_{phase}_scan.nii.gz'; scan=cache.get(scanname)
        source='existing_cache'; archive_member=f'{pid}/{scanname}'; archive_path=''
        if scan is None or not scan.exists():
            scan=recovered/'ct'/scanname; scan.parent.mkdir(exist_ok=True)
            if not scan.exists():
                found=False
                for ap,z in valid.items():
                    if archive_member in z.namelist():
                        scan.write_bytes(z.read(archive_member)); source='valid_archive_recovered'; archive_path=str(ap); found=True; break
                if not found:
                    for ap,index in salvage.items():
                        if archive_member not in index: continue
                        info=index[archive_member]
                        with ap.open('rb') as f: f.seek(info['offset']); b=f.read(info['size'])
                        raw=b if info['method']==0 else zlib.decompress(b,-15)
                        if zlib.crc32(raw)&0xffffffff!=info['crc']: raise RuntimeError(f'CRC failure in salvage {archive_member}')
                        scan.write_bytes(raw); source='crc_verified_truncated_archive_recovery'; archive_path=str(ap); found=True; break
                if not found: scan=None; source='missing_in_local_complete_entries'
            else: source='revision_cache_resumed'
        m=meta[(meta.patient_id==pid)&(meta.ct_phase==phase)]
        row=dict(sample_id=stem,patient_id=pid,ct_phase=phase,phase_name={0:'unenhanced',1:'arterial',2:'portal_venous',3:'delayed'}[phase],lesion_index=int(li),metadata_match=len(m)==1,source_ct=scanname,ct_path=str(scan or ''),mask_path=str(mask),cache_source=source,archive_path=archive_path,archive_member=archive_member,in_old_318=stem in set(old.sample_id),eligible=False,exclusion_reason='')
        try:
            if scan is None: raise FileNotFoundError('CT unavailable in local cache or recoverable archive')
            ch=headers(scan); mh=headers(mask); image_mask=sitk.ReadImage(str(mask)); arr=sitk.GetArrayFromImage(image_mask); nz=np.argwhere(arr>0)
            row.update(ct_sha256=digest(scan),mask_sha256=digest(mask),ct_size=str(ch['size']),mask_size=str(mh['size']),ct_spacing=str(ch['spacing']),mask_spacing=str(mh['spacing']),ct_origin=str(ch['origin']),mask_origin=str(mh['origin']),ct_direction=str(ch['direction']),mask_direction=str(mh['direction']),mask_labels=str(np.unique(arr).tolist()),mask_voxels=len(nz),physical_volume_mm3=float(len(nz)*np.prod(mh['spacing'])))
            if len(nz)==0: raise ValueError('Empty mask')
            delta={k:float(np.max(np.abs(np.asarray(ch[k])-np.asarray(mh[k])))) for k in ['spacing','origin','direction']}
            same_size=ch['size']==mh['size']; strict=same_size and all(v==0 for v in delta.values()); tiny=same_size and all(v<=1e-4 for v in delta.values())
            corners=[tuple(int(v) for v in xyz) for xyz in [nz.min(0)[::-1],nz.max(0)[::-1]]]
            points=[image_mask.TransformIndexToPhysicalPoint(c) for c in corners]
            transform=np.asarray(ch['direction']).reshape(3,3)@np.diag(ch['spacing'])
            indices=[np.linalg.solve(transform,np.asarray(pt)-ch['origin']) for pt in points]
            in_bounds=all(np.all(x>=-0.5) and np.all(x<=np.asarray(ch['size'])-0.5) for x in indices)
            row.update(geometry_exact=strict,geometry_tiny_tolerance=tiny,geometry_max_spacing_difference=delta['spacing'],geometry_max_origin_difference=delta['origin'],geometry_max_direction_difference=delta['direction'],roi_in_image_physical_bounds=in_bounds)
            if not tiny: raise ValueError('Nontrivial CT/mask geometry requires independent review before eligibility')
            if not in_bounds: raise ValueError('ROI outside image physical bounds')
            row['eligible']=bool(len(m)==1)
            if len(m)!=1: row['exclusion_reason']='No unique official metadata match'
        except Exception as exc: row['exclusion_reason']=str(exc)
        rows.append(row)
        if (i+1)%20==0:
            pd.DataFrame(rows).to_csv(OUT/'phase_pairing_audit_checkpoint.csv',index=False)
            print(f'AUDIT {i+1}/{len(masknames)}',flush=True)
    audit=pd.DataFrame(rows); audit.to_csv(OUT/'phase_pairing_audit.csv',index=False)
    patient=[]
    for _,c in clinical.iterrows():
        pid=c.patient_id; a=audit[audit.patient_id==pid]; arterial=a[a.ct_phase==1]
        expected=int(meta[(meta.patient_id==pid)&(meta.ct_phase==1)].tumor_count.sum())
        eligible=expected>0 and len(arterial)==expected and bool(arterial.eligible.all())
        reason='Included' if eligible else ('No source arterial-phase tumour mask' if expected==0 else '; '.join(sorted(set(arterial.loc[~arterial.eligible,'exclusion_reason'].tolist()))) or 'Mask count differs from official metadata')
        patient.append(dict(patient_id=pid,in_old_195=pid in set(old.patient_id),old_exclusion_reason='Included' if pid in set(old.patient_id) else 'Not extracted in old local cache',source_masks=len(a),source_arterial_masks=expected,correctly_paired_arterial_masks=int(arterial.eligible.sum()),phase_analysis_eligible=eligible,primary_exclusion_reason=reason,os_days=float(c.survival_time),death_event=int(c.death)))
    pat=pd.DataFrame(patient); pat.to_csv(OUT/'phase_analysis_eligibility.csv',index=False)
    summary=dict(source_patients=len(clinical),source_masks=len(masknames),official_mask_count=int(meta.tumor_count.sum()),old_patients=old.patient_id.nunique(),old_pairs=len(old),phase_counts=audit.groupby('ct_phase').size().to_dict(),correctly_paired=int(audit.eligible.sum()),eligible_arterial_patients=int(pat.phase_analysis_eligible.sum()),eligible_arterial_events=int(pat.loc[pat.phase_analysis_eligible,'death_event'].sum()),exclusion_counts=pat.primary_exclusion_reason.value_counts().to_dict(),independent_radiologist_review='source_intake_only_no_new_review_assumed',DICOM_files=0)
    dump(OUT/'audit_summary.json',summary)
    clinical.to_csv(OUT/'source_clinical_233.csv',index=False)
    dump(REV/'00_protocol_and_provenance/CHECKPOINT_DATA_AUDIT.json',dict(completed=True,timestamp=datetime.datetime.now().isoformat(),summary=summary))
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=='__main__': main()
