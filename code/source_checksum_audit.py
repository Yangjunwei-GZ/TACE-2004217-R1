"""Hash source archives and confirm actual cache bytes against official ZIP member CRCs."""
from pathlib import Path
import os,hashlib,json,zipfile,io,zlib
import pandas as pd
ROOT=Path(os.environ.get('TACE_PROJECT_ROOT',str(Path.cwd()))); REV=ROOT/'revision_R1_20261014'

def main():
    metadata=json.loads((ROOT/'results_corrected_radiomics/20260816_212013/stage6_1_qims_representative_ct_addendum_20260826_162621/00_input_authority/WAW_TACE_ZENODO_METADATA_20260826.json').read_text())
    expected={x['key']:x['checksum'].removeprefix('md5:') for x in metadata['files']}
    rows=[]
    files=[p for p in (ROOT/'data_raw/WAW_TACE').iterdir() if p.name in expected]
    for p in files:
        sha=hashlib.sha256(); md5=hashlib.md5()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(8*1024*1024),b''): sha.update(b); md5.update(b)
        row=dict(path=str(p),bytes=p.stat().st_size,sha256=sha.hexdigest(),MD5=md5.hexdigest(),official_MD5=expected[p.name],MD5_match=md5.hexdigest()==expected[p.name],official_record='https://zenodo.org/records/12741586')
        rows.append(row); pd.DataFrame(rows).to_csv(REV/'01_data_audit/source_archive_checksum_audit.csv',index=False); print('HASHED',p.name,row['MD5_match'],flush=True)
    indices={}
    for p in files:
        if p.name.startswith('ct_scans') and zipfile.is_zipfile(p):
            z=zipfile.ZipFile(p); indices.update({Path(x.filename).name:(str(p),x.CRC) for x in z.infolist() if x.filename.endswith('.nii.gz')})
    tail=REV/'01_data_audit/part1_remote_tail.bin'
    if tail.exists():
        z=zipfile.ZipFile(io.BytesIO(tail.read_bytes())); indices.update({Path(x.filename).name:('official_part1_remote_central_directory',x.CRC) for x in z.infolist() if x.filename.endswith('.nii.gz')})
    a=pd.read_csv(REV/'01_data_audit/phase_pairing_audit.csv'); checks=[]
    for path in a.ct_path.dropna().unique():
        if not path or not Path(path).exists(): continue
        crc=0
        with Path(path).open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''): crc=zlib.crc32(b,crc)
        origin,authority=indices[Path(path).name]
        checks.append(dict(ct_path=path,crc32=crc&0xffffffff,official_crc32=authority,source_archive=origin,matched=(crc&0xffffffff)==authority))
    pd.DataFrame(checks).to_csv(REV/'01_data_audit/cache_to_official_archive_CRC_audit.csv',index=False)
    print('CT CRC matches',sum(x['matched'] for x in checks),'/',len(checks),flush=True)

if __name__=='__main__':main()
