"""Recover only missing official ZIP members; validate per-member CRC and preserve source ZIP."""
from pathlib import Path
import concurrent.futures, io, json, struct, subprocess, zipfile, zlib, hashlib
import pandas as pd

ROOT=Path(str(Path.cwd())); REV=ROOT/'revision_R1_20261014'; AUDIT=REV/'01_data_audit'
TAIL_START=10358430958
URL='https://zenodo.org/records/12741586/files/ct_scans_1_4_wawtace_09_05_24.zip?download=1'

def main():
    index=zipfile.ZipFile(io.BytesIO((AUDIT/'part1_remote_tail.bin').read_bytes()))
    entries={Path(x.filename).name:x for x in index.infolist() if x.filename.endswith('.nii.gz')}
    a=pd.read_csv(AUDIT/'phase_pairing_audit.csv'); missing=sorted(a.loc[a.exclusion_reason.str.contains('CT unavailable',na=False),'source_ct'].unique())
    out=REV/'02_extraction/source_recovered/ct'; out.mkdir(parents=True,exist_ok=True)
    local=ROOT/'data_raw/WAW_TACE/ct_scans_1_4_wawtace_09_05_24.zip'
    ranges=REV/'02_extraction/part1_range_downloads'; ranges.mkdir(exist_ok=True)
    def one(name):
        info=entries[name]; offset=info.header_offset+TAIL_START; target=out/name
        result=dict(name=name,member=info.filename,offset=offset,crc_expected=info.CRC,bytes_expected=info.file_size,source='')
        if target.exists():
            b=target.read_bytes()
            if zlib.crc32(b)&0xffffffff==info.CRC: result['source']='existing_revision_crc_verified'; result['success']=True; return result
        with local.open('rb') as f:
            f.seek(offset); header=f.read(30)
            if header[:4]==b'PK\x03\x04':
                fields=struct.unpack('<4s5H3I2H',header); nl,el=fields[-2:]; f.seek(offset+30+nl+el); b=f.read(info.compress_size)
                raw=b if info.compress_type==0 else zlib.decompress(b,-15)
                if zlib.crc32(raw)&0xffffffff==info.CRC:
                    target.write_bytes(raw); result.update(source='local_bytes_recovered_using_official_central_directory',success=True); return result
        packet=ranges/(name+'.zipmember'); h=ranges/(name+'.headers.txt')
        stop=offset+30+len(info.filename.encode())+info.compress_size+1023
        cmd=['curl.exe','-L','--fail','--silent','--show-error','--retry','2','--retry-delay','20','--max-time','240','--range',f'{offset}-{stop}','--max-filesize',str(info.compress_size+4096),'--dump-header',str(h),'--output',str(packet),URL]
        cp=subprocess.run(cmd,capture_output=True,text=True)
        result['stderr']=cp.stderr
        if cp.returncode: result.update(success=False,error=f'curl exit {cp.returncode}'); return result
        payload=packet.read_bytes()
        if payload[:4]!=b'PK\x03\x04': result.update(success=False,error='Range does not start at ZIP local header'); return result
        fields=struct.unpack('<4s5H3I2H',payload[:30]); nl,el=fields[-2:]; b=payload[30+nl+el:30+nl+el+info.compress_size]
        raw=b if info.compress_type==0 else zlib.decompress(b,-15)
        if len(raw)!=info.file_size or zlib.crc32(raw)&0xffffffff!=info.CRC: result.update(success=False,error='Member length or CRC mismatch'); return result
        target.write_bytes(raw); result.update(source='official_remote_range_crc_verified',success=True,sha256=hashlib.sha256(raw).hexdigest()); return result
    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for result in pool.map(one,missing):
            records.append(result); (AUDIT/'part1_recovery_log.json').write_text(json.dumps(records,indent=2),encoding='utf-8'); print(result['name'],result['source'],result['success'],flush=True)
    print('Recovered',sum(x['success'] for x in records),'/',len(records),flush=True)

if __name__=='__main__': main()
