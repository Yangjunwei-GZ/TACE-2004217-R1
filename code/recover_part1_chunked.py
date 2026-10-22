"""Bounded resumable member-range recovery when single Zenodo streams are throttled."""
from pathlib import Path
import io,json,struct,zipfile,zlib,subprocess,concurrent.futures,hashlib,time,threading,argparse
import pandas as pd
ROOT=Path(str(Path.cwd())); REV=ROOT/'revision_R1_20261014'; AUDIT=REV/'01_data_audit'
URL='https://zenodo.org/records/12741586/files/ct_scans_1_4_wawtace_09_05_24.zip?download=1'; START=10358430958

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--only'); args=parser.parse_args()
    z=zipfile.ZipFile(io.BytesIO((AUDIT/'part1_remote_tail.bin').read_bytes())); entries={Path(i.filename).name:i for i in z.infolist() if i.filename.endswith('.nii.gz')}
    a=pd.read_csv(AUDIT/'phase_pairing_audit.csv'); missing=sorted(a.loc[a.exclusion_reason.str.contains('CT unavailable',na=False),'source_ct'].unique(),key=lambda n:('_1_scan' not in n,n))
    if args.only: missing=[n for n in missing if n==args.only]
    prefix='part1_member_retry' if args.only else 'part1'
    folder=REV/'02_extraction/part1_chunks'; folder.mkdir(exist_ok=True); output=REV/'02_extraction/source_recovered/ct'; output.mkdir(parents=True,exist_ok=True)
    jobs=[]; plan={}; chunk=524288
    for name in missing:
        info=entries[name]; offset=info.header_offset+START; length=30+len(info.filename.encode())+info.compress_size
        parts=[]
        for j,start in enumerate(range(offset,offset+length,chunk)):
            end=min(start+chunk,offset+length)-1; p=folder/f'{name}.{j:04d}.part'; parts.append((p,end-start+1)); jobs.append((p,start,end))
        plan[name]=(info,parts)
    throttle=threading.Lock(); state={'next':0.0,'cooldown':0.0}
    def permit():
        with throttle:
            now=time.monotonic(); wait=max(state['next'],state['cooldown'])-now
            if wait>0: time.sleep(wait)
            state['next']=time.monotonic()+.65
    def fetch(job):
        p,start,end=job; expected=end-start+1
        if p.exists() and p.stat().st_size==expected: return dict(file=p.name,status='resumed_complete',bytes=expected)
        cmd=['curl.exe','-L','--fail','--silent','--show-error','--connect-timeout','25','--max-time','120','--range',f'{start}-{end}','--max-filesize',str(expected),'--output',str(p),URL]
        for attempt in range(3):
            permit(); c=subprocess.run(cmd,capture_output=True,text=True)
            if c.returncode==0: break
            if '429' in c.stderr:
                with throttle: state['cooldown']=max(state['cooldown'],time.monotonic()+65)
            elif attempt<2: time.sleep(5)
        return dict(file=p.name,status='complete' if c.returncode==0 and p.exists() and p.stat().st_size==expected else 'failed',bytes=p.stat().st_size if p.exists() else 0,error=c.stderr[-600:])
    records=[]; t=time.time(); success=0; recovered=[]
    remaining={name:len(parts) for name,(info,parts) in plan.items()}
    def assemble(name):
        info,parts=plan[name]
        if not all(p.exists() and p.stat().st_size==n for p,n in parts): return dict(name=name,success=False,reason='incomplete chunks')
        b=b''.join(p.read_bytes() for p,n in parts)
        if b[:4]!=b'PK\x03\x04': return dict(name=name,success=False,reason='invalid local header')
        fields=struct.unpack('<4s5H3I2H',b[:30]); nl,el=fields[-2:]; payload=b[30+nl+el:30+nl+el+info.compress_size]; raw=payload if info.compress_type==0 else zlib.decompress(payload,-15)
        if len(raw)!=info.file_size or zlib.crc32(raw)&0xffffffff!=info.CRC: return dict(name=name,success=False,reason='CRC/length mismatch')
        target=output/name
        if not target.exists(): target.write_bytes(raw)
        return dict(name=name,success=True,sha256=hashlib.sha256(raw).hexdigest(),crc32=info.CRC,bytes=len(raw))
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as pool:
        fs=[pool.submit(fetch,j) for j in jobs]
        for f in concurrent.futures.as_completed(fs):
            r=f.result(); records.append(r); success+=r['status']!='failed'
            name=r['file'].rsplit('.',2)[0]; remaining[name]-=1
            if remaining[name]==0:
                verified=assemble(name); recovered.append(verified)
                (AUDIT/f'{prefix}_final_recovery.json').write_text(json.dumps(recovered,indent=2),encoding='utf-8')
                print('MEMBER',name,verified['success'],flush=True)
            if len(records)%20==0:
                (AUDIT/f'{prefix}_chunk_recovery_checkpoint.json').write_text(json.dumps(dict(completed=len(records),successful=success,total=len(jobs),elapsed=time.time()-t,last=records[-20:]),indent=2),encoding='utf-8'); print(f'CHUNKS {len(records)}/{len(jobs)} successful {success}',flush=True)
    (AUDIT/f'{prefix}_chunk_recovery_records.json').write_text(json.dumps(records,indent=2),encoding='utf-8'); (AUDIT/f'{prefix}_final_recovery.json').write_text(json.dumps(recovered,indent=2),encoding='utf-8'); print('MEMBERS recovered',sum(x['success'] for x in recovered),'/',len(recovered),flush=True)

if __name__=='__main__':main()
