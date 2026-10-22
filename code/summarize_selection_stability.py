"""Describe saved training-fit selection; never train, select a new rule or rewrite checkpoints."""
from pathlib import Path
import argparse,csv,gzip,json,hashlib
TARGETS=['original_glszm_SizeZoneNonUniformity','original_glszm_LargeAreaHighGrayLevelEmphasis']
def summarize(validation,output):
    validation=Path(validation);output=Path(output)
    if output.exists():raise ValueError('Use a new output directory')
    with (validation/'feature_recurrence_records.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    logs={(r['mode'],int(r['run']),int(r['fold']),r['feature']):r for r in rows}
    if len(logs)!=len(rows):raise ValueError('Duplicate feature log identities')
    seen=set();fits=[];hashes={}
    for p in sorted((validation/'checkpoints').glob('*.json.gz')):
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        with gzip.open(p,'rt',encoding='utf-8') as f:r=json.load(f)
        for fold in r['folds']:
            selected=[];dev=fold['radiomics_development']
            for name,beta in zip(dev['features'],dev['coefficients']):
                k=(r['mode'],r['run_id'],fold['fold'],name);row=logs[k];flag=abs(beta)>1e-8
                if abs(float(row['coefficient'])-beta)>1e-14 or (row['selected'].lower()=='true')!=flag:raise ValueError('Checkpoint/log disagreement')
                if k in seen:raise ValueError('Duplicate checkpoint fit identity')
                seen.add(k)
                if flag:selected.append(name)
            fits.append({'mode':r['mode'],'run':r['run_id'],'fold':fold['fold'],**{t:t in selected for t in TARGETS},'both_selected':all(t in selected for t in TARGETS),'selected_feature_count':len(selected)})
    if seen!=set(logs):raise ValueError('Unmatched checkpoint/log records')
    summary=[]
    for mode in ['point','bootstrap']:
        f=[r for r in fits if r['mode']==mode]
        if not f:raise ValueError('Missing validation mode')
        for t in TARGETS+['both_selected']:
            count=sum(r[t] for r in f);summary.append({'mode':mode,'feature':t,'selected_fits':count,'total_fits':len(f),'percent':count/len(f)*100})
    output.mkdir(parents=True)
    for name,data in [('fit_records.csv',fits),('selection_summary.csv',summary)]:
        with (output/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    result={'definition':'abs(fitted radiomics coefficient)>1e-8','new_model_fits':0,'fit_count':len(fits),'checkpoints':hashes,'summary':summary,'interpretation':'Dependent-fit descriptive recurrence; joint selection is not identical complete-model reproduction. No binomial interval.'}
    (output/'selection_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--validation',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=summarize(a.validation,a.output);print(json.dumps({'status':'PASS_CHECKPOINT_LOG_MATCH','fits':r['fit_count'],'checkpoints':len(r['checkpoints']),'new_model_fits':0}))
