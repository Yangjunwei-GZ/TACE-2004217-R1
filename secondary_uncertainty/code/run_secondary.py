"""Compute matched secondary metrics from the new run actual fitted survival evidence."""
from pathlib import Path
import sys,traceback,time,collections
sys.path.insert(0,str(Path(__file__).resolve().parent))
from secondary_common import *
from secondary_metrics import *

def main():
    start=now()
    # A new computation always needs a new directory. Never erase an audit trail
    # or leave a completion certificate from an earlier attempt beside new work.
    protected_outputs=[ROOT/'results/STATISTICS_COMPLETE.json',ROOT/'results/COMPUTATION_STATUS.json',ROOT/'results/task_ledger.csv',ROOT/'logs/TASK_EVENTS.jsonl',ROOT/'logs/PROGRAM_FAILURES.jsonl']
    if any(p.exists() for p in protected_outputs):
        raise FileExistsError('REFUSED_EXISTING_OUTPUT: use run_supplement.py --output NEW_DIRECTORY; previous results and certificates were not modified')
    (ROOT/'logs').mkdir(parents=True,exist_ok=True)
    plan=readj(ROOT/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json')
    freeze=readj(ROOT/'protocol/PLAN_FREEZE_RECORD.json')
    if sha(ROOT/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json')!=freeze['plan_sha256']:raise RuntimeError('Frozen plan changed')
    m=load_matrix();grid=np.asarray(plan['ibs']['grid'],float);bound=readj(REPRO/'CURRENT_SCIENTIFIC_BINDINGS.json')
    out=ROOT/'results';(out/'reconstructed_predictions').mkdir(parents=True,exist_ok=True)
    failure_path=ROOT/'logs/PROGRAM_FAILURES.jsonl';failure_path.write_text('',encoding='utf-8')
    events_path=ROOT/'logs/TASK_EVENTS.jsonl';events_path.write_text('',encoding='utf-8')
    fold_rows=[];run_rows=[];cal_rows=[];tasks=[];failures=[];reconstruction=[];predictions=[];partitions=[]
    expected=[(mode,k) for mode,n in [('point',10),('bootstrap',500)] for k in range(n)]
    def logevent(event):
        with events_path.open('a',encoding='utf-8') as f:f.write(json.dumps(event)+'\n')
    cfg=readj(CONFIG)
    for mode,k in expected:
        p=CP/f'{mode}_{k:04d}.json.gz';fid=None
        base={'mode':mode,'run':k,'seed':cfg['seed_'+mode]+k,'checkpoint_sha256':None,'expected_checkpoint_sha256':bound['checkpoints'].get(p.name),'matrix_sha256':sha(MATRIX),'config_file_sha256':sha(CONFIG),'plan_sha256':freeze['plan_sha256']}
        logevent(dict(base,event='ATTEMPTED',time=now()))
        try:
            base['checkpoint_sha256']=sha(p)
            r=checkpoint(p)
            if sha(p)!=bound['checkpoints'][p.name]:raise ValueError('Checkpoint changed')
            if (r['mode'],r['run_id'],r['seed'])!=(mode,k,base['seed']):raise ValueError('Checkpoint run/seed identity mismatch')
            d=draw(m,r);probs={name:np.full((len(d),len(grid)),np.nan) for name in MODELS};fold_assignment=np.full(len(d),-1);supported_rows=np.zeros((len(d),len(grid)),bool)
            baselines=[];localfold=[];max_risk_err=0.;max_prob_err=0.;prob_count=0
            for f,tr,va in verified_folds(m,r):
                fid=f['fold'];rd=f['radiomics_development']
                for frame in [tr,va]:frame['radiomics_score']=rad_score(frame,rd)
                for role,frame in [('train',tr),('validation',va)]:
                    for row in frame.itertuples():partitions.append(dict(mode=mode,run=k,fold=fid,role=role,patient_id=row.patient_id,draw_row=int(frame.loc[row.Index,'_draw_row']),source_row=int(frame.loc[row.Index,'_source_row'])))
                fold_assignment[va['_draw_row'].to_numpy(int)]=fid
                for name in MODELS:
                    s=f['models'][name];coef=np.asarray(s['coefficients']);offset=s['fitted_baseline'].get('prediction_offset',0.);rtr=apply(tr,s['transform'])@coef-offset;rva=apply(va,s['transform'])@coef-offset
                    err=float(np.max(abs(rva-np.asarray(s['risk']))));max_risk_err=max(max_risk_err,err)
                    if err>1e-7:raise ValueError('Frozen risk reconstruction differs')
                    prob,baseline,supported,saved_error=saved_survival(s,rtr,tr.death_event,tr.os_months,rva,grid)
                    max_prob_err=max(max_prob_err,saved_error);prob_count+=int(supported.sum())*len(va)
                    supported_rows[va['_draw_row'].to_numpy(int)]=supported
                    probs[name][va['_draw_row'].to_numpy(int)]=prob
                    baselines.append(dict(fold=fid,model=name,**baseline))
                    for saved in s.get('time_metrics',[]):
                        if saved['time']!='IBS' and saved['status']=='evaluated':
                            j=int(saved['time'])-1;err=float(np.max(abs((1-prob[:,j])-np.asarray(saved['predicted_risk']))))
                            max_prob_err=max(max_prob_err,err);prob_count+=len(va)
                            if err>1e-7:raise ValueError('Saved held-out survival probability mismatch')
                    metrics=ipcw_metrics(tr.death_event,tr.os_months,va.death_event,va.os_months,rva,prob,grid)
                    brier=[x['Brier'] for x in metrics]
                    ibs=float(np.trapz(brier,grid)/(grid[-1]-grid[0])) if all(x is not None for x in brier) else None
                    for entry in metrics:
                        localfold.append(dict(mode=mode,run=k,seed=r['seed'],fold=fid,model=name,**entry))
                    localfold.append(dict(mode=mode,run=k,seed=r['seed'],fold=fid,model=name,time='IBS_1_36',IBS=ibs,IBS_reason='' if ibs is not None else 'incomplete_fixed_grid',grid_start=1,grid_end=36,grid_n=36,n=len(va)))
            if (fold_assignment<0).any():raise ValueError('Incomplete OOF prediction assignment')
            for name in MODELS:probability_check(probs[name])
            fold_rows.extend(localfold)
            for name in MODELS:
                for t in [12,24,36,'IBS_1_36']:
                    for metric in (['IBS'] if isinstance(t,str) else ['AUC','Brier']):
                        v=[x[metric] for x in localfold if x['model']==name and x['time']==t]
                        reasons=[x.get(metric+'_reason','') for x in localfold if x['model']==name and x['time']==t and x[metric] is None]
                        assert len(v)==5
                        run_rows.append(dict(mode=mode,run=k,seed=r['seed'],model=name,time=t,metric=metric,value=float(np.mean(v)) if all(x is not None for x in v) else None,evaluable_folds=sum(x is not None for x in v),reason=';'.join(sorted(set(reasons)))))
                for t in [12,24,36]:
                    risk=1-probs[name][:,t-1]
                    if supported_rows[:,t-1].all():cal,u,w,reason=calibration(d.patient_id,d.death_event,d.os_months,risk,t)
                    else:cal,u,w,reason=[],None,None,'out_of_training_followup_support'
                    for c in cal:cal_rows.append(dict(mode=mode,run=k,model=name,time=t,**c))
                    for metric,value in [('Calibration_MAE_unweighted',u),('Calibration_MAE_n_weighted',w)]:
                        run_rows.append(dict(mode=mode,run=k,seed=r['seed'],model=name,time=t,metric=metric,value=value,evaluable_folds=5,reason=reason))
                    for i,row in d.iterrows():predictions.append(dict(mode=mode,run=k,seed=r['seed'],model=name,time=t,patient_id=row.patient_id,draw_row=i,source_row=int(row['_source_row']),fold=int(fold_assignment[i]),death_event=int(row.death_event),os_months=float(row.os_months),survival_probability=float(probs[name][i,t-1]) if supported_rows[i,t-1] else None,probability_supported=bool(supported_rows[i,t-1])))
            for name in MODELS:probs[name][~supported_rows]=np.nan
            np.savez_compressed(out/'reconstructed_predictions'/f'{mode}_{k:04d}.npz',grid=grid,patient_id=d.patient_id.to_numpy(str),source_row=d['_source_row'].to_numpy(int),draw_row=d['_draw_row'].to_numpy(int),fold=fold_assignment,event=d.death_event.to_numpy(bool),followup=d.os_months.to_numpy(float),A=probs[A],B=probs[B],supported_by_training_followup=supported_rows)
            with gzip.open(out/'reconstructed_predictions'/f'{mode}_{k:04d}_baseline.json.gz','wt',encoding='utf-8') as f:json.dump(baselines,f)
            rec=dict(base,status='SUCCESS',folds=5,models=2,max_linear_risk_difference=max_risk_err,max_saved_survival_difference=max_prob_err,saved_probabilities_checked=prob_count,coefficient_refits=0,fitted_baselines_reused=10)
            tasks.append(rec);reconstruction.append(rec);logevent(dict(rec,event='SUCCESS',time=now()))
        except Exception as exc:
            record=dict(base,fold=fid,status='FAILED',exception_type=type(exc).__name__,exception=str(exc),traceback=traceback.format_exc(),time=now())
            with failure_path.open('a',encoding='utf-8') as f:f.write(json.dumps(record)+'\n')
            failures.append(record);tasks.append(dict(base,status='FAILED'));logevent(dict(record,event='FAILED'))
        if k%50==0:print(mode,k,'successful',sum(x['status']=='SUCCESS' for x in tasks),'failures',len(failures),flush=True)
    pd.DataFrame(tasks).to_csv(out/'task_ledger.csv',index=False)
    pd.DataFrame(fold_rows).to_csv(out/'fold_metrics.csv.gz',index=False,compression='gzip')
    pd.DataFrame(cal_rows).to_csv(out/'calibration_group_details.csv.gz',index=False,compression='gzip')
    pd.DataFrame(run_rows).to_csv(out/'run_metrics.csv',index=False)
    pd.DataFrame(predictions).to_csv(out/'patient_OOF_fixed_time_predictions.csv.gz',index=False,compression='gzip')
    pd.DataFrame(partitions).to_csv(out/'original_train_validation_partitions.csv.gz',index=False,compression='gzip')
    pd.DataFrame(reconstruction).to_csv(ROOT/'audit/reconstruction_agreement.csv',index=False)
    status={'status':'COMPUTED_PENDING_INDEPENDENT_QA' if not failures else 'FAILED','started_at':start,'finished_at':now(),'planned':len(expected),'attempted':len(tasks),'successful':sum(x['status']=='SUCCESS' for x in tasks),'failed':len(failures),'coefficient_refits':0,'feature_selection_refits':0,'fitted_baselines_reused':sum(x.get('fitted_baselines_reused',0) for x in tasks),'independent_fixed_beta_baseline_checks':sum(x.get('fitted_baselines_reused',0) for x in tasks)}
    writej(out/'COMPUTATION_STATUS.json',status)
    completeness_gate(expected,tasks,failures)
    print(json.dumps(status,indent=2))
if __name__=='__main__':
    raise SystemExit('Use run_new_secondary.py compute with explicit new revision/runtime/output; direct engine execution is disabled.')
