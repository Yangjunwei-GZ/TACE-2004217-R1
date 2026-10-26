"""Parameterised, fail-closed adapter for the newly frozen TACE cohort.

Commands: freeze (before secondary metric computation), compute, audit.
No model is trained here. Only this revision's actual fitted checkpoints are read.
"""
from pathlib import Path
import argparse, collections, hashlib, importlib, json, os, platform, sys, traceback
from datetime import datetime
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_key,'1')

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))  # Required by the locked embeddable CPython ._pth isolation.
A='clinical_diameter_mean_hu'; B=A+'_radiomics'
EXPECTED=[(mode,k) for mode,n in [('point',10),('bootstrap',500)] for k in range(n)]

def now():return datetime.now().astimezone().isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def readj(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def writej(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
def code_hash(directory):
    hashes={p.name:sha(p) for p in sorted(Path(directory).glob('*.py'))}
    return hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
def paths(args):
    r=Path(args.revision).resolve();o=Path(args.output).resolve();rt=Path(args.runtime).resolve()
    if o==r or r in o.parents:raise ValueError('Secondary output must be outside the protected primary revision tree')
    return r,o,rt,r/'03_model_development/analysis_matrix.csv',r/'00_protocol_and_provenance/revision_config.json'
def require_locked_config(cfg):
    expected={'point_repeats':10,'bootstrap_replicates':500,'outer_folds':5,'inner_folds':3,
        'seed_point':20261014,'seed_bootstrap':20261015,'primary_comparison':[A,B],
        'prediction_months':[12,24,36],'ibs_monthly_grid':[1,36], 'baseline_penalizer':.1}
    for k,v in expected.items():
        if cfg.get(k)!=v:raise ValueError('Frozen protocol mismatch: '+k)
def freeze(args):
    import pandas as pd
    r,o,rt,matrix,cfgpath=paths(args);cfg=readj(cfgpath);require_locked_config(cfg)
    if (o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json').exists():raise FileExistsError('Plan already exists; use another new output directory')
    m=pd.read_csv(matrix,dtype={'patient_id':str})
    if m.patient_id.duplicated().any() or m.patient_id.isna().any():raise ValueError('Patient matrix must have one unique nonmissing patient identifier per row')
    p=readj(HERE/'ORIGINAL_SECONDARY_METHODS_REFERENCE.json')
    p.update(frozen_at=now(),stage='Revision-stage supplementary secondary evaluation applied to the newly corrected cohort; not original-study prespecification',
        parent_method_plan_sha256=sha(HERE/'ORIGINAL_SECONDARY_METHODS_REFERENCE.json'),
        matrix_sha256=sha(matrix),config_sha256=sha(cfgpath),configuration_sha256=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest(),
        runtime_source_sha256=code_hash(rt),secondary_code_sha256=code_hash(HERE),
        actual_patient_count=len(m),actual_events=int(m.death_event.sum()),patient_ids=m.patient_id.tolist(),
        source={'revision':str(r),'runtime':str(rt),'matrix':str(matrix),'config':str(cfgpath)},
        model_rules='Use this new run saved outer-training transforms, selected radiomics coefficients including Coxnet offset, fitted Cox/KM baseline and A/B monthly survival evidence. Independently check fitted Breslow against exact training rows at fixed coefficients. No selection, tuning or regression coefficient fitting in this secondary stage. Primary models are newly trained by the separate scientific execution.',
        protected='Do not change current or historical CT, masks, clinical decisions, patient/ROI lists, matrices, regression fits, random seeds, primary metrics or checkpoint bytes. Do not import prior 118-patient secondary results. Published v1.0.0 remains unchanged.',
        verification='Validate complete current primary task audit and all checkpoint identities, draws, outer/inner partitions, training outcomes, model specifications, source/config/matrix SHA; check saved fitted probability against a separate fixed-beta reconstruction; recompute all secondary paired intervals; verify source hashes again after computation.')
    p.pop('source_audit_sha256',None)
    p['ibs']['support_basis']='The same fixed monthly 1..36-month window is retained without inspecting model performance. New-run follow-up and censoring support are audited before metrics. Unsupported folds/replicates are retained as NA; the window is not adjusted to increase evaluability or significance.'
    p['sensitivity']='Ten current point repeats are reported separately as split-sensitivity mean and SD. Single-repeat bootstrap intervals are not attached to the ten-repeat mean. Prior cohort results and variable-window IBS are not reused.'
    p['bootstrap']='Current 500 full-development patient-bootstrap checkpoints, original seeds 20261015..20261514 and grouped folds. Patient multiplicities are retained; duplicates cannot cross training and validation. Percentile 2.5/97.5 linear quantiles of matched single-repeat statistics; same common evaluable set for A, B and B-A. Internal randomized-CV uncertainty, not external performance of the final formula.'
    for sub in ('protocol','audit','logs','results'): (o/sub).mkdir(parents=True,exist_ok=True)
    writej(o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json',p)
    writej(o/'protocol/PLAN_FREEZE_RECORD.json',{'frozen_at':now(),'plan_sha256':sha(o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json'),'command':sys.argv,'actual_environment':environment()})
    print(json.dumps({'status':'SECONDARY_METHODS_FROZEN','patients':len(m),'events':int(m.death_event.sum()),'plan':str(o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json')},ensure_ascii=False))

def environment():
    from importlib.metadata import version
    v={}
    for package in ('numpy','pandas','scipy','scikit-learn','scikit-survival','lifelines'):
        try:v[package]=version(package)
        except Exception:v[package]='not installed'
    return dict(python=sys.version,executable=sys.executable,platform=platform.platform(),packages=v)

def current_plan(args):
    r,o,rt,matrix,cfgpath=paths(args);p=readj(o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json');f=readj(o/'protocol/PLAN_FREEZE_RECORD.json')
    if sha(o/'protocol/SUPPLEMENTARY_ANALYSIS_PLAN.json')!=f['plan_sha256']:raise ValueError('Supplementary frozen plan changed')
    for actual,expected,label in [(sha(matrix),p['matrix_sha256'],'matrix'),(sha(cfgpath),p['config_sha256'],'configuration'),(code_hash(rt),p['runtime_source_sha256'],'runtime'),(code_hash(HERE),p['secondary_code_sha256'],'secondary code')]:
        if actual!=expected:raise ValueError('Frozen '+label+' fingerprint changed')
    if str(r)!=p['source']['revision'] or str(rt)!=p['source']['runtime']:raise ValueError('Run paths differ from frozen supplemental bindings')
    require_locked_config(readj(cfgpath));return p

def validate_current_inputs(args):
    """Read all existing files without invoking a primary audit that writes them."""
    import numpy as np
    import pandas as pd
    from sksurv.nonparametric import CensoringDistributionEstimator
    from sksurv.util import Surv
    from secondary_common import checkpoint,draw
    r,o,rt,matrix,cfgpath=paths(args);p=current_plan(args)
    integrity_path=r/'04_validation/VALIDATION_INTEGRITY.json';integrity=readj(integrity_path)
    if integrity.get('status')!='COMPLETE':raise RuntimeError('Current primary validation is not COMPLETE')
    vp_path=r/'04_validation/VALIDATION_PLAN.json';vp=readj(vp_path)
    inputs={'analysis_matrix_sha256':p['matrix_sha256'],'configuration_file_sha256':p['config_sha256'],'configuration_sha256':p['configuration_sha256']}
    if vp.get('inputs')!=inputs or integrity.get('input_fingerprints')!=inputs:raise ValueError('Primary audit and plan do not bind the current matrix/config')
    if vp.get('execution_source_sha256')!=p['runtime_source_sha256']:raise ValueError('Primary runtime source differs from frozen secondary input binding')
    if vp.get('matrix',{}).get('patient_ids')!=p['patient_ids']:raise ValueError('Primary plan patient identities differ')
    sys.path.insert(0,str(rt));vi=importlib.import_module('validation_integrity')
    directory=r/'04_validation/checkpoints';expected={f'{mode}_{k:04d}.json.gz' for mode,k in EXPECTED}
    observed={x.name for x in directory.glob('*.json.gz')}
    if observed!=expected:raise ValueError('Checkpoint set missing/unexpected: '+json.dumps({'missing':sorted(expected-observed),'unexpected':sorted(observed-expected)}))
    manifest={};support=[];m=pd.read_csv(matrix,dtype={'patient_id':str})
    for mode,k in EXPECTED:
        fp=directory/f'{mode}_{k:04d}.json.gz';cp=checkpoint(fp)
        issues=vi.checkpoint_issues(cp,vp,fp)
        if issues:raise ValueError('Current checkpoint failed full protocol audit: '+fp.name+' '+json.dumps(issues,ensure_ascii=False))
        manifest[fp.name]=sha(fp);d=draw(m,cp)
        for fold in cp['folds']:
            tr=d.iloc[fold['train_draw_row_indices']];va=d.iloc[fold['validation_draw_row_indices']]
            censor=CensoringDistributionEstimator().fit(Surv.from_arrays(tr.death_event.to_numpy(bool),tr.os_months.to_numpy(float)))
            for t in (12,24,36):
                try:g=float(censor.predict_proba([t])[0]);reason='' if np.isfinite(g) and g>0 else 'nonpositive_censoring_survival'
                except ValueError:g=None;reason='censoring_support_insufficient'
                support.append(dict(mode=mode,run=k,seed=cp['seed'],fold=fold['fold'],time=t,train_max=float(tr.os_months.max()),validation_max=float(va.os_months.max()),followup_supported=bool(t<min(tr.os_months.max(),va.os_months.max())),censor_survival_t=g,censoring_reason=reason))
    bound={'matrix_sha256':p['matrix_sha256'],'config_sha256':p['config_sha256'],'checkpoints':manifest,'primary_integrity_sha256':sha(integrity_path),'primary_plan_sha256':sha(vp_path),'runtime_source_sha256':code_hash(rt),'recorded_at':now()}
    writej(o/'audit/source_binding/CURRENT_SCIENTIFIC_BINDINGS.json',bound)
    pd.DataFrame(support).to_csv(o/'audit/NEW_RUN_FOLLOWUP_SUPPORT.csv',index=False)
    writej(o/'audit/NEW_RUN_SOURCE_AUDIT.json',dict(status='COMPLETE_CURRENT_PRIMARY_INPUTS_VERIFIED',checked_at=now(),patients=len(m),events=int(m.death_event.sum()),checkpoints=len(manifest),folds=sum(1 for _ in support)//3,source=bound,actual_environment=environment()))
    return bound

def audit(args):
    """Independent summary recomputation and source-byte preservation gate."""
    import numpy as np
    import pandas as pd
    from secondary_metrics import completeness_gate
    r,o,rt,matrix,cfgpath=paths(args);p=current_plan(args)
    out=o/'results';tasks=pd.read_csv(out/'task_ledger.csv').to_dict('records')
    failures=[json.loads(x) for x in (o/'logs/PROGRAM_FAILURES.jsonl').read_text().splitlines() if x]
    completeness_gate(EXPECTED,tasks,failures)
    events=[json.loads(x) for x in (o/'logs/TASK_EVENTS.jsonl').read_text().splitlines() if x]
    attempts=[(e['mode'],e['run']) for e in events if e['event']=='ATTEMPTED']
    successes=[(e['mode'],e['run']) for e in events if e['event']=='SUCCESS']
    if collections.Counter(attempts)!=collections.Counter(EXPECTED) or collections.Counter(successes)!=collections.Counter(EXPECTED):raise ValueError('Planned/attempted/succeeded event sets differ')
    df=pd.read_csv(out/'run_metrics.csv',dtype={'time':str});summary=pd.read_csv(out/'uncertainty_summary.csv',dtype={'time':str});paired=pd.read_csv(out/'paired_bootstrap_metrics.csv',dtype={'time':str})
    if len(df)!=510*26 or len(summary)!=13 or len(paired)!=500*13:raise ValueError('Metric row coverage differs from the frozen 13 endpoints x two models')
    def same(x,y):
        if pd.isna(x) and pd.isna(y):return
        if not np.isclose(float(x),float(y),rtol=1e-11,atol=1e-12):raise ValueError('Independent numeric recomputation mismatch: '+str((x,y)))
    for row in summary.to_dict('records'):
        g=df[(df.metric==row['metric'])&(df.time==row['time'])]
        if g.duplicated(['mode','run','model']).any():raise ValueError('Duplicate per-run endpoint')
        w=g.pivot(index=['mode','run'],columns='model',values='value')
        if set(w.index)!=set(EXPECTED):raise ValueError('Missing per-run endpoint')
        b=w.loc['bootstrap'];valid=b[A].notna()&b[B].notna();n=int(valid.sum());pt=w.loc[('point',0)]
        if row['bootstrap_joint_evaluable']!=n or row['bootstrap_NA']!=500-n:raise ValueError('Non-evaluable denominator mismatch')
        good=n>=450 and pt[[A,B]].notna().all()
        if row['CI_status']!=('estimated' if good else 'not_reliably_estimable'):raise ValueError('CI estimability threshold mismatch')
        for name,key in ((A,'A'),(B,'B'),('delta','delta')):
            val=pt[B]-pt[A] if name=='delta' else pt[name]
            same(row['delta_B_minus_A' if name=='delta' else key+'_point'],val)
            z=(b.loc[valid,B]-b.loc[valid,A]).to_numpy() if name=='delta' else b.loc[valid,name].to_numpy()
            ci=np.quantile(z,[.025,.975],method='linear') if good else [np.nan,np.nan]
            same(row[key+'_CI_lower'],ci[0]);same(row[key+'_CI_upper'],ci[1])
        pp=paired[(paired.metric==row['metric'])&(paired.time==row['time'])].sort_values('run')
        if pp.run.tolist()!=list(range(500)) or pp.seed.tolist()!=list(range(20261015,20261515)):raise ValueError('Paired rows lost original run/seed')
        if not np.allclose(pp.A,b.sort_index()[A],equal_nan=True) or not np.allclose(pp.B,b.sort_index()[B],equal_nan=True) or not np.allclose(pp.delta,pp.B-pp.A,equal_nan=True):raise ValueError('Paired table differs from per-run estimates')
    na=pd.read_csv(out/'NON_EVALUABLE_RUNS.csv',dtype={'time':str})
    if len(na)!=int(df.value.isna().sum()) or na.reason.isna().any():raise ValueError('NA table and reasons do not reconcile')
    folds=pd.read_csv(out/'fold_metrics.csv.gz',dtype={'time':str})
    for (mode,k,fold,model),f in folds.groupby(['mode','run','fold','model'],sort=False):
        ib=f[f.time=='IBS_1_36'].iloc[0];monthly=f[f.time!='IBS_1_36'].copy();monthly['t']=monthly.time.astype(float);monthly=monthly.sort_values('t')
        if monthly.t.tolist()!=list(map(float,range(1,37))) or len(f)!=37:raise ValueError('IBS monthly grid changed/missing')
        expected_ibs=float(np.trapz(monthly.Brier,monthly.t)/35) if monthly.Brier.notna().all() else np.nan
        same(ib.IBS,expected_ibs)
    bound=readj(o/'audit/source_binding/CURRENT_SCIENTIFIC_BINDINGS.json')
    for name,value in bound['checkpoints'].items():
        if sha(r/'04_validation/checkpoints'/name)!=value:raise ValueError('Source checkpoint changed after computation: '+name)
    if sha(r/'04_validation/VALIDATION_INTEGRITY.json')!=bound['primary_integrity_sha256'] or sha(r/'04_validation/VALIDATION_PLAN.json')!=bound['primary_plan_sha256']:raise ValueError('Source validation certificate/plan changed')
    result={'status':'STATISTICS_COMPLETE_AUTOMATED_QA','verified_at':now(),'patients':p['actual_patient_count'],'events':p['actual_events'],'planned':510,'attempted':len(attempts),'successful':len(successes),'failed':0,'coefficient_refits_in_secondary_stage':0,'regression_refits_in_secondary_stage':0,'fitted_baselines_reused':5100,'interval_rows':len(summary),'intervals_estimable':int((summary.CI_status=='estimated').sum()),'not_reliably_estimable':int((summary.CI_status!='estimated').sum()),'NA_model_endpoint_rows':len(na),'source_hashes_unchanged':True,'note':'Automated numerical and integrity audit; report limitations and every NA denominator. This certificate does not assert external validation or scientific benefit.'}
    writej(o/'audit/SECONDARY_INDEPENDENT_RECOMPUTATION.json',result)
    writej(out/'STATISTICS_COMPLETE.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2));return result

def configure(args):
    r,o,rt,_,_=paths(args)
    os.environ['TACE_SECONDARY_OUTPUT']=str(o);os.environ['TACE_NEW_REVISION']=str(r);os.environ['TACE_NEW_RUNTIME']=str(rt)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['freeze','compute','audit'])
    parser.add_argument('--revision',required=True,help='New revision folder containing 00_protocol_and_provenance,03_model_development,04_validation')
    parser.add_argument('--runtime',required=True,help='Actual Python model runtime whose source SHA appears in new checkpoints')
    parser.add_argument('--output',required=True,help='Separate new secondary output directory, outside protected primary revision tree')
    args=parser.parse_args();configure(args)
    try:
        if args.command=='freeze':freeze(args)
        elif args.command=='audit':audit(args)
        else:
            o=Path(args.output)
            if any((o/x).exists() for x in ['results/task_ledger.csv','results/COMPUTATION_STATUS.json','results/STATISTICS_COMPLETE.json','logs/TASK_EVENTS.jsonl']):raise FileExistsError('Existing computation/audit preserved; create and freeze a new output directory')
            validate_current_inputs(args)
            import run_secondary,summarize_secondary
            run_secondary.main();summarize_secondary.main();audit(args)
    except Exception as exc:
        o=Path(args.output);entry={'time':now(),'command':sys.argv,'status':'REFUSED_EXISTING_OUTPUT' if isinstance(exc,FileExistsError) else 'FAILED','exception_type':type(exc).__name__,'exception':str(exc),'traceback':traceback.format_exc(),'previous_outputs_preserved':True}
        (o/'logs').mkdir(parents=True,exist_ok=True)
        with (o/'logs/ENTRY_FAILURES.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(entry,ensure_ascii=False)+'\n')
        print(traceback.format_exc(),file=sys.stderr);return 1
    return 0
if __name__=='__main__':sys.exit(main())
