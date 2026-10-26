"""New-run bindings only. No old patient matrix, result or checkpoint is imported."""
from pathlib import Path
import os, sys, json, gzip, hashlib
from datetime import datetime
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sksurv.util import Surv

A='clinical_diameter_mean_hu'; B=A+'_radiomics'; MODELS=[A,B]
ROOT=Path(os.environ.get('TACE_SECONDARY_OUTPUT', str(Path(__file__).resolve().parent.parent/'audit/secondary_not_configured')))
SOURCE=Path(os.environ.get('TACE_NEW_REVISION','__EXPLICIT_NEW_REVISION_REQUIRED__'))
REPRO=ROOT/'audit/source_binding'
CP=SOURCE/'04_validation/checkpoints'
MATRIX=SOURCE/'03_model_development/analysis_matrix.csv'
CONFIG=SOURCE/'00_protocol_and_provenance/revision_config.json'
RUNTIME=Path(os.environ.get('TACE_NEW_RUNTIME','__EXPLICIT_NEW_RUNTIME_REQUIRED__'))

def now(): return datetime.now().astimezone().isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def readj(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def writej(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def load_matrix():return pd.read_csv(MATRIX,dtype={'patient_id':str,'sample_id':str})
def labels(d):return Surv.from_arrays(d.death_event.to_numpy(bool),d.os_months.to_numpy(float))
def apply(d,tr):
    x=d[tr['columns']].to_numpy(float)
    x=np.where(np.isnan(x),np.asarray(tr['median']),x)
    return ((x-np.asarray(tr['mean']))/np.asarray(tr['sd']))[:,np.asarray(tr['keep'],bool)]
def checkpoint(path):
    with gzip.open(path,'rt',encoding='utf-8') as handle:return json.load(handle)
def draw(matrix,r):
    idx=np.arange(len(matrix)) if r['mode']=='point' else np.random.default_rng(r['seed']).integers(0,len(matrix),len(matrix))
    if r.get('draw_source_row_indices')!=idx.tolist():raise ValueError('Stored draw indices disagree with original seed')
    d=matrix.iloc[idx].copy().reset_index(drop=True)
    if r.get('draw_patient_ids')!=d.patient_id.tolist():raise ValueError('Draw patient identity mismatch')
    d['_source_row']=idx;d['_draw_row']=np.arange(len(d))
    return d
def verified_folds(matrix,r):
    d=draw(matrix,r); ids=np.asarray(sorted(d.patient_id.unique()))
    cfg=readj(CONFIG);splits=list(KFold(cfg['outer_folds'],shuffle=True,random_state=r['seed']).split(ids))
    if len(r['folds'])!=len(splits):raise ValueError('Missing or extra outer folds')
    for k,((it,iv),f) in enumerate(zip(splits,r['folds'])):
        tr=d[d.patient_id.isin(ids[it])].copy();va=d[d.patient_id.isin(ids[iv])].copy()
        if f['fold']!=k or f['train_ids']!=tr.patient_id.tolist() or f['validation_ids']!=va.patient_id.tolist():raise ValueError('Original seeded fold identity mismatch')
        if set(tr.patient_id)&set(va.patient_id):raise ValueError('Patient leakage across training and validation')
        if f.get('train_draw_row_indices')!=tr['_draw_row'].tolist() or f.get('validation_draw_row_indices')!=va['_draw_row'].tolist():raise ValueError('Draw-row partition mismatch')
        if f.get('training_events')!=tr.death_event.tolist() or f.get('training_times')!=tr.os_months.tolist():raise ValueError('Training outcomes mismatch')
        yield f,tr,va
def binding_paths():return readj(REPRO/'CURRENT_SCIENTIFIC_BINDINGS.json')

def rad_score(d,development):
    # Coxnet.predict subtracts an offset. It is not generally X beta.
    offset=development['fitted_baseline'].get('prediction_offset',0.)
    return apply(d,development['transform'])@np.asarray(development['coefficients'])-offset

def saved_survival(s,train_risk,event,time,val_risk,grid):
    """Use the actually fitted baseline and check its serialized probabilities.

    Fixed-beta Breslow reconstruction is an independent numerical check, not a
    regression refit. Unsupported times remain NA in exported predictions.
    A finite mathematical extension is used only inside a metric routine that
    rejects those times before computing a score; it is never reported as a
    supported prediction.
    """
    from secondary_metrics import reconstruct_survival,probability_check
    baseline=s.get('fitted_baseline'); monthly=s.get('monthly_survival')
    if not isinstance(baseline,dict) or not isinstance(monthly,dict):
        raise ValueError('New checkpoint lacks fitted baseline/monthly probability evidence; do not reuse an old checkpoint')
    x=np.asarray(baseline['times'],float); idx=np.searchsorted(x,grid,side='right')-1
    keep=idx>=0
    if baseline['kind']=='fitted_Breslow':
        h=np.zeros(len(grid));h[keep]=np.asarray(baseline['cumulative_hazard'])[idx[keep]]
        prob=np.exp(-np.outer(np.exp(val_risk),h))
        independently_reconstructed,_=reconstruct_survival(train_risk,event,time,val_risk,grid)
        max_error=float(np.max(abs(prob-independently_reconstructed)))
        if max_error>1e-7:raise ValueError('Saved baseline disagrees with fixed-beta Breslow reconstruction')
    elif baseline['kind']=='training_KM_null_model':
        surv=np.ones(len(grid));surv[keep]=np.asarray(baseline['baseline_survival'])[idx[keep]]
        prob=np.tile(surv,(len(val_risk),1));max_error=0.
        from sksurv.nonparametric import kaplan_meier_estimator
        tx,sy=kaplan_meier_estimator(np.asarray(event,bool),np.asarray(time,float))
        ji=np.searchsorted(tx,grid,side='right')-1;vv=np.ones(len(grid));valid=ji>=0;vv[valid]=sy[ji[valid]]
        if not np.allclose(prob,np.tile(vv,(len(val_risk),1)),atol=1e-12):raise ValueError('Null-model baseline differs from training KM')
    else:raise ValueError('Unknown saved baseline kind')
    probability_check(prob)
    if float(baseline['max_training_followup'])!=float(np.max(time)):raise ValueError('Baseline support disagrees with training times')
    supported=np.asarray(grid)<=np.max(time)
    if monthly['months']!=list(map(int,grid)) or monthly['supported_by_training_followup']!=supported.tolist():raise ValueError('Saved monthly grid/support mismatch')
    stored=np.asarray([[np.nan if x is None else x for x in row] for row in monthly['survival_probability']],float)
    if stored.shape!=prob.shape or np.any(np.isfinite(stored[:,~supported])) or not np.isfinite(stored[:,supported]).all():raise ValueError('Saved monthly probability shape/NA mismatch')
    saved_error=float(np.max(abs(stored[:,supported]-prob[:,supported]))) if supported.any() else 0.
    if saved_error>1e-12:raise ValueError('Saved monthly probability differs from fitted baseline')
    return prob,dict(**baseline,independent_fixed_beta_error=max_error,saved_probability_error=saved_error),supported,saved_error
