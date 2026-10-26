"""Frozen secondary estimators. All censoring estimates are training-only."""
import numpy as np
import pandas as pd
from scipy.stats import norm
from sksurv.nonparametric import CensoringDistributionEstimator
from sksurv.linear_model.coxph import BreslowEstimator
from sksurv.util import Surv

class NotEvaluable(Exception):pass

def probability_check(prob):
    prob=np.asarray(prob,float)
    if not np.isfinite(prob).all():raise ValueError('Nonfinite survival probability')
    if np.min(prob)<-1e-12 or np.max(prob)>1+1e-12:raise ValueError('Probability outside [0,1]')
    if prob.ndim==2 and np.any(np.diff(prob,axis=1)>1e-12):raise ValueError('Survival increases with time')

def reconstruct_survival(train_risk,event,time,val_risk,grid):
    base=BreslowEstimator().fit(np.asarray(train_risk),np.asarray(event,bool),np.asarray(time,float))
    x=base.cum_baseline_hazard_.x;y=base.cum_baseline_hazard_.y
    idx=np.searchsorted(x,grid,side='right')-1
    H=np.where(idx<0,0.,y[np.maximum(idx,0)])
    prob=np.exp(-np.outer(np.exp(val_risk),H))
    probability_check(prob)
    return prob,{'time':x.tolist(),'cumulative_hazard':y.tolist(),'grid_cumulative_hazard':H.tolist()}

def ipcw_metrics(train_event,train_time,val_event,val_time,risk,prob,grid):
    """Explicit IPCW formulas; do not require t >= min validation follow-up."""
    e=np.asarray(val_event,bool);y=np.asarray(val_time,float);risk=np.asarray(risk,float)
    grid=np.asarray(grid,float);prob=np.asarray(prob,float);probability_check(prob)
    censor=CensoringDistributionEstimator().fit(Surv.from_arrays(np.asarray(train_event,bool),np.asarray(train_time,float)))
    rows=[]
    for j,t in enumerate(grid):
        row={'time':float(t),'AUC':None,'Brier':None,'AUC_reason':'','Brier_reason':''}
        case=(y<=t)&e;control=y>t
        row.update(n=len(y),cases=int(case.sum()),controls=int(control.sum()),censored_before=int(((y<=t)&~e).sum()))
        if t>=min(np.max(train_time),y.max()):
            row['AUC_reason']=row['Brier_reason']='followup_support_insufficient';rows.append(row);continue
        try:
            gt=float(censor.predict_proba([t])[0]);gy=censor.predict_proba(y[case]) if case.any() else np.array([])
        except ValueError:
            row['AUC_reason']=row['Brier_reason']='censoring_support_insufficient';rows.append(row);continue
        if not np.isfinite(gt) or gt<=0 or not np.isfinite(gy).all() or (gy<=0).any():
            row['AUC_reason']=row['Brier_reason']='nonpositive_censoring_survival';rows.append(row);continue
        row['censor_survival_t']=gt
        w=1/gy
        row['Brier']=float(((prob[case,j]**2*w).sum()+((1-prob[control,j])**2/gt).sum())/len(y))
        if not case.any():row['AUC_reason']='no_event_by_horizon'
        elif not control.any():row['AUC_reason']='no_control_beyond_horizon'
        else:
            dif=risk[case,None]-risk[None,control]
            concord=(dif>1e-8).astype(float)+.5*(np.abs(dif)<=1e-8)
            row['AUC']=float((concord*w[:,None]).sum()/(w.sum()*control.sum()))
        rows.append(row)
    return rows

def km_event(time,event,t):
    """Ordinary KM and pointwise Greenwood log-log CI, events before censor ties."""
    y=np.asarray(time,float);e=np.asarray(event,bool);surv=1.;green=0.
    for u in np.unique(y[(y<=t)&e]):
        n=(y>=u).sum();d=((y==u)&e).sum()
        if n==d:surv=0.;break
        surv*=1-d/n;green+=d/(n*(n-d))
    if t>y.max() and surv>0:raise NotEvaluable('group_followup_ends_before_horizon')
    if surv in (0.,1.):lo=hi=surv;boundary=True
    else:
        se=np.sqrt(green)/abs(np.log(surv));z=norm.ppf(.975);ll=np.log(-np.log(surv))
        lo=float(np.exp(-np.exp(ll+z*se)));hi=float(np.exp(-np.exp(ll-z*se)));boundary=False
    return {'observed_km_risk':1-surv,'observed_CI_lower':1-hi,'observed_CI_upper':1-lo,'boundary_interval':boundary}

def calibration(ids,event,time,risk,t):
    g=pd.DataFrame({'patient_id':ids,'event':np.asarray(event,bool),'followup':time,'predicted_risk':risk})
    g['bin']=pd.qcut(g.predicted_risk,5,labels=False,duplicates='drop')+1
    rows=[]
    for k,x in g.groupby('bin'):
        row={'bin':int(k),'n':len(x),'unique_n':x.patient_id.nunique(),'mean_predicted_risk':float(x.predicted_risk.mean()),'median_predicted_risk':float(x.predicted_risk.median()),'events_by_t':int(((x.followup<=t)&x.event).sum()),'events_total':int(x.event.sum()),'at_risk_t':int((x.followup>=t).sum()),'max_followup_months':float(x.followup.max()),'status':'evaluated','reason':''}
        try:
            row.update(km_event(x.followup,x.event,t));row['absolute_calibration_error']=abs(row['mean_predicted_risk']-row['observed_km_risk'])
        except NotEvaluable as exc:row.update(status='NA',reason=str(exc))
        rows.append(row)
    if len(rows)!=5:return rows,None,None,'fewer_than_five_distinct_groups'
    if any(r['status']!='evaluated' for r in rows):return rows,None,None,'group_followup_ends_before_horizon'
    unweighted=float(np.mean([r['absolute_calibration_error'] for r in rows]))
    weighted=float(np.average([r['absolute_calibration_error'] for r in rows],weights=[r['n'] for r in rows]))
    return rows,unweighted,weighted,''

def completeness_gate(expected,records,failures):
    ids=[(x['mode'],x['run']) for x in records]
    if len(ids)!=len(set(ids)) or set(ids)!=set(expected) or failures or any(x['status']!='SUCCESS' for x in records):
        raise RuntimeError('PARTIAL/FAILED: missing, duplicate or failed execution; no completion certificate')
