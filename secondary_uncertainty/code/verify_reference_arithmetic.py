"""Check current reference arithmetic; not a rerun of model development."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
A='clinical_diameter_mean_hu';B=A+'_radiomics'
def verify(reference):
    root=Path(reference)
    summary=pd.read_csv(root/'uncertainty_summary.csv',dtype={'time':str})
    paired=pd.read_csv(root/'paired_bootstrap_metrics.csv',dtype={'time':str})
    runs=pd.read_csv(root/'run_metrics.csv',dtype={'time':str})
    assert len(summary)==13 and len(paired)==6500 and len(runs)==13260
    details=[]
    def equal(a,b):
        assert (pd.isna(a) and pd.isna(b)) or np.isclose(float(a),float(b),rtol=1e-11,atol=1e-12),(a,b)
    for row in summary.to_dict('records'):
        p=paired[(paired.metric==row['metric'])&(paired.time==row['time'])].sort_values('run')
        assert p.run.tolist()==list(range(500)) and p.seed.tolist()==list(range(20261015,20261515))
        assert np.allclose(p.delta,p.B-p.A,equal_nan=True)
        valid=p.A.notna()&p.B.notna();n=int(valid.sum())
        assert n==int(row['bootstrap_joint_evaluable']) and 500-n==int(row['bootstrap_NA'])
        assert (p.joint_evaluable.to_numpy(bool)==valid.to_numpy()).all()
        r=runs[(runs.metric==row['metric'])&(runs.time==row['time'])]
        assert not r.duplicated(['mode','run','model']).any()
        w=r.pivot(index=['mode','run'],columns='model',values='value')
        assert set(w.index)=={(mode,k) for mode,nr in [('point',10),('bootstrap',500)] for k in range(nr)}
        point=w.loc[('point',0)]
        equal(row['A_point'],point[A]);equal(row['B_point'],point[B]);equal(row['delta_B_minus_A'],point[B]-point[A])
        ci_ok=n>=450 and point[[A,B]].notna().all()
        assert row['CI_status']==('estimated' if ci_ok else 'not_reliably_estimable')
        for key in ['A','B','delta']:
            lo,hi=np.quantile(p.loc[valid,key],[.025,.975],method='linear') if ci_ok else [np.nan,np.nan]
            equal(row[key+'_CI_lower'],lo);equal(row[key+'_CI_upper'],hi)
        details.append(dict(metric=row['metric'],time=row['time'],joint_evaluable=n,NA=500-n,CI_status=row['CI_status']))
    na=pd.read_csv(root/'NON_EVALUABLE_RUNS.csv')
    assert len(na)==int(runs.value.isna().sum())
    return dict(status='PASS_REFERENCE_ARITHMETIC_ONLY',summary_rows=len(summary),paired_rows=len(paired),per_run_rows=len(runs),NA_model_endpoint_rows=len(na),details=details,model_training_performed=False,original_checkpoint_bodies_required_for_this_arithmetic_check=False,scope='Checks numerical intervals from current reference records; does not reproduce full model development')
def main():
    p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    result=verify(a.reference);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
