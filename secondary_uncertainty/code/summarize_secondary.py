from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from secondary_common import *
from secondary_metrics import completeness_gate
def main():
    out=ROOT/'results';tasks=pd.read_csv(out/'task_ledger.csv').to_dict('records')
    failures=[json.loads(x) for x in (ROOT/'logs/PROGRAM_FAILURES.jsonl').read_text().splitlines() if x]
    completeness_gate([(mode,k) for mode,n in [('point',10),('bootstrap',500)] for k in range(n)],tasks,failures)
    df=pd.read_csv(out/'run_metrics.csv');rows=[];paired=[];na=[]
    for (metric,t),g in df.groupby(['metric','time'],sort=False):
        wide=g.pivot(index=['mode','run'],columns='model',values='value');bs=wide.loc['bootstrap'];pt=wide.loc[('point',0)]
        joint=bs[MODELS].notna().all(axis=1);valid=bs.loc[joint];n=int(joint.sum());good=n>=450 and pt[MODELS].notna().all()
        point=wide.loc['point'];delta=point[B]-point[A]
        row={'metric':metric,'time':t,'A_point':pt[A],'B_point':pt[B],'delta_B_minus_A':pt[B]-pt[A],'bootstrap_planned':500,'bootstrap_joint_evaluable':n,'bootstrap_NA':500-n,'CI_status':'estimated' if good else 'not_reliably_estimable','point_estimator':'point_0000; five-fold mean for AUC/Brier/IBS; patient-OOF quintiles for calibration','A_ten_repeat_mean':point[A].mean(),'B_ten_repeat_mean':point[B].mean(),'delta_ten_repeat_mean':delta.mean(),'delta_ten_repeat_SD':delta.std(),'point_repeats_joint_evaluable':int(point[MODELS].notna().all(axis=1).sum())}
        for name,key in [(A,'A'),(B,'B'),('delta','delta')]:
            z=valid[B]-valid[A] if key=='delta' else valid[name]
            lo,hi=np.quantile(z,[.025,.975]) if good else [np.nan,np.nan]
            row[key+'_CI_lower']=lo;row[key+'_CI_upper']=hi
        rows.append(row)
        for i,x in bs.iterrows():paired.append({'metric':metric,'time':t,'run':i,'seed':20261015+i,'A':x[A],'B':x[B],'delta':x[B]-x[A],'joint_evaluable':bool(joint.loc[i])})
        for x in g[g.value.isna()].to_dict('records'):na.append(x)
    pd.DataFrame(rows).to_csv(out/'uncertainty_summary.csv',index=False)
    pd.DataFrame(paired).to_csv(out/'paired_bootstrap_metrics.csv',index=False)
    pd.DataFrame(na,columns=df.columns).to_csv(out/'NON_EVALUABLE_RUNS.csv',index=False)
    if na:pd.DataFrame(na).groupby(['mode','model','metric','time','reason'],dropna=False).size().reset_index(name='count').to_csv(out/'NON_EVALUABLE_REASON_COUNTS.csv',index=False)
    cal=pd.read_csv(out/'calibration_group_details.csv.gz');cal[(cal['mode']=='point')&(cal.run==0)].to_csv(out/'point0_calibration_with_KM_CI.csv',index=False)
    print(pd.DataFrame(rows)[['metric','time','A_point','B_point','delta_B_minus_A','delta_CI_lower','delta_CI_upper','bootstrap_joint_evaluable','CI_status']].to_string(index=False))
if __name__=='__main__':raise SystemExit('Use run_new_secondary.py compute/audit; direct summary publication is disabled.')
