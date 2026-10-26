"""Synthetic adapter checks; none is a scientific cohort analysis."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import copy, json, tempfile, unittest, subprocess
import numpy as np
import pandas as pd
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.util import Surv
import secondary_common as common
from secondary_common import A,B,draw,verified_folds,rad_score,saved_survival
import run_new_secondary as entry
import run_secondary as execution
import summarize_secondary as summary
import sys
RUNTIME=Path(__file__).resolve().parent.parent/'code/model_runtime'
sys.path.insert(0,str(RUNTIME))
from survival_evidence import fitted_baseline,monthly_survival

AUDIT=Path(__file__).resolve().parent.parent/'audit/secondary_synthetic_tests'
AUDIT.mkdir(parents=True,exist_ok=True)

class TestNewAdapter(unittest.TestCase):
    def setUp(self):
        self.e=np.array([1,0,1,1,0,1,0,1,0,1],bool)
        self.t=np.array([2,4,5,8,10,13,16,20,30,45.])
        self.x=np.column_stack([np.linspace(-1,1,10),np.sin(np.arange(10))])
        self.fit=CoxPHSurvivalAnalysis(alpha=.1).fit(self.x,Surv.from_arrays(self.e,self.t))
        self.base=fitted_baseline(self.fit,self.t)
        self.risk=self.fit.predict(self.x[:4]);self.grid=np.arange(1,37,dtype=float)
        self.record={'fitted_baseline':self.base,'monthly_survival':monthly_survival(self.base,self.risk)}
    def test_saved_probability_and_independent_baseline_agree(self):
        prob,base,supported,err=saved_survival(self.record,self.fit.predict(self.x),self.e,self.t,self.risk,self.grid)
        supported_by_step_domain=self.grid>=self.base['times'][0]
        expected=np.array([f(self.grid[supported_by_step_domain]) for f in self.fit.predict_survival_function(self.x[:4])])
        np.testing.assert_allclose(prob[:,supported_by_step_domain],expected,atol=1e-12)
        # sksurv StepFunction extends its first stored value before the first
        # observation. The Cox baseline must instead be H0(t)=0 there.
        np.testing.assert_allclose(prob[:,~supported_by_step_domain],1.,atol=1e-12)
        self.assertTrue(supported.all());self.assertLess(err,1e-12)
    def test_corrupt_saved_probability_rejected(self):
        r=copy.deepcopy(self.record);r['monthly_survival']['survival_probability'][0][11]=.123456
        with self.assertRaisesRegex(ValueError,'differs'):saved_survival(r,self.fit.predict(self.x),self.e,self.t,self.risk,self.grid)
    def test_wrong_month_grid_rejected(self):
        r=copy.deepcopy(self.record);r['monthly_survival']['months'][-1]=37
        with self.assertRaisesRegex(ValueError,'grid/support'):saved_survival(r,self.fit.predict(self.x),self.e,self.t,self.risk,self.grid)
    def test_old_checkpoint_without_probability_not_reused(self):
        with self.assertRaisesRegex(ValueError,'lacks'):saved_survival({},self.fit.predict(self.x),self.e,self.t,self.risk,self.grid)
    def test_unsupported_months_retained_as_null_evidence(self):
        t=self.t/2;model=CoxPHSurvivalAnalysis(alpha=.1).fit(self.x,Surv.from_arrays(self.e,t))
        risk=model.predict(self.x[:4]);base=fitted_baseline(model,t)
        r={'fitted_baseline':base,'monthly_survival':monthly_survival(base,risk)}
        _,_,supported,_=saved_survival(r,model.predict(self.x),self.e,t,risk,self.grid)
        self.assertEqual(int(supported.sum()),22);self.assertIsNone(r['monthly_survival']['survival_probability'][0][22])
    def test_coxnet_offset_in_radscore(self):
        x=self.x+np.array([2.,4.]);fit=CoxnetSurvivalAnalysis(alphas=[.1],l1_ratio=.5,fit_baseline_model=True).fit(x,Surv.from_arrays(self.e,self.t))
        tr={'columns':['f1','f2'],'median':[0.,0.],'mean':[0.,0.],'sd':[1.,1.],'keep':[True,True]}
        rd={'coefficients':fit.coef_.ravel().tolist(),'transform':tr,'fitted_baseline':fitted_baseline(fit,self.t)}
        d=pd.DataFrame(x,columns=['f1','f2'])
        self.assertGreater(abs(float(fit.offset_[0])),1e-7)
        np.testing.assert_allclose(rad_score(d,rd),fit.predict(x),atol=1e-12)
    def test_original_seed_and_grouped_duplicate_partitions(self):
        from sklearn.model_selection import KFold
        m=pd.DataFrame({'patient_id':[str(i) for i in range(30)],'death_event':[1,0]*15,'os_months':np.arange(1,31.)})
        seed=20261015;ix=np.random.default_rng(seed).integers(0,30,30);d=m.iloc[ix].reset_index(drop=True)
        cp={'mode':'bootstrap','seed':seed,'draw_source_row_indices':ix.tolist(),'draw_patient_ids':d.patient_id.tolist(),'folds':[]}
        unique=np.array(sorted(d.patient_id.unique()))
        for k,(it,iv) in enumerate(KFold(5,shuffle=True,random_state=seed).split(unique)):
            tr=d[d.patient_id.isin(unique[it])];va=d[d.patient_id.isin(unique[iv])]
            cp['folds'].append(dict(fold=k,train_ids=tr.patient_id.tolist(),validation_ids=va.patient_id.tolist(),train_draw_row_indices=tr.index.tolist(),validation_draw_row_indices=va.index.tolist(),training_events=tr.death_event.tolist(),training_times=tr.os_months.tolist()))
        with tempfile.TemporaryDirectory(dir=AUDIT) as temp:
            path=Path(temp)/'config.json';path.write_text('{"outer_folds":5}')
            with patch.object(common,'CONFIG',path):
                folds=list(verified_folds(m,cp));self.assertEqual(len(folds),5)
                for _,tr,va in folds:self.assertFalse(set(tr.patient_id)&set(va.patient_id))
                corrupted=copy.deepcopy(cp);corrupted['folds'][0]['validation_draw_row_indices'][0]=999
                with self.assertRaisesRegex(ValueError,'partition'):list(verified_folds(m,corrupted))
        cp['draw_source_row_indices'][0]=(ix[0]+1)%30
        with self.assertRaisesRegex(ValueError,'seed'):draw(m,cp)
    def test_paired_percentile_uses_common_450_and_preserves_NA(self):
        with tempfile.TemporaryDirectory(dir=AUDIT) as temp:
            root=Path(temp);out=root/'results';out.mkdir();(root/'logs').mkdir()
            tasks=[{'mode':m,'run':k,'status':'SUCCESS'} for m,k in entry.EXPECTED]
            pd.DataFrame(tasks).to_csv(out/'task_ledger.csv',index=False);(root/'logs/PROGRAM_FAILURES.jsonl').write_text('')
            rows=[]
            for mode,k in entry.EXPECTED:
                av=.3+k*.0001;bv=.28+k*.00015
                for name,value in [(A,av),(B,bv)]:
                    rows.append(dict(mode=mode,run=k,model=name,time=12,metric='Brier',value=np.nan if mode=='bootstrap' and name==B and k<50 else value,reason='synthetic_NA' if mode=='bootstrap' and name==B and k<50 else ''))
            pd.DataFrame(rows).to_csv(out/'run_metrics.csv',index=False)
            pd.DataFrame([dict(mode='point',run=0,model=A,time=12,bin=1)]).to_csv(out/'calibration_group_details.csv.gz',index=False)
            with patch.object(summary,'ROOT',root):summary.main()
            s=pd.read_csv(out/'uncertainty_summary.csv').iloc[0];self.assertEqual(s.bootstrap_joint_evaluable,450);self.assertEqual(s.bootstrap_NA,50)
            delta=np.array([(.28+k*.00015)-(.3+k*.0001) for k in range(50,500)])
            np.testing.assert_allclose([s.delta_CI_lower,s.delta_CI_upper],np.quantile(delta,[.025,.975]),atol=1e-12)
            self.assertFalse((out/'STATISTICS_COMPLETE.json').exists())
    def test_config_and_frozen_input_drift_refused(self):
        cfg={'point_repeats':10,'bootstrap_replicates':500,'outer_folds':5,'inner_folds':3,'seed_point':20261014,'seed_bootstrap':20261015,'primary_comparison':[A,B],'prediction_months':[12,24,36],'ibs_monthly_grid':[1,36],'baseline_penalizer':.1}
        with tempfile.TemporaryDirectory(dir=AUDIT) as temp:
            root=Path(temp);rev=root/'revision';out=root/'secondary';runtime=root/'runtime';runtime.mkdir();(runtime/'runtime.py').write_text('x = 1\n')
            (rev/'03_model_development').mkdir(parents=True);(rev/'00_protocol_and_provenance').mkdir()
            matrix=rev/'03_model_development/analysis_matrix.csv';matrix.write_text('patient_id,death_event,os_months\n1,1,12\n2,0,20\n')
            cp=rev/'00_protocol_and_provenance/revision_config.json';cp.write_text(json.dumps(cfg))
            args=SimpleNamespace(revision=str(rev),runtime=str(runtime),output=str(out))
            entry.freeze(args);entry.current_plan(args)
            matrix.write_text(matrix.read_text()+'3,0,30\n')
            with self.assertRaisesRegex(ValueError,'matrix'):entry.current_plan(args)
        cfg['seed_bootstrap']+=1
        with self.assertRaisesRegex(ValueError,'seed_bootstrap'):entry.require_locked_config(cfg)
    def test_real_isolated_cli_imports_before_input_gate(self):
        with tempfile.TemporaryDirectory(dir=AUDIT) as temp:
            root=Path(temp)
            command=[sys.executable,'-B','-X','utf8',str(Path(entry.__file__)), 'compute',
                     '--revision',str(root/'missing_revision'),'--runtime',str(RUNTIME),
                     '--output',str(root/'missing_secondary_plan')]
            r=subprocess.run(command,cwd=str(AUDIT.parent),capture_output=True,text=True,encoding='utf-8')
            self.assertEqual(r.returncode,1)
            self.assertNotIn('ModuleNotFoundError',r.stderr)
            self.assertIn('SUPPLEMENTARY_ANALYSIS_PLAN.json',r.stderr)
            self.assertFalse((root/'missing_secondary_plan/results/STATISTICS_COMPLETE.json').exists())

if __name__=='__main__':unittest.main(verbosity=2)
