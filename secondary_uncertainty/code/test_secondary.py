from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from secondary_common import *
from secondary_metrics import *
from sksurv.metrics import brier_score,cumulative_dynamic_auc,integrated_brier_score
from sksurv.linear_model import CoxPHSurvivalAnalysis
from lifelines import KaplanMeierFitter

class TestSecondary(unittest.TestCase):
 def setUp(self):
    self.et=np.array([1,0,1,1,0,1,0,1,0,1],bool);self.tt=np.array([2,4,5,8,10,13,16,20,30,45.])
    self.ev=np.array([1,0,1,0,1,1],bool);self.tv=np.array([3,7,11,18,28,40.]);self.r=np.array([2,1.5,1,.5,0,-.5])
    self.grid=np.array([6.,12.,24.]);self.prob=np.exp(-np.exp(self.r[:,None])*self.grid[None,:]/60)
 def test_metric_agreement_with_library(self):
    rows=ipcw_metrics(self.et,self.tt,self.ev,self.tv,self.r,self.prob,self.grid)
    yt=Surv.from_arrays(self.et,self.tt);yv=Surv.from_arrays(self.ev,self.tv)
    bs=brier_score(yt,yv,self.prob,self.grid)[1];auc=cumulative_dynamic_auc(yt,yv,self.r,self.grid)[0]
    np.testing.assert_allclose([r['Brier'] for r in rows],bs,atol=1e-12)
    np.testing.assert_allclose([r['AUC'] for r in rows],auc,atol=1e-12)
    self.assertAlmostEqual(np.trapz(bs,self.grid)/(self.grid[-1]-self.grid[0]),integrated_brier_score(yt,yv,self.prob,self.grid),12)
 def test_brier_probability_direction_uncensored(self):
    e=np.ones(6,bool);t=np.array([2.,4,6,8,10,12]);p=np.array([[0],[0],[0],[1],[1],[1]])
    rows=ipcw_metrics(e,t,e,t,-t,p,[7.]);self.assertEqual(rows[0]['Brier'],0);self.assertEqual(rows[0]['AUC'],1)
    wrong=ipcw_metrics(e,t,e,t,t,1-p,[7.]);self.assertEqual(wrong[0]['Brier'],1);self.assertEqual(wrong[0]['AUC'],0)
 def test_no_cases_AUC_NA_brier_evaluable(self):
    rows=ipcw_metrics(self.et,self.tt,self.ev,self.tv,self.r,np.ones((6,1)),[1.])
    self.assertIsNone(rows[0]['AUC']);self.assertEqual(rows[0]['AUC_reason'],'no_event_by_horizon');self.assertEqual(rows[0]['Brier'],0)
 def test_censor_zero_NA(self):
    rows=ipcw_metrics([1,0,0],[1.,2,3],[1,0,1],[1.,3,5],[1.,0,-1],[[.2],[.5],[.8]],[4.])
    self.assertIsNone(rows[0]['Brier']);self.assertTrue(rows[0]['Brier_reason'])
 def test_probability_guards(self):
    for x in [np.array([[1.1,.5]]),np.array([[.2,.3]]),np.array([[np.nan,.2]])]:
        with self.assertRaises(ValueError):probability_check(x)
 def test_baseline_reconstruction_matches_fitted_model(self):
    x=np.column_stack([np.linspace(-1,1,10),np.sin(np.arange(10))]);fit=CoxPHSurvivalAnalysis(alpha=.1).fit(x,Surv.from_arrays(self.et,self.tt))
    pred,_=reconstruct_survival(x@fit.coef_,self.et,self.tt,x[:4]@fit.coef_,self.grid)
    expected=np.array([f(self.grid) for f in fit.predict_survival_function(x[:4])]);np.testing.assert_allclose(pred,expected,atol=1e-12)
 def test_km_loglog_matches_lifelines(self):
    km=KaplanMeierFitter().fit(self.tv,self.ev)
    for t in self.grid:
        z=km_event(self.tv,self.ev,t);ci=km.confidence_interval_survival_function_.loc[:t].iloc[-1]
        self.assertAlmostEqual(z['observed_km_risk'],1-km.predict(t),12)
        self.assertAlmostEqual(z['observed_CI_lower'],1-ci.iloc[1],12);self.assertAlmostEqual(z['observed_CI_upper'],1-ci.iloc[0],12)
 def test_group_unsupported_NA(self):
    with self.assertRaises(NotEvaluable):km_event([1,2,3],[1,0,0],4)
 def test_499_failures_block_completion(self):
    with self.assertRaises(RuntimeError):completeness_gate([('bootstrap',i) for i in range(500)],[{'mode':'bootstrap','run':0,'status':'SUCCESS'}],[{'run':i} for i in range(1,500)])
 def test_missing_duplicate_or_failed_block_completion(self):
    expected=[('point',0),('point',1)]
    for records in [[{'mode':'point','run':0,'status':'SUCCESS'}],[{'mode':'point','run':0,'status':'SUCCESS'}]*2,[{'mode':'point','run':0,'status':'SUCCESS'},{'mode':'point','run':1,'status':'FAILED'}]]:
        with self.assertRaises(RuntimeError):completeness_gate(expected,records,[])
 def test_complete_tasks_pass(self):completeness_gate([('point',0)],[{'mode':'point','run':0,'status':'SUCCESS'}],[])

if __name__=='__main__':unittest.main(verbosity=2)
