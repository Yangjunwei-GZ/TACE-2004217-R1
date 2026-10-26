"""Focused package tests: not new medical evidence or model validation."""
from pathlib import Path
import argparse,datetime,json,platform,sys,unittest
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parent))
class Results(unittest.TextTestResult):
    def __init__(self,*a,**k):super().__init__(*a,**k);self.records=[]
    def addSuccess(self,t):super().addSuccess(t);self.records.append(dict(test=t.id(),status='PASS'))
    def addFailure(self,t,e):super().addFailure(t,e);self.records.append(dict(test=t.id(),status='FAIL',exception=str(e[1])))
    def addError(self,t,e):super().addError(t,e);self.records.append(dict(test=t.id(),status='ERROR',exception=str(e[1])))
def main():
    p=argparse.ArgumentParser();p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromName(n) for n in ['test_public_inputs','test_public_author_gate']])
    r=unittest.TextTestRunner(verbosity=2,resultclass=Results).run(suite)
    d=dict(status='PASS' if r.wasSuccessful() else 'FAIL',scope='Focused public-entry and synthetic fault injection; no scientific training',tests=r.testsRun,failures=len(r.failures),errors=len(r.errors),records=r.records,python=sys.version,platform=platform.platform(),model_fits=0,completed_at=datetime.datetime.now().astimezone().isoformat())
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(d,indent=2),encoding='utf-8')
    return 0 if r.wasSuccessful() else 1
if __name__=='__main__':sys.exit(main())
