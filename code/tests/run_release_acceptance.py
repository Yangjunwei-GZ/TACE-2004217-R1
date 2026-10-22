from pathlib import Path
import argparse,datetime,json,platform,sys,unittest
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parent))
class Results(unittest.TextTestResult):
 def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.records=[]
 def addSuccess(self,test):super().addSuccess(test);self.records.append({'test':test.id(),'status':'PASS'})
 def addFailure(self,test,err):super().addFailure(test,err);self.records.append({'test':test.id(),'status':'FAIL','exception':str(err[1])})
 def addError(self,test,err):super().addError(test,err);self.records.append({'test':test.id(),'status':'ERROR','exception':str(err[1])})
def main():
 p=argparse.ArgumentParser();p.add_argument('--report',required=True,type=Path);a=p.parse_args()
 suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromName(n) for n in ['test_public_inputs','test_validation_integrity']])
 r=unittest.TextTestRunner(verbosity=2,resultclass=Results).run(suite)
 d={'scope':'Unit and SYNTHETIC failure/validation fixtures, not real scientific model fitting','tests_run':r.testsRun,'failures':len(r.failures),'errors':len(r.errors),'status':'PASS' if r.wasSuccessful() else 'FAIL','python':sys.version,'platform':platform.platform(),'actual_package':str(HERE.parent.parent),'cwd':str(Path.cwd()),'records':r.records,'completed_at':datetime.datetime.now().astimezone().isoformat()}
 a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(d,indent=2),encoding='utf-8');return 0 if r.wasSuccessful() else 1
if __name__=='__main__':sys.exit(main())
