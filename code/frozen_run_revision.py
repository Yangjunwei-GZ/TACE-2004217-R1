"""Frozen-matrix or official-image route. Training requires explicit --train."""
from pathlib import Path
import argparse,datetime,json,os,shutil,subprocess,sys,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
from frozen_input_pipeline import prepare,verify_outputs,verify_bindings,check_environment,dump
from analysis_input_guard import assert_inputs_unchanged,write_completion_marker
from validation_integrity import audit_validation,verify_paired_summary
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project-root',required=True,type=Path)
 mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--from-matrix',action='store_true');mode.add_argument('--source-root',type=Path)
 action=p.add_mutually_exclusive_group();action.add_argument('--preflight-only',action='store_true');action.add_argument('--extract-only',action='store_true');action.add_argument('--train',action='store_true')
 p.add_argument('--workers',type=int,default=2);p.add_argument('--rscript',default=shutil.which('Rscript'));a=p.parse_args()
 root=a.project_root.resolve();code=Path(__file__).resolve().parent;stage='prepare';log=[];created=False
 try:
  if a.workers<1:raise ValueError('workers must be positive')
  if a.extract_only and a.from_matrix:raise ValueError('--extract-only requires --source-root')
  if a.train and (not a.rscript or not (Path(a.rscript).is_file() or shutil.which(a.rscript))):raise RuntimeError('Training requires Rscript with survival; supply --rscript')
  if a.source_root or a.train:environment=check_environment()
  else:environment={'python':sys.version,'scope':'Matrix identity preflight only; dependency lock enforced before extraction/training'}
  if root.exists() and any(root.iterdir()):raise ValueError('Use a new empty run directory')
  created=True;root,rev,binding=prepare(root,a.source_root);dump(root/'ACTUAL_ENVIRONMENT.json',environment)
  env=dict(os.environ,TACE_PROJECT_ROOT=str(root),PYTHONDONTWRITEBYTECODE='1',PYTHONPATH='')
  def run(script,*args):
   nonlocal stage
   stage=script;verify_bindings(binding);command=[sys.executable,str(code/script),*map(str,args)];start=datetime.datetime.now().astimezone().isoformat()
   print('RUN',script,*args,flush=True);r=subprocess.run(command,env=env,cwd=root)
   log.append({'stage':stage,'command':command,'started_at':start,'exit_code':r.returncode});dump(root/'ACTUAL_EXECUTION_LOG.json',log)
   if r.returncode:raise RuntimeError('Stage failed with exit code '+str(r.returncode)+': '+stage)
   verify_bindings(binding)
  if a.preflight_only:
   print(json.dumps({'status':'INPUT_PREFLIGHT_ONLY','model_fits':0,'features_reextracted':False}));return
  if a.source_root:run('raw_feature_extract.py','--workers',a.workers)
  audit=verify_outputs(root,bool(a.source_root));verify_bindings(binding)
  from public_execution_gate import prepare_execution_gate
  audit=prepare_execution_gate(root,bool(a.source_root))
  if not a.train:
   print(json.dumps({'status':'INPUTS_VERIFIED_NO_TRAINING','model_fits':0,'features_reextracted':bool(a.source_root)}));return
  for script,args in [('diagnostics_and_tests.py',[]),('validate_pipeline.py',['--mode','point','--workers',a.workers]),('validate_pipeline.py',['--mode','bootstrap','--workers',a.workers])]:
   assert_inputs_unchanged(root,audit,script);run(script,*args)
  stage='validation_integrity';audit_validation(rev);run('validate_pipeline.py','--mode','summary');verify_paired_summary(rev)
  stage='PH_global';r=subprocess.run([a.rscript,str(code/'ph_global.R'),str(rev/'03_model_development/FINAL_REVISION_PATIENT_SCORES.csv'),str(rev/'03_model_development/PH_GLOBAL_RISK_SET_DIAGNOSTICS.csv')],env=env,cwd=root);log.append({'stage':stage,'exit_code':r.returncode});dump(root/'ACTUAL_EXECUTION_LOG.json',log)
  if r.returncode:raise RuntimeError('R survival diagnostics failed')
  assert_inputs_unchanged(root,audit,'reporting');run('publication_data.py')
  if a.source_root:run('publication_figures.py')
  verify_bindings(binding);assert_inputs_unchanged(root,audit,'completion');marker=write_completion_marker(root,audit)
  from decision_provenance import decision_provenance
  marker.update(scope='MODEL_VALIDATION_ONLY',**decision_provenance(root/'PUBLIC_INPUT_ADOPTION.csv'),medical_evidence_requirements_satisfied=False,submission_ready=False,secondary_uncertainty_complete=False,publication_materials_complete=False)
  dump(root/'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json',marker)
  print('MODEL VALIDATION COMPLETE: new run fingerprints; secondary statistics and document synchronization remain separate. Decision attribution does not certify a signed independent review or submission readiness.')
 except Exception as exc:
  if created and root.exists():dump(root/'EXECUTION_FAILURE.json',{'status':'FAILED','stage':stage,'exception':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc(),'time':datetime.datetime.now().astimezone().isoformat(),'completion_authorized':False})
  print('FAILED '+stage+': '+str(exc),file=sys.stderr);return 1
 return 0
if __name__=='__main__':sys.exit(main() or 0)
