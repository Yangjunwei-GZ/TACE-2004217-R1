"""Apply published frozen decisions; verify exact identities before any model fit."""
from pathlib import Path
import csv,hashlib,importlib.metadata,json,math,os,platform,shutil,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from release_identity import verify_release, PACKAGE
from analysis_input_guard import verify_analysis_inputs
REVISION='revision_R1_20261014'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f)
  if not r.fieldnames or len(r.fieldnames)!=len(set(r.fieldnames)):raise ValueError('Missing/duplicate CSV headers: '+str(p))
  return list(r)
def write(p,rs,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields or list(rs[0]));w.writeheader();w.writerows(rs)
def dump(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
def true(v):return str(v).lower()=='true'
def unique(rs,key,label):
 ids=[r[key] for r in rs]
 if not ids or any(not i for i in ids) or len(ids)!=len(set(ids)):raise ValueError('Missing/duplicate '+label+' IDs')
 return {r[key]:r for r in rs}
def verify_rosters(package=PACKAGE):
 package=Path(package);f=package/'frozen_inputs'
 images=rows(f/'FINAL_IMAGE_AND_ALIAS_ROSTER.csv');reps=rows(f/'FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv');candidates=rows(f/'FINAL_CANDIDATE_DISPOSITION.csv');mapping=rows(f/'SOURCE_FILE_MAP.csv')
 im=unique(images,'sample_id','ROI');rp=unique(reps,'patient_id','representative');cp=unique(candidates,'patient_id','candidate');mp=unique(mapping,'sample_id','mapping')
 active={k:v for k,v in im.items() if true(v['active_analysis_record'])};expected={k for k,v in cp.items() if v['eligibility']=='include'}
 if set(rp)!=expected or {v['patient_id'] for v in active.values()}!=expected or set(im)!=set(mp):raise ValueError('Patient/ROI/candidate mapping sets differ')
 for sid,r in im.items():
  m=mp[sid]
  for a,b in [('patient_id','patient_id'),('ct_sha256','ct_sha256'),('original_mask_sha256','original_mask_sha256'),('final_mask_sha256','final_mask_sha256'),('action','frozen_action'),('alias_of','frozen_alias_of'),('active_analysis_record','active_analysis_record')]:
   if r[a]!=m[b]:raise ValueError('Wrong frozen CT/mask/decision mapping: '+sid+' '+a)
  if sid in active and (r['action'] not in ['retain','replace','voxel_edit'] or r['alias_of'] not in ['','未提供']):raise ValueError('Inactive/alias record cannot be an active ROI: '+sid)
  if m['final_mask_kind']=='packaged_final':
   p=(package/m['final_mask_package_path']).resolve()
   if package.resolve() not in p.parents or not p.is_file() or sha(p)!=m['final_mask_sha256']:raise ValueError('Missing or mismatched final mask: '+sid)
 for pid,r in rp.items():
  if r['sample_id'] not in active or active[r['sample_id']]['patient_id']!=pid:raise ValueError('Representative is not an active ROI: '+pid)
  ranked=sorted((x for x in active.values() if x['patient_id']==pid),key=lambda x:(-float(x['physical_volume_mm3']),x['sample_id']))[0]
  if ranked['sample_id']!=r['sample_id']:raise ValueError('Frozen representative differs from the physical-volume rule: '+pid)
  if r['lesions_number']!=cp[pid]['source_clinical_lesion_count']:raise ValueError('Clinical tumour count differs from frozen source: '+pid)
 return images,reps,candidates,mapping
def compare_csv(reference,actual,atol,rtol,ignore=('ct_path','mask_path')):
 a,b=rows(reference),rows(actual)
 if len(a)!=len(b) or not a or list(a[0])!=list(b[0]):raise ValueError('CSV row count/schema mismatch: '+str(actual))
 maximum_abs=maximum_rel=0.0;mismatches=[]
 for i,(x,y) in enumerate(zip(a,b)):
  for k in x:
   if k in ignore:continue
   if k in ['patient_id','sample_id'] and x[k]!=y[k]:raise ValueError('ID or row order changed at '+str(i)+' '+k)
   if x[k]==y[k]:continue
   if x[k]=='' or y[k]=='':mismatches.append((i,k,x[k],y[k]));continue
   try:
    u,v=float(x[k]),float(y[k]);delta=abs(v-u)
    if not math.isfinite(u) or not math.isfinite(v):raise ValueError('nonfinite')
    maximum_abs=max(maximum_abs,delta);maximum_rel=max(maximum_rel,delta/max(abs(u),1e-12))
    if delta>atol+rtol*abs(u):mismatches.append((i,k,x[k],y[k]))
   except ValueError:mismatches.append((i,k,x[k],y[k]))
 if mismatches:raise ValueError('Frozen numeric/nonpath input mismatch: '+json.dumps(mismatches[:15]))
 for r in b:
  for k,v in r.items():
   if k.startswith('original_') and not math.isfinite(float(v)):raise ValueError('Nonfinite radiomics feature: '+k)
 return {'rows':len(a),'columns':len(a[0]),'max_absolute_error':maximum_abs,'max_relative_error':maximum_rel,'atol':atol,'rtol':rtol,'row_order_and_IDs_identical':True,'ignored_relocated_columns':list(ignore)}
def new_run(root):
 root=Path(root).resolve()
 if root==PACKAGE or PACKAGE in root.parents or root in PACKAGE.parents:raise ValueError('Use a new run directory outside the package')
 if root.exists() and any(root.iterdir()):raise ValueError('Use a new empty run directory; never resume old checkpoints here')
 root.mkdir(parents=True,exist_ok=True);rev=root/REVISION
 for n in ['00_protocol_and_provenance','01_data_audit','02_extraction','03_model_development','04_validation','05_figures_and_tables']:(rev/n).mkdir(parents=True)
 shutil.copytree(PACKAGE/'code/protocol',rev/'00_protocol_and_provenance',dirs_exist_ok=True)
 shutil.copytree(PACKAGE/'code/reporting_inputs/01_data_audit',rev/'01_data_audit',dirs_exist_ok=True)
 from legacy_audit_input import ensure_legacy_pipeline_audit
 ensure_legacy_pipeline_audit(root)
 return root,rev
def check_environment():
 if platform.python_version()!='3.10.11':raise RuntimeError('Use tested locked CPython 3.10.11')
 actual={}
 for line in (PACKAGE/'code/environment-python310-executed.lock.txt').read_text().splitlines():
  if '==' in line and not line.startswith('#'):
   n,v=line.split('==',1)
   if n.lower() in ['pip','setuptools','wheel']:continue
   actual[n]=importlib.metadata.version(n)
   if actual[n]!=v:raise RuntimeError('Locked dependency differs: '+n+' expected '+v+' got '+actual[n])
 return {'python':sys.version,'platform':platform.platform(),'packages':actual}
def source_files(source_root,mapping):
 root=Path(source_root).resolve()
 if not root.is_dir():raise ValueError('source-root must be an extracted official-data directory')
 required=set()
 for r in mapping:
  if true(r['active_analysis_record']):
   for key in ['ct','original_mask']:required.add((r[key+'_filename'],r[key+'_sha256']))
 names={x[0] for x in required};found={n:[] for n in names}
 for directory,dirs,files in os.walk(root):
  dirs[:]=[d for d in dirs if d not in {'.git','__pycache__'}]
  for n in names.intersection(files):
   p=(Path(directory)/n).resolve()
   if root not in p.parents:raise ValueError('Source symlink escapes source-root')
   found[n].append(p)
 hashes={};resolved={};rejected=[]
 for name,expected in sorted(required):
  matches=[]
  for p in sorted(found[name]):
   if p not in hashes:hashes[p]=sha(p)
   if hashes[p]==expected:matches.append(p)
   else:rejected.append({'filename':name,'actual_sha256':hashes[p],'expected_sha256':expected})
  if not matches:raise ValueError('Missing or wrong source bytes: '+name+' expected SHA256 '+expected)
  resolved[(name,expected)]=matches[0]
 return resolved,rejected
def prepare(root,source_root=None):
 release=verify_release();images,reps,candidates,mapping=verify_rosters();root,rev=new_run(root)
 source_bindings={};kind='matrix' if source_root is None else 'raw'
 if source_root is None:
  shutil.copy2(PACKAGE/'code/data/analysis_matrix.csv',rev/'03_model_development/analysis_matrix.csv')
  shutil.copytree(PACKAGE/'code/reporting_inputs/02_extraction',rev/'02_extraction',dirs_exist_ok=True)
 else:
  sources,rejected=source_files(source_root,mapping);mp={r['sample_id']:r for r in mapping}
  pair=rows(rev/'01_data_audit/phase_pairing_audit.csv');byid={r['sample_id']:r for r in pair}
  for r in images:
   if not true(r['active_analysis_record']):continue
   m=mp[r['sample_id']];ct=sources[(m['ct_filename'],m['ct_sha256'])];original=sources[(m['original_mask_filename'],m['original_mask_sha256'])]
   source_bindings[str(ct)]=m['ct_sha256'];source_bindings[str(original)]=m['original_mask_sha256']
   if m['final_mask_kind']=='packaged_final':
    source=PACKAGE/m['final_mask_package_path'];mask=root/'final_masks'/source.name;mask.parent.mkdir(exist_ok=True);shutil.copy2(source,mask)
   else:mask=original
   source_bindings[str(mask)]=m['final_mask_sha256'];r.update(ct_path=str(ct),original_mask_path=str(original),final_mask_path=str(mask))
   q=byid[r['sample_id']];q.update(ct_path=str(ct),mask_path=str(mask),eligible='True')
   if q['patient_id']!=r['patient_id'] or q['ct_sha256']!=r['ct_sha256'] or q['mask_sha256']!=r['final_mask_sha256']:raise ValueError('Frozen analysis pairing disagrees: '+r['sample_id'])
  active={r['sample_id'] for r in images if true(r['active_analysis_record'])}
  for q in pair:
   if q['ct_phase']=='1':q['eligible']=str(q['sample_id'] in active)
  write(rev/'01_data_audit/phase_pairing_audit.csv',pair)
  dump(root/'SOURCE_RESOLUTION_AUDIT.json',{'source_root':str(Path(source_root).resolve()),'resolved_unique_required_files':len(sources),'nonmatching_same_name_files':rejected,'application':'Frozen masks copied exactly; no new medical or geometric decision.'})
 for name,data in [('FINAL_IMAGE_AND_ALIAS_ROSTER.csv',images),('FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv',reps),('FINAL_CANDIDATE_DISPOSITION.csv',candidates)]:write(root/name,data)
 binding={'release_manifest_sha256':release['manifest_sha256'],'files':{str(root/n):sha(root/n) for n in ['FINAL_IMAGE_AND_ALIAS_ROSTER.csv','FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv','FINAL_CANDIDATE_DISPOSITION.csv']},'source_files':source_bindings,'route':kind}
 for folder in ['00_protocol_and_provenance','01_data_audit']:binding['files'].update({str(p):sha(p) for p in (rev/folder).rglob('*') if p.is_file()})
 dump(root/'FROZEN_RUN_INPUT_BINDINGS.json',binding)
 if kind=='raw':
  geometry_check(images);verify_bindings(binding)
  dump(root/'INPUT_PREFLIGHT.json',{'status':'PASS_RAW_SOURCE_PREFLIGHT_NOT_EXTRACTION','patients':len(reps),'active_ROIs':len(active),'final_masks_applied':sum(true(r['active_analysis_record']) and r['final_mask_kind']=='packaged_final' for r in mapping),'model_fits':0})
 else:verify_outputs(root,False)
 return root,rev,binding
def geometry_check(images):
 import SimpleITK as sitk
 import numpy as np
 for r in images:
  if not true(r['active_analysis_record']):continue
  ct=sitk.ImageFileReader();ct.SetFileName(r['ct_path']);ct.ReadImageInformation();ma=sitk.ReadImage(r['final_mask_path'])
  if ct.GetSize()!=ma.GetSize() or any(not np.allclose(a,b,atol=1e-4,rtol=0) for a,b in [(ct.GetSpacing(),ma.GetSpacing()),(ct.GetOrigin(),ma.GetOrigin()),(ct.GetDirection(),ma.GetDirection())]):raise ValueError('CT/final mask geometry mismatch: '+r['sample_id'])
  n=int(np.count_nonzero(sitk.GetArrayFromImage(ma)>0));volume=n*float(np.prod(ma.GetSpacing()))
  if n==0 or not math.isclose(volume,float(r['physical_volume_mm3']),abs_tol=1e-8,rel_tol=1e-12):raise ValueError('Empty/changed ROI or volume: '+r['sample_id'])
def verify_bindings(binding):
 current=verify_release()
 if current['manifest_sha256']!=binding['release_manifest_sha256']:raise ValueError('Release changed during execution')
 for p,value in {**binding['files'],**binding['source_files']}.items():
  if not Path(p).is_file() or sha(p)!=value:raise ValueError('Input changed during execution: '+p)
def verify_outputs(root,raw):
 root=Path(root);rev=root/REVISION;protocol=json.loads((PACKAGE/'frozen_inputs/NUMERICAL_COMPARISON_PROTOCOL.json').read_text())
 atol=protocol['feature_and_matrix_atol'] if raw else 0;rtol=protocol['feature_and_matrix_rtol'] if raw else 0;comparisons={}
 for name,ref,actual in [('features',PACKAGE/'code/reporting_inputs/02_extraction/arterial_lesion_features.csv',rev/'02_extraction/arterial_lesion_features.csv'),('matrix',PACKAGE/'code/data/analysis_matrix.csv',rev/'03_model_development/analysis_matrix.csv')]:comparisons[name]=compare_csv(ref,actual,atol,rtol)
 audit=verify_analysis_inputs(root,verify_image_hashes=raw)
 dump(root/'INPUT_REPRODUCTION_COMPARISON.json',{'status':'PASS','route':'raw' if raw else 'matrix','comparisons':comparisons,'new_model_fits':0,'features_reextracted':raw})
 dump(root/'INPUT_PREFLIGHT.json',{'status':'PASS_EXTRACTED_INPUTS' if raw else 'PASS_FROZEN_MATRIX_INPUTS','patients':audit['actual_patient_count'],'active_ROIs':audit['actual_active_unique_mask_count'],'model_fits':0})
 return audit
