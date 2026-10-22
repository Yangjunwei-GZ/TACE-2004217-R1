from pathlib import Path
import os,json,hashlib,sys,importlib.metadata,platform
import radiomics,SimpleITK
root=Path(os.environ.get('TACE_PROJECT_ROOT',str(Path.cwd())))/'revision_R1_20261014'
wheel=root.parent/'.runtime/wheelhouse_py310_fix/pyradiomics_fix-3.0.1.dev4-cp310-cp310-win_amd64.whl'
data={'python':platform.python_version(),'executable':sys.executable,'pyradiomics_distribution':importlib.metadata.version('pyradiomics-fix'),'module_version_fallback':radiomics.__version__,'module_file':radiomics.__file__,'SimpleITK':SimpleITK.Version_VersionString(),'wheel_sha256':hashlib.sha256(wheel.read_bytes()).hexdigest(),'pydeps_imported':any('.pydeps' in p for p in sys.path)}
(root/'02_extraction/RUNTIME_PROVENANCE.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
# Record previously verified default values in the revision YAML only; no extraction is executed.
p=root/'02_extraction/extraction_parameters.yaml'
text=p.read_text(encoding='utf-8')
if 'voxelArrayShift:' not in text:text+='  voxelArrayShift: 0\n  correctMask: false\n  geometryTolerance: null\n# Mask interpolation: nearest neighbour in physical coordinates.\n# SimpleITK coordinate and direction tolerances: 1e-6, not overridden.\n'
p.write_text(text,encoding='utf-8')
print(json.dumps(data,indent=2),flush=True)
