"""Final-panel V90 method: isolated source adaptation, unchanged hyperparameters and emitter."""
from pathlib import Path
from functools import lru_cache
import importlib.util,json
import numpy as np
from v57_pie import read,write,sha
from v88_mixture_quantiles import gaussian
from v63_identity_guard import guard
R=Path(__file__).resolve().parent;O=R/'experiments/final_v91';A=O/'analysis'
def check():
 p=read(A/'plan.json')
 for name,h in p['code_sha256'].items():assert sha(R/name)==h,name
 for name,h in p['source_receipt_sha256'].items():assert sha(R/name)==h,name
 return p
@lru_cache(1)
def data():
 rec=read(A/'source_provenance.json');assert sha(A/'source_cache.npz')==rec['cache_sha256']
 with np.load(A/'source_cache.npz',allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def isolated(name,file,transform=None):
 spec=importlib.util.spec_from_file_location(name,R/file);m=importlib.util.module_from_spec(spec)
 if transform is None:spec.loader.exec_module(m)
 else:exec(compile(transform((R/file).read_text(encoding='utf-8')),str(R/file),'exec'),m.__dict__)
 return m
@lru_cache(1)
def assay():
 m=isolated('final_assay_isolated','v62_conserved_response.py');m.data=data;return m
@lru_cache(1)
def posterior():
 m=isolated('final_posterior_isolated','v86_joint_posterior.py',lambda s:s.replace("ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'","A/'source_cache.npz'").replace('Whole recipient study excluded.','All public source families retained for blind final inference; no final treatment outcomes supplied.'));m.A=A;m.OUT=O;m.data=data;m.check=check;m.component_likelihood=gaussian;return m
@lru_cache(1)
def teacher():
 rec=read(A/'cd4_collection.json')
 for name,h in rec['combined'].items():assert sha(R/'reference/final_v91/cd4'/name)==h,name
 m=isolated('final_cd4_isolated','v24_methods.py',lambda s:s.replace('Within-condition donor evidence; Rest development overlaps public summary biology and is calibration, not blind donor validation. Fresh-condition response rows excluded.','All public conditions retained for blind final production; donor evidence is reused public biology, not independent final validation.'));m.A=A;m.read_cache=data;m.load_plan=check;original=m.load_subset
 def subset(directory,name):return original(R/'reference/final_v91/cd4',name)
 m.load_subset=subset;return m,m.make_teacher(True)
def predict(query,bulk):
 check();m,t=teacher();d=data();base,diag=m.predict(d,t,query,production=True);v24=base['cdf50_de'];parent=assay();raw=parent.predict(query,v24)['assay_center50'];proposal=np.clip(v24.astype(float)+1.5*(raw-v24.astype(float)),-3,3).astype(np.float32);own=np.array([list(d['genes']).index(x) for x in d['targets']]);panel=np.isin(d['genes'],d['targets']);proposal[np.arange(300),own]=0;v66,tr66=guard(v24,proposal,bulk,own,panel);out,tr90=posterior().predict(query,v66,bulk,'production');v90=out['interval25'];assert np.isfinite(v90).all() and (v90[np.arange(300),own]==0).all();return v90,dict(v24=diag,v66=tr66,v90=tr90['interval25']),dict(v24=v24,v66=v66,posterior25=out['posterior25'],shared25=out['shared25'],interval25_shuffled=out['interval25_shuffled'])
