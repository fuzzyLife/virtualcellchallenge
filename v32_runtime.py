"""Numerically unchanged baseline caching for V32-V34, no scientific reselection."""
from functools import lru_cache
import json,hashlib
import numpy as np
import v32_methods as m
import v25_methods as v25
import v24_methods as v24

@lru_cache(None)
def teacher(production=False):
 name='teacher_all' if production else 'teacher_development';p=m.ROOT/'experiments/v24_cd4_response/analysis/models'/f'{name}.npz';r=json.loads(p.with_suffix('.json').read_text());assert m.sha(p)==r['sha256'];assert r['plan_sha256']==m.sha(p.parents[1]/'plan.json')
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
original=v25.baseline
memo={}
def baseline(data,query,excluded=(),production=False,parent=None):
 if parent is not None:return parent.copy()
 key=hashlib.sha256(np.asarray(query).tobytes()+str(sorted(excluded)).encode()).hexdigest();folder=m.A/'baseline_cache';folder.mkdir(exist_ok=True);file=folder/(key+'.npz');side=file.with_suffix('.json')
 if key not in memo:
  if side.exists():
   rec=json.loads(side.read_text());assert rec['source_sha256']==m.sha(m.ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz') and rec['teacher_sha256']==m.sha(m.ROOT/'experiments/v24_cd4_response/analysis/models/teacher_all.npz') and m.sha(file)==rec['sha256']
   with np.load(file,allow_pickle=False) as z:memo[key]=z['effect'].copy()
  else:
   memo[key]=original(data,query,excluded,production);np.savez_compressed(file,effect=memo[key]);m.write(side,{'sha256':m.sha(file),'source_sha256':m.sha(m.ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'teacher_sha256':m.sha(m.ROOT/'experiments/v24_cd4_response/analysis/models/teacher_all.npz'),'excluded':list(excluded),'query_sha256':hashlib.sha256(np.asarray(query).tobytes()).hexdigest()})
 return memo[key].copy()
def install():
 v25.make_teacher=teacher;v25.baseline=baseline;m.old.baseline=baseline