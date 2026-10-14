"""Isolated V18 response expert definitions; exact V16 population generator retained."""
from pathlib import Path
import json
import numpy as np
import build_submission as b
ROOT=b.ROOT;OUT=ROOT/'experiments/v18_pie_expert';A=OUT/'analysis'
NAMES=['v16','pie25','confidence','amplitude','shuffled']
FAMILIES=['H1','HCT116','HEK293T']
def shuffled_indices(n):
 order=np.random.default_rng(18018).permutation(n);out=np.empty(n,dtype=int);out[order]=np.roll(order,1);assert np.all(out!=np.arange(n));return out

def blend(base,pie,confidence,support,own):
 base=np.asarray(base,np.float32);pie=np.asarray(pie,np.float32);confidence=np.asarray(confidence,np.float32);support=np.asarray(support,bool).copy()
 assert base.shape==pie.shape==confidence.shape==support.shape and np.isfinite(base).all() and np.isfinite(pie).all() and np.isfinite(confidence).all()
 assert ((confidence>=0)&(confidence<=1)).all();support[own]=False
 w=.25*np.clip((confidence-.5)/.5,0,1)*support
 out=base*(1-w)+pie*w;plain=np.where(support,.75*base+.25*pie,base)
 mask=np.arange(len(base))!=own;scale=float(abs(out[mask]).sum()/max(abs(base[mask]).sum(),1e-12));matched=base*scale
 for value in [out,plain,matched]:value[own]=base[own]
 return {'v16':base.copy(),'pie25':plain.astype(np.float32),'confidence':out.astype(np.float32),'amplitude':matched.astype(np.float32)}

def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V18 dependency changed: '+name
 return p

def load_predictions(family,stage='development'):
 f=A/'inference'/f'{family}_{stage}.npz';rec=json.loads(f.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(f) and rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(f,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
