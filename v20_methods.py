"""Stack baseline-calibrated mean response and conditional program residual emission."""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
import build_submission as b
from v8_methods import gated
from v16_methods import emit as v16_emit
ROOT=b.ROOT;OUT=ROOT/'experiments/v20_stack_transfer';A=OUT/'analysis';R=ROOT/'reference/v20_stack'
FAMILIES=['H1','HCT116','HEK293T'];NAMES=['v16','stack_mean','stack_joint','raw_mean','shuffled','reweight','amplitude']

def align_counts(x,genes,vocab):
 mapping={g:i for i,g in enumerate(vocab)};keep=np.flatnonzero(np.isin(genes,vocab));z=x[:,keep].tocsr();index=np.array([mapping[str(genes[i])] for i in keep]);out=sparse.csr_matrix((z.data.copy(),index[z.indices],z.indptr.copy()),shape=(x.shape[0],len(vocab)));out.sort_indices();return out

def summarize_pair(control,treated,query,model,measured,supported_vocab,vocab,genes,control_theta=None,treated_theta=None):
 """LFC and uncertainty on the same covered vocabulary mass; missing outputs stay masked."""
 common=np.flatnonzero(supported_vocab);mapping={g:i for i,g in enumerate(genes)};comp=np.array([mapping[vocab[i]] for i in common]);a=np.asarray(control[:,common].toarray() if sparse.issparse(control) else control[:,common],float);t=np.asarray(treated[:,common].toarray() if sparse.issparse(treated) else treated[:,common],float);la=np.maximum(a.sum(1,keepdims=True),1);lt=np.maximum(t.sum(1,keepdims=True),1);extra_a=np.zeros(len(common));extra_t=extra_a.copy()
 if control_theta is not None:
  extra_a=np.mean((a+a*a/np.maximum(control_theta[:,common],1e-8))/(la*la),axis=0)/len(a);extra_t=np.mean((t+t*t/np.maximum(treated_theta[:,common],1e-8))/(lt*lt),axis=0)/len(t)
 a/=la;t/=lt
 ma=a.mean(0);mt=t.mean(0);lfc=np.log2((mt+1e-9)/(ma+1e-9));variance=((a.var(0)/len(a)+extra_a)/(ma+1e-9)**2+(t.var(0)/len(t)+extra_t)/(mt+1e-9)**2)/np.log(2)**2;shrink=(lfc*lfc/(lfc*lfc+variance+.0025))*.04/(.04+variance);q=query[comp];q=q/max(q.sum(),1e-12);raw=np.log2((mt+1e-9)/(q+1e-9));g=len(genes);out={k:np.zeros(g,np.float32) for k in ['lfc','shrink','raw']};mask=np.zeros(g,bool);mask[comp]=True;out['lfc'][comp]=np.clip(lfc,-3,3);out['shrink'][comp]=shrink;out['raw'][comp]=np.clip(raw,-3,3)
 # Program inputs are paired natural-log CP10k shifts, centered to separate joint residual from mean.
 positions=model['positions'];lookup={int(c):i for i,c in enumerate(comp)};local=np.array([lookup.get(int(v),-1) for v in positions]);valid=local>=0;d=np.zeros((400,len(positions)));d[:,valid]=np.log1p(t[:,local[valid]]*1e4)-np.log1p(a[:,local[valid]]*1e4);direction=(d.mean(0)/model['scale'])@model['basis'];d-=d.mean(0);program=(d/model['scale'])@model['basis']
 if valid.mean()<.8:program[:]=0;direction[:]=0
 out.update(supported=mask,program_delta=program.astype(np.float32),direction=direction.astype(np.float32),program_coverage=float(valid.mean()));return out

def target_effects(base,expert,ti,query,own):
 mask=expert['supported'][ti].copy();mask[own]=False;delta=.25*expert['shrink'][ti]*np.clip(expert['lfc'][ti]-base,-.5,.5);delta[~mask]=0;calibrated=np.clip(base+gated(delta,np.nan_to_num(query),own),-3,3).astype(np.float32)
 shuffled_delta=.25*expert['shuffled_shrink'][ti]*np.clip(expert['shuffled_lfc'][ti]-base,-.5,.5);shuffled_delta[~mask]=0;shuffled=np.clip(base+gated(shuffled_delta,np.nan_to_num(query),own),-3,3).astype(np.float32);raw=base.copy();raw[mask]=gated(expert['raw'][ti].copy(),np.nan_to_num(query),own)[mask]
 ratio=float(abs(calibrated).sum()/max(abs(base).sum(),1e-12));amplitude=(base*ratio).astype(np.float32)
 result={'v16':base.copy(),'stack_mean':calibrated,'stack_joint':calibrated,'raw_mean':raw,'shuffled':shuffled,'reweight':calibrated,'amplitude':amplitude}
 for value in result.values():value[own]=base[own]
 return result

def conditional_rows(query_inputs,template_rows):
 lookup={int(r):i for i,r in enumerate(query_inputs['query_rows'])};score=query_inputs['query_scores'][[lookup[int(r)] for r in template_rows]];neighbors=cdist(score,query_inputs['cached_scores']).argmin(1);return score,neighbors

def emit(template,effect,own,seed,kind,prepared,atlas,expert=None,ti=None,query_inputs=None,template_rows=None):
 if kind=='stack_joint' and (not np.any(expert['program_delta'][ti]) or not expert['supported'][ti].any()):return v16_emit(template,np.arange(400),effect,own,seed,'residual10',prepared,atlas)
 if kind not in ['stack_joint','reweight']:return v16_emit(template,np.arange(400),effect,own,seed,'residual10',prepared,atlas)
 depth,p,bulk=prepared;libs=np.asarray(depth,dtype=np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();rng=np.random.default_rng(seed);expected=libs[:,None]*wp;poisson=rng.poisson(expected);empirical=template.toarray().astype(float)*factor
 score,nearest=conditional_rows(query_inputs,template_rows);positions=query_inputs['positions']
 if kind=='stack_joint':
  delta=expert['program_delta'][ti][nearest]@query_inputs['basis'].T;gain=.25*np.clip(delta*query_inputs['scale']/np.maximum(query_inputs['derivative'],.1),-1,1);gain[:,~expert['supported'][ti,positions]]=0
  if own is not None:gain[:,positions==own]=0
  empirical[:,positions]*=np.exp2(gain);effective=400.
 else:
  direction=expert['direction'][ti];projection=score@direction
  if np.linalg.norm(direction)<1e-8:return v16_emit(template,np.arange(400),effect,own,seed,'residual10',prepared,atlas)
  projection=(projection-projection.mean())/max(float(projection.std()),1e-8);strength=min(float(np.linalg.norm(direction)/4),.5);weight=np.exp(strength*np.clip(projection,-2,2));weight/=weight.sum();draw=np.random.default_rng(seed+20020).choice(400,400,replace=True,p=weight);empirical=template[draw].toarray().astype(float)*factor;effective=float(1/(weight**2).sum())
 empirical*=np.divide(libs,empirical.sum(1))[:,None];initial=.25*expected+.65*poisson+.1*empirical;initial/=initial.sum(1,keepdims=True);count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200);np.testing.assert_array_equal(count.sum(1,dtype=np.int64),libs);assert np.issubdtype(count.dtype,np.integer) and (count>=0).all();rp=(count/libs[:,None]).mean(0);rb=count.sum(0,dtype=np.int64)/libs.sum();diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'changed_transcript_fraction':float(abs(count-template.toarray()).sum()/2/libs.sum()),'global_moment_refit':True,'state_effective_sample_size':effective};return sparse.csr_matrix(count),diag

def load_plan():
 plan=json.loads((A/'plan.json').read_text())
 for name,h in plan['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V20 dependency changed: '+name
 for item in plan['input_files']:assert b.sha256(ROOT/item['path'])==item['sha256']
 assert b.sha256(A/'assets.json')==plan['assets_manifest_sha256']
 return plan

def load_query(family):
 file=A/'inputs'/f'{family}_query.npz';rec=json.loads((A/'inputs'/f'{family}.json').read_text());assert b.sha256(file)==rec['query_arrays_sha256']
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def load_expert(family,stage):
 file=A/'experts'/f'{family}_{stage}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
