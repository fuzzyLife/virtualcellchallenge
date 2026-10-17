"""Experimental RNAi phenotype fingerprints, no CRISPR outcome input."""
import json
from collections import Counter
from pathlib import Path
import h5py,numpy as np,pandas as pd
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits
import build_submission as b
OUT=b.ROOT/'experiments/v43_lincs_response';A=OUT/'analysis';R=OUT/'reference'
def main():
 rec=json.loads((A/'collection.json').read_text())
 for item in rec['files']:assert b.sha256(R/item['file'])==item['sha256']
 un=json.loads((A/'unpacked.json').read_text());fpath=R/'GSE106127_CGS_PRIME_n33839x978.gctx'
 assert b.sha256(fpath)==next(x['sha256'] for x in un['files'] if x['file']==fpath.name)
 meta=pd.read_csv(R/'GSE106127_CGS_PRIME_meta.txt.gz',sep='\t');gene=pd.read_csv(R/'GSE106127_gene_info.txt.gz',sep='\t')
 with h5py.File(fpath,'r') as f:
  ids=f['0/META/COL/id'].asstr()[:];gids=f['0/META/ROW/id'].asstr()[:];values=f['0/DATA/0/matrix'][:]
 meta=meta.set_index('sig_id').loc[ids].reset_index();gm=gene.set_index(gene.pr_gene_id.astype(str));syms=gm.loc[gids,'pr_gene_symbol'].astype(str).to_numpy()
 _,genes,panel=b.read_spec(b.ROOT);unique=Counter(syms);keep=np.array([unique[s]==1 and s not in panel for s in syms]);syms=syms[keep];values=values[:,keep].astype(float)
 ok=(meta.pert_type=='trt_sh.cgs')&(meta.distil_nsample>=3)&np.isfinite(meta.distil_cc_q75)&(meta.distil_cc_q75>0)
 contexts=sorted(meta.loc[ok,'cell_id'].unique());names=np.array(sorted(meta.loc[ok,'pert_iname'].unique()),str);lookup={t:i for i,t in enumerate(names)}
 num=np.zeros((len(names),len(syms)));den=np.zeros(len(names));sign=np.zeros_like(num);n=np.zeros(len(names),int);used=[]
 for c in contexts:
  rows=np.flatnonzero(ok&(meta.cell_id==c));w=np.clip(meta.distil_cc_q75.to_numpy(float)[rows],0,1)**2;center=np.average(values[rows],axis=0,weights=w)
  for i,ww in zip(rows,w):
   t=str(meta.pert_iname.iloc[i]);ti=lookup[t];v=np.clip(values[i]-center,-8,8);v[syms==t]=0
   num[ti]+=ww*v;sign[ti]+=ww*np.sign(v);den[ti]+=ww;n[ti]+=1;used.append(str(ids[i]))
 good=n>=3;names=names[good];n=n[good];mu=num[good]/den[good,None];agree=abs(sign[good])/den[good,None];finger=np.clip(mu,-5,5)*agree
 center=finger.mean(0);scale=np.maximum(finger.std(0),.1);x=(finger-center)/scale
 u,s,vt=randomized_svd(x,n_components=64,n_iter=7,random_state=43043)
 vectors=u*s;vectors/=np.maximum(np.linalg.norm(vectors,axis=1,keepdims=True),1e-12)
 assert np.isfinite(vectors).all() and len(set(names))==len(names);assert not np.isin(syms,panel).any()
 file=A/'experimental_fingerprints.npz';assert not file.exists();np.savez_compressed(file,keys=names,vectors=vectors.astype(np.float32),landmarks=np.asarray(syms,str),contexts=np.array(contexts),context_support=n,center=center,scale=scale,components=vt,singular_values=s)
 b.write_json(A/'fingerprints.json',{'arrays_sha256':b.sha256(file),'script_sha256':b.sha256(Path(__file__)),'collection_sha256':b.sha256(A/'collection.json'),'targets':len(names),'panel_targets':int(np.isin(panel,names).sum()),'landmarks':len(syms),'dimensions':64,'contexts':contexts,'signature_ids':used,'source':'LINCS GSE106127 CGS PRIME RNAi, >=3 hairpins and positive consistency, >=3 contexts/target','transformation':'Cell-line center using all eligible interventions, quality^2 weighted consensus, downweight sign disagreement, standardize landmarks and PCA64/L2-normalize. All 300 competition target genes removed from output landmarks; each remaining perturbed gene own landmark zeroed before aggregation.','outcome_scope':'PCA uses independent RNAi phenotypes only; no held CRISPR response or competition outcomes. z-score units are never treated as counts or LFC. Cross-assay mapping is fitted separately.'})
 print('V43 fingerprints',len(names),'targets',int(np.isin(panel,names).sum()),'panel targets',len(syms),'landmarks',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()


