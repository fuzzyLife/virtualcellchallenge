"""Extract valid LINCS RNAi consensus signatures; never reinterpret z as counts."""
import json,gzip
from pathlib import Path
from collections import Counter
import h5py,numpy as np,pandas as pd
import build_submission as b
ROOT=b.ROOT;OUT=ROOT/'experiments/v43_lincs_response';A=OUT/'analysis';R=OUT/'reference'
def main():
 rec=json.loads((A/'collection.json').read_text(encoding='utf-8'))
 for r in rec['files']:assert b.sha256(R/r['file'])==r['sha256']
 meta=pd.read_csv(R/'GSE106127_CGS_PRIME_meta.txt.gz',sep='\t');gene=pd.read_csv(R/'GSE106127_gene_info.txt.gz',sep='\t');_,genes,targets=b.read_spec(ROOT)
 path=R/'GSE106127_CGS_PRIME_n33839x978.gctx'
 with h5py.File(path,'r') as f:
  ids=f['0/META/COL/id'].asstr()[:];gids=f['0/META/ROW/id'].asstr()[:];values=f['0/DATA/0/matrix'][:]
 assert values.shape==(len(ids),len(gids)) and np.isfinite(values).all();assert len(set(ids))==len(ids)
 meta=meta.set_index('sig_id').loc[ids].reset_index();gm=gene.set_index(gene.pr_gene_id.astype(str));symbols=gm.loc[gids,'pr_gene_symbol'].astype(str).to_numpy();counts=Counter(symbols);lookup={g:i for i,g in enumerate(genes)};positions=np.array([lookup.get(s,-1) if counts[s]==1 else -1 for s in symbols]);keep=positions>=0;pos=positions[keep];values=values[:,keep].astype(np.float32)
 eligible=(meta.pert_type=='trt_sh.cgs')&(meta.distil_nsample>=3)&np.isfinite(meta.distil_cc_q75)&(meta.distil_cc_q75>0)
 contexts=sorted(meta.loc[eligible,'cell_id'].unique());tp={t:i for i,t in enumerate(targets)};z=np.full((len(contexts),len(targets),len(pos)),np.nan,np.float32);quality=np.zeros((len(contexts),len(targets)),np.float32);source=[];centers=[]
 for ci,c in enumerate(contexts):
  allids=np.flatnonzero(eligible&(meta.cell_id==c));w=np.clip(meta.distil_cc_q75.to_numpy(float)[allids],0,1)**2
  center=np.average(values[allids],axis=0,weights=w);centers.append(center)
  for i in allids:
   t=str(meta.pert_iname.iloc[i])
   if t not in tp:continue
   ti=tp[t];assert np.isnan(z[ci,ti]).all(),'Duplicate gene/cell signature requires explicit aggregation'
   z[ci,ti]=np.clip(values[i]-center,-8,8);quality[ci,ti]=float(np.clip(meta.distil_cc_q75.iloc[i],0,1)**2)
   source.append({'signature_id':str(ids[i]),'matrix_row':int(i),'cell':c,'target':t,'hairpins':int(meta.distil_nsample.iloc[i]),'sister_cc_q75':float(meta.distil_cc_q75.iloc[i]),'time':str(meta.pert_itime.iloc[i])})
 mask=np.isfinite(z);w=quality[:,:,None]*mask;den=w.sum(0);mean=np.divide(np.where(mask,z,0)*w,den[None,:,:],out=np.zeros_like(z),where=den[None,:,:]>0).sum(0)
 agreement=np.divide(abs((np.where(mask,np.sign(z),0)*w).sum(0)),den,out=np.zeros_like(mean),where=den>0);n=mask.sum(0);signal=np.clip(mean,-5,5)*agreement;signal[n<3]=np.nan
 for ti,t in enumerate(targets):signal[ti,genes[pos]==t]=np.nan
 file=A/'response_features.npz';assert not file.exists();np.savez_compressed(file,signal=signal.astype(np.float32),mean_z=mean,agreement=agreement,n_contexts=n,z_by_context=z,quality=quality,centers=np.array(centers),positions=pos,genes=np.asarray(genes,str),targets=np.asarray(targets,str),contexts=np.array(contexts,str))
 b.write_json(A/'extraction.json',{'arrays_sha256':b.sha256(file),'collection_sha256':b.sha256(A/'collection.json'),'script_sha256':b.sha256(Path(__file__)),'contexts':contexts,'matching_targets_with_3_contexts':int(np.any(np.isfinite(signal),axis=1).sum()),'output_landmark_overlap':len(pos),'RNAi_consensus_signatures':len(source),'rows':source,'units':'PRIME-moderated z scores. Center each landmark within source cell line using all eligible RNAi interventions, then quality-weighted across contexts and downweight sign disagreement. No conversion to LFC or counts.','quality':'RNAi trt_sh.cgs only, at least3 hairpin signatures and positive sister consistency. Overexpression/vector controls excluded. At least3 source cell contexts for support; different lines not treated as independent donor replicates.','gene_mapping':'Exact unique symbol match only; no guessed aliases or imputed output genes. Own target omitted from transferable features.','PRIME_note':'Source publication discusses PC1 removal for CGS; this does not eliminate shRNA seed effects. Cross-assay calibration must be validated.'})
 print('LINCS extract',len(source),'RNAi signatures',len(pos),'landmarks',int(np.any(np.isfinite(signal),axis=1).sum()),'targets',flush=True)
if __name__=='__main__':main()

