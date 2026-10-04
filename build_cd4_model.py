"""Select public CD4 transfer settings on disjoint donors and fit a hybrid model.

Validation is within CD4, not a VCC score or a cell-line holdout. Anonymous
contexts use a condition-averaged response; no cell-line identity is inferred.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import build_submission as build
ROOT=Path(__file__).resolve().parent

def load_subset(directory,name):
 meta=pd.read_csv(directory/f'{name}.csv')
 with np.load(directory/f'{name}.npz',allow_pickle=False) as z:data={key:z[key].copy() for key in z.files}
 if not np.array_equal(meta.source_row,data['source_rows']):raise ValueError('Subset metadata row mismatch')
 if data['log_fc'].shape!=(len(meta),len(data['genes'])):raise ValueError('Subset axes mismatch')
 return meta,data

def eligible(meta):
 return ((meta.n_cells_target>=100)&meta.ontarget_significant.astype(bool)&
         ~meta.distal_offtarget_flag.astype(bool)&~meta.neighboring_gene_KD.astype(bool)&
         ~meta.single_guide_estimate.astype(bool)&~meta.low_target_gex.astype(bool)).to_numpy()

def align_disjoint(left,right):
 lm,ld=left;rm,rd=right
 key=['target_contrast','culture_condition']
 if lm.duplicated(key).any() or rm.duplicated(key).any():raise ValueError('Duplicate donor contrast')
 lm=lm.assign(position=np.arange(len(lm)));rm=rm.assign(position=np.arange(len(rm)))
 rows=lm.merge(rm,on=key,suffixes=('_left','_right'))
 li=rows.position_left.to_numpy();ri=rows.position_right.to_numpy()
 good=eligible(lm)[li]&eligible(rm)[ri];rows=rows.loc[good];li=li[good];ri=ri[good]
 if len(rows)<20:raise ValueError('Insufficient independent donor contrasts')
 shared=sorted(set(ld['genes'])&set(rd['genes']))
 if len(shared)!=len(set(shared)):raise ValueError('Duplicate gene mapping')
 lmap={g:i for i,g in enumerate(ld['genes'])};rmap={g:i for i,g in enumerate(rd['genes'])}
 la=np.array([lmap[g] for g in shared]);ra=np.array([rmap[g] for g in shared])
 lfc_l=ld['log_fc'][li][:,la].astype(float);lfc_r=rd['log_fc'][ri][:,ra].astype(float)
 q_l=ld['adj_p_value'][li][:,la];q_r=rd['adj_p_value'][ri][:,ra]
 valid=np.isfinite(lfc_l)&np.isfinite(lfc_r)&(np.asarray(shared)[None,:]!=rows.target_contrast_gene_name_left.to_numpy()[:,None])
 lfc_l=np.clip(np.nan_to_num(lfc_l,nan=0.,posinf=0.,neginf=0.),-3,3)
 lfc_r=np.clip(np.nan_to_num(lfc_r,nan=0.,posinf=0.,neginf=0.),-3,3)
 return rows,lfc_l,lfc_r,q_l,q_r,valid

def select_settings(directory,analysis):
 report=json.loads((directory.parent.parent/'analysis/model_audit/cd4_extraction.json').read_text())
 names=[Path(item['subset']).stem for item in report if Path(item['subset']).stem!='merged']
 if len(names)!=2 or set(names[0].split('_'))&set(names[1].split('_')):raise ValueError('Donors are not disjoint')
 rows,left,right,lq,rq,valid=align_disjoint(load_subset(directory,names[0]),load_subset(directory,names[1]))
 denominator=np.maximum(valid.sum(axis=1),1)
 baseline=float(np.mean(((left**2+right**2)*valid).sum(axis=1)/denominator))
 records=[]
 for threshold in [.01,.05,.1,.2]:
  lp=left*(lq<threshold);rp=right*(rq<threshold)
  for alpha in [0.,.3,.5,.7,1.]:
   error=float(np.mean((((alpha*lp-right)**2+(alpha*rp-left)**2)*valid).sum(axis=1)/denominator))
   pred=((lq<threshold)|(rq<threshold))&valid
   sign_agree=np.sign(left)==np.sign(right)
   records.append({'fdr_threshold':threshold,'transfer_shrinkage':alpha,'relative_lfc_mse':error/max(baseline,1e-12),
      'selected_direction_pairs':int(pred.sum()),'selected_direction_agreement':float(sign_agree[pred].mean()) if pred.any() else None,
      'contrasts':len(rows),'targets':int(rows.target_contrast_gene_name_left.nunique()),'genes':valid.shape[1]})
 best=min(records,key=lambda row:(row['relative_lfc_mse'],row['fdr_threshold'],row['transfer_shrinkage']))
 if best['transfer_shrinkage']==0:raise ValueError('Independent donors did not support nonzero transfer')
 pd.DataFrame(records).to_csv(analysis/'donor_validation.csv',index=False)
 summary={'disjoint_donor_pairs':names,'selected':best,'criterion':'Mean bidirectional clipped log2FC squared error per eligible target-condition, excluding own target; zero-response anchor.',
          'limitations':['Within-CD4 donor validation is not cross-cell-type validation.',
          'Donor subsets are enriched for larger responses and do not represent all 300 targets.',
          'Grid selection error is exploratory; no third independent validation split.']}
 build.write_json(analysis/'donor_validation.json',summary)
 print('Selected donor settings',best,flush=True)
 return best

def fit_model(root,genes,targets,analysis,k562_shrinkage=.7,noise_z=1.96):
 directory=root/'reference/cd4_probe'
 settings=select_settings(directory,analysis)
 meta,data=load_subset(directory,'merged');symbols=data['genes']
 if len(set(symbols))!=len(symbols):raise ValueError('CD4 symbols are not unique')
 # K562 supplies only targets or measured genes without an eligible CD4 estimate.
 fallback_dir=analysis/'k562_fallback';fallback_dir.mkdir(exist_ok=True)
 k562=build.reference_model(root,genes,targets,fallback_dir,k562_shrinkage,noise_z)
 effects=np.zeros((1,len(targets),len(genes)),dtype=float)
 effects[0][:,k562[3]]=k562[0][0]*k562_shrinkage
 available=k562[1].copy();accepted=eligible(meta)
 gene_index={gene:i for i,gene in enumerate(genes)}
 si=np.array([i for i,gene in enumerate(symbols) if gene in gene_index]);di=np.array([gene_index[symbols[i]] for i in si])
 rows=[]
 for ti,target in enumerate(targets):
  indices=np.flatnonzero((meta.target_contrast_gene_name.to_numpy()==target)&accepted)
  if len(indices):
   raw=data['log_fc'][indices][:,si].astype(float)
   q=data['adj_p_value'][indices][:,si]
   values=np.where(np.isfinite(raw)&np.isfinite(q)&(q<settings['fdr_threshold']),np.clip(raw,-3,3),0.)
   # Average all eligible conditions; a condition-specific hit is attenuated.
   weights=np.sqrt(meta.n_cells_target.to_numpy()[indices])
   effects[0,ti,di]=settings['transfer_shrinkage']*np.average(values,axis=0,weights=weights)
   available[0,ti]=True
  rows.append({'target_gene':target,'eligible_cd4_conditions':len(indices),'cd4_source_cells':int(meta.n_cells_target.to_numpy()[indices].sum()),
   'k562_available':bool(k562[1][0,ti]),'status':'CD4 with K562 unmeasured-gene fallback' if len(indices) else 'K562 fallback' if k562[1][0,ti] else 'target-only fallback',
   'nonzero_downstream_effects':int(np.count_nonzero(effects[0,ti,genes!=target]))})
 pd.DataFrame(rows).to_csv(analysis/'reference_target_coverage.csv',index=False)
 summary={'model':'Donor-calibrated significant CD4 condition consensus, with K562 fallback',
  'targets_with_reference':int(available.sum()),'targets_with_cd4':int(sum(row['eligible_cd4_conditions']>0 for row in rows)),
  'cd4_shared_genes':len(si),'k562_shared_genes':len(k562[3]),'union_shared_genes':len(set(di)|set(k562[3])),
  'cd4_settings':settings,'k562_fallback_shrinkage':k562_shrinkage,'noise_z_assumption':noise_z,
  'condition_weighting':'Square-root source cell counts; no anonymous context assignment',
  'quality_filters':'At least 100 cells, significant on-target, two guides, no flagged cis/distal off-target or low target expression',
  'basal_profile_available':False,'validation_scope':'Disjoint public CD4 donors only, not VCC or cross-context performance.',
  'source_subsets':json.loads((root/'analysis/model_audit/cd4_extraction.json').read_text())}
 build.write_json(analysis/'reference_model.json',summary)
 np.savez_compressed(analysis/'reference_model_arrays.npz',effects=effects.astype(np.float32),available=available,
  basal=np.zeros((1,len(genes)),dtype=np.float32),dest_positions=np.arange(len(genes)),cell_lines=np.array(['CD4/K562 hybrid']),genes=genes,targets=targets)
 return effects,available,np.zeros((1,len(genes))),np.arange(len(genes)),{'shrinkage':1.,'temperature':0.},['CD4/K562 hybrid'],summary

if __name__=='__main__':
 out=ROOT/'experiments/v2_cd4_hybrid/analysis';out.mkdir(parents=True,exist_ok=True)
 _,genes,targets=build.read_spec(ROOT);fit_model(ROOT,genes,targets,out)