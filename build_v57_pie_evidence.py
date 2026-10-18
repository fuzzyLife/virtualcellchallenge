"""Stream all saved training rows; project exact PIE donor evidence onto our targets.

Preserves the original global clip, normalization and split key. No evaluation
labels enter evidence. H1/Orion diagnostics subtract matching query responses.
"""
import sys,gc,argparse
import numpy as np
from v57_pie import ROOT,A,CACHE,R,SOURCES,MATCH,sha,read,write
sys.path.insert(0,str(R/'code/src'))
from pie.data import evidence as e
from pie.data.preprocessed import PreprocessedDir,PreprocessedMeta
from pie.data.splits import resolve_split,load_split

def dirs():
 return [PreprocessedDir(p,PreprocessedMeta.model_validate_json((p/'meta.json').read_text(encoding='utf-8'))) for p in sorted((CACHE/'data').iterdir()) if (p/'meta.json').exists()]

def main():
 assert read(A/'download_receipts.json')['passed']
 import pandas as pd
 targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy(str);nt=len(targets);ti={t:i for i,t in enumerate(targets)}
 stats=read(ROOT/'analysis/post_v17_research/PIE_replogle_xdataset_data_stats.json')['stats']
 splitpath=ROOT/'analysis/v56_pie_release/xd_train.json';assert sha(splitpath)==stats['train_json_sha256']
 split=load_split(splitpath);dd=dirs();kept=resolve_split(split,dd)
 cfg=e.EvidenceConfig(seed=0,chunk=2048,lfc_clip_percentile=95.)
 key=e.evidence_cache_key(dd,sha(splitpath),cfg,.05);assert key==stats['evidence_key'],(key,stats['evidence_key'])
 bd=[e._open_build_dir(d,kept[d.dataset]) for d in dd]
 rng=np.random.default_rng(0);sample=np.concatenate([e._sample_abs_log2fc(d,2048,rng) for d in bd])
 if len(sample)>e._LFC_SAMPLE_MAX:sample=rng.choice(sample,size=e._LFC_SAMPLE_MAX,replace=False)
 clip=float(np.percentile(sample,95)) if len(sample) else 1.;clip=clip if clip>0 else 1.
 rules=e.Rules(.05,1.,clip);print('Exact author clip',clip,'sample',len(sample),flush=True)
 genes=sorted({g for d in bd for g in d.genes});gi={g:i for i,g in enumerate(genes)};ng=len(genes)
 totals={k:0. for k in ['dp_sum','dp_sumsq','dp_cnt','lfc_sum','lfc_sumsq','lfc_cnt']}
 out=CACHE/'evidence_projection';out.mkdir(exist_ok=True)
 acc={k:np.zeros((nt,ng),dtype=np.float64 if k in e.SUM_ARRAYS else np.int64) for k in e.GLOBAL_ARRAYS}
 by={k:np.zeros((len(dd),nt,ng),np.int64) for k in e.PER_DATASET_ARRAYS}
 excluded={name:{k:np.zeros((nt,ng),dtype=np.float64 if k in e.SUM_ARRAYS else np.int64) for k in e.GLOBAL_ARRAYS} for name in MATCH}
 counts={};audit_rows=[]
 for slot,d in enumerate(bd):
  cols=np.asarray([gi[g] for g in d.genes]);counts[d.dataset]=len(d.kept)
  for start in range(0,len(d.kept),256):
   rows=d.kept[start:start+256]
   own=e._own_row_stats(fc=d.fc[rows],fdr=d.fdr[rows],delta=d.dp[rows],tested=d.tested[rows],rules=rules)
   assert np.isfinite(own['dp_sum']).all()
   for k in totals:totals[k]+=float(own[k].sum(dtype=np.float64))
   for local,r in enumerate(rows):
    ctx,target=d.donor_keys[start+local]
    if target not in ti:continue
    t=ti[target]
    for k in e.GLOBAL_ARRAYS:acc[k][t,cols]+=own[k][local]
    for k,src in e.PER_DATASET_SOURCE.items():by[k][slot,t,cols]+=own[src][local]
    names=[n for n,(ds,cc) in MATCH.items() if ds==d.dataset and (cc is None or cc==ctx)]
    if names:
     packed=e._pack_rows({k:own[k][local:local+1] for k in e.GLOBAL_ARRAYS});rounded=e._unpack_rows(*packed)
     for name in names:
      for k in e.GLOBAL_ARRAYS:excluded[name][k][t,cols]+=rounded[k][0]
    audit_rows.append({'dataset':d.dataset,'context':ctx,'target':target,'source_row':int(r),'excluded_for':names})
   if start%2048==0:print('Evidence stream',d.dataset,start,'/',len(d.kept),flush=True)
   del own
  gc.collect()
 scales=e._scales_from_totals(totals)
 for k,dt in e.GLOBAL_ARRAYS.items():np.save(out/(k+'.npy'),acc[k].astype(dt))
 for k,dt in e.PER_DATASET_ARRAYS.items():np.save(out/(k+'.npy'),by[k].astype(dt))
 for name in excluded:
  for k in e.GLOBAL_ARRAYS:np.save(out/(name+'_'+k+'.npy'),excluded[name][k])
 files={p.name:sha(p) for p in out.glob('*.npy')}
 rec={'passed':True,'quantity':'Query-target projection of genuine saved training evidence; not a full author cache','original_evidence_key':key,'train_json_sha256':sha(splitpath),'genes':genes,'targets':targets.tolist(),'datasets':[d.dataset for d in bd],'clip':clip,'scales':scales.to_dict(),'full_training_totals':totals,'training_rows':counts,'target_donor_records':audit_rows,'files':files,'global_scaling':'All training rows, all contributing genes. Float64 streaming sums; equivalent sufficient statistics up to summation roundoff.','query_exclusion':'H1 removes matching target records from all VCC25 contexts; Orion removes matching target records in its own line. CD4 has no matching training source. Context response blocks missing for novel anonymous query contexts.'}
 write(out/'projection.json',rec);write(A/'evidence_receipt.json',{k:v for k,v in rec.items() if k not in ['genes','targets','target_donor_records','files']})
 write(A/'evidence_donor_provenance.json',{'train_split_sha256':sha(splitpath),'records':audit_rows})
 print('Evidence projected',scales.to_dict(),'donors',len(audit_rows),flush=True)

def blocks(target,source,axis):
 p=CACHE/'evidence_projection';meta=read(p/'projection.json');i=meta['targets'].index(target);mp={g:j for j,g in enumerate(meta['genes'])};cols=np.array([mp[g] for g in axis])
 st={k:np.asarray(np.load(p/(k+'.npy'),mmap_mode='r')[i,cols],np.float64 if k in e.SUM_ARRAYS else np.int64) for k in e.GLOBAL_ARRAYS}
 by={k:np.asarray(np.load(p/(k+'.npy'),mmap_mode='r')[:,i,cols],np.int64) for k in e.PER_DATASET_ARRAYS}
 removed=0
 if source in MATCH:
  own={k:np.asarray(np.load(p/(source+'_'+k+'.npy'),mmap_mode='r')[i,cols]) for k in e.GLOBAL_ARRAYS}
  removed=int(own['dp_cnt'].max(initial=0));e._subtract_donors(st,by,own,meta['datasets'].index(MATCH[source][0]))
 for k in e.COUNT_ARRAYS:assert (st[k]>=0).all()
 result=e._derive_blocks(st,by,e.Scales.from_dict(meta['scales']))
 result.update({k:np.zeros((len(axis),4),np.float32) for k in e.CONTEXT_BLOCK_KEYS})
 assert all(np.isfinite(v).all() for v in result.values())
 return result,{'query_records_excluded':removed,'donor_genes':int((st['dp_cnt']>0).sum()),'maximum_donors':int(st['dp_cnt'].max(initial=0))}

if __name__=='__main__':main()
