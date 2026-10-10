"""Verified, locally staged RDS extraction and deduplicated full-axis Flex source cache."""
import argparse,subprocess,tempfile,hashlib,json,shutil,gc
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from collect_v12_resources import digest
import build_submission as build
ROOT=Path(__file__).resolve().parent;D=ROOT/'reference/v12_resources';A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def extract(name):
 r=json.loads((ROOT/'analysis/post_v11_research/data/viperturb_record.json').read_text());f=next(f for f in r['files'] if f['key']==name)
 source=D/'viperturb'/name
 if not source.exists():source=D/'inspected_bins'/name
 assert source.stat().st_size==f['size']
 out=D/'flex_extract'/name.replace('.RDS','');out.mkdir(parents=True,exist_ok=True)
 if (out/'extraction_summary.txt').exists():return
 # Stage once and verify publisher checksum DURING the copy. Retain a local
 # temp copy after an extraction failure, allowing reproducible repair/retry.
 cache=Path(tempfile.gettempdir())/'v12_flex_rds_cache';cache.mkdir(exist_ok=True);local=cache/name
 if local.exists():assert local.stat().st_size==f['size'] and 'md5:'+digest(local,'md5')==f['checksum']
 else:
  pending=local.with_suffix('.partial');md5=hashlib.md5()
  with source.open('rb') as src,pending.open('wb') as dst:
   while chunk:=src.read(8*1024**2):dst.write(chunk);md5.update(chunk)
  assert pending.stat().st_size==f['size'] and 'md5:'+md5.hexdigest()==f['checksum'];pending.replace(local)
 cmd=[str(D/'tools/R/bin/x64/Rscript.exe'),'--vanilla',str(ROOT/'extract_v12_flex.R'),str(local),str(out),str(ROOT/'pert_counts.csv'),name,str(D/'guest_list_protospacer_reference.csv')]
 subprocess.run(cmd,check=True)
 # Only the exact verified file created in our named temp cache is removed.
 assert local.parent.resolve()==cache.resolve() and local.name==name;local.unlink()
 print('Extracted verified bin',name,flush=True)

def combine():
 from fit_v7_response import read_cache,summarize_raw
 d=read_cache();genes=d['genes'];targets=d['targets'];blocks=[];frames=[];manifest=[];assay=None
 for name in ['genome_wide_binA','genome_wide_binB','genome_wide_binC']:
  folder=D/'flex_extract'/name;meta=pd.read_csv(folder/'selected_cells.csv');sg=pd.read_csv(folder/'gene_names.csv').gene_name.to_numpy(dtype=str)
  shape=list(map(int,(folder/'shape.txt').read_text().splitlines()));v=np.fromfile(folder/'counts_data.bin',dtype='<i4');idx=np.fromfile(folder/'counts_indices.bin',dtype='<i4');ptr=np.fromfile(folder/'counts_indptr.bin',dtype='<i4')
  assert len(v)==len(idx)==shape[2] and len(ptr)==shape[0]+1 and ptr[-1]==len(v)
  x=sparse.csr_matrix((v,idx,ptr),shape=tuple(shape[:2]));np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(),meta.total_counts)
  lookup={g:i for i,g in enumerate(genes)};dest=np.array([lookup.get(g,-1) for g in sg]);coo=x.tocoo();keep=dest[coo.col]>=0
  projected=sparse.csr_matrix((coo.data[keep],(coo.row[keep],dest[coo.col[keep]])),shape=(len(meta),len(genes)))
  if assay is None:assay=np.isin(genes,sg)
  else:np.testing.assert_array_equal(assay,np.isin(genes,sg))
  blocks.append(projected);frames.append(meta);manifest.append({'bin':name,'cells':len(meta),'genes':len(sg),'measured_challenge_genes':int(np.isin(genes,sg).sum()),'metadata_sha256':digest(folder/'selected_cells.csv')});del x,coo;gc.collect()
 x=sparse.vstack(blocks,format='csr');meta=pd.concat(frames,ignore_index=True);seen={};keep=[];duplicates=0
 for i,cell in enumerate(meta.cell_id):
  if cell in seen:
   first=seen[cell];assert meta.iloc[i].gene_target==meta.iloc[first].gene_target=='Non-Targeting','Duplicate perturbed cell across bins'
   assert meta.iloc[i]['sample']==meta.iloc[first]['sample'] and meta.iloc[i].total_counts==meta.iloc[first].total_counts
   assert (x[i]-x[first]).nnz==0,'Duplicated NTC counts differ';duplicates+=1
  else:seen[cell]=i;keep.append(i)
 x=x[keep];meta=meta.iloc[keep].reset_index(drop=True);meta['partition']=0
 meta['aliquot_sample']=meta['sample'];meta['sample']=meta['aliquot_sample'].str.extract(r'(L[12])$',expand=False)
 assert meta['sample'].notna().all() and set(meta['sample'])=={'L1','L2'}
 assert meta[meta.gene_target.eq('Non-Targeting')].groupby('sample').size().min()>=20
 assert meta.cell_id.is_unique and set(targets)<=set(meta.gene_target)
 counts=meta.gene_target.value_counts();summary=summarize_raw(x,meta,genes,targets,assay,2)
 own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);quality=(summary['knockdown']>=.25)&(summary['control'][np.arange(len(targets)),own]>5)
 summary['raw_valid']=summary['valid'].copy();summary['valid'] &= quality[:,None];summary['q'][~summary['valid']]=1;summary['target_qc_pass']=quality
 out=D/'flex_source';out.mkdir(exist_ok=True);sparse.save_npz(out/'selected_counts.npz',x);meta.to_csv(out/'selected_cells.csv',index=False)
 np.savez_compressed(out/'summaries.npz',**summary)
 pd.DataFrame({'target_gene':targets,'cells':[int(counts.get(t,0)) for t in targets],'own_measured':[bool(assay[np.flatnonzero(genes==t)[0]]) for t in targets], 'usable_outputs':summary['valid'].sum(axis=1),'target_qc_pass':quality,'knockdown':summary['knockdown']}).to_csv(A/'flex_target_coverage.csv',index=False)
 record={'doi':'10.5281/zenodo.18460279','license':'CC-BY-4.0','bins':manifest,'duplicated_controls_removed':duplicates,'unique_cells':len(meta),'unique_controls':int(meta.gene_target.eq('Non-Targeting').sum()),'all_300_targets_present':True,'measured_challenge_genes':int(assay.sum()),'summary_sha256':digest(out/'summaries.npz'), 'counts_sha256':digest(out/'selected_counts.npz'),'metadata_sha256':digest(out/'selected_cells.csv'),'script_sha256':digest(Path(__file__)),'r_script_sha256':digest(ROOT/'extract_v12_flex.R'),'preprocessing_plan_sha256':digest(A/'flex_source_preprocessing_plan.json'),'qc_plan_sha256':digest(A/'flex_source_qc_plan.json'),'control_identity_correction_sha256':digest(A/'flex_control_identity_correction.json'),'normalization':'Actual full RNA-assay library before projection; matched per-lane NTCs; original aliquots retained; normal-q summaries are proxies only.'}
 build.write_json(A/'flex_extraction.json',record);print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--bin');p.add_argument('--combine',action='store_true');a=p.parse_args()
 if a.bin:extract(a.bin)
 if a.combine:combine()
