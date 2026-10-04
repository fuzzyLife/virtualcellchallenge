"""Combine the competition controls and specifications into one lossless CSV.

Default layout: genes as rows; context-prefixed control cells as numeric columns.
Cell IDs, non-targeting guide IDs and technical totals are encoded in headers.
Only the supplied competition controls are exported; there are no observed
knockdown responses. The cache separates HDF5 extraction from streaming CSV
writing, permitting a writer with only numpy and pandas installed.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path
import time
from urllib.parse import quote
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent

def log(message):print(time.strftime('%H:%M:%S')+' '+message,flush=True)
def digest(path):
 value=hashlib.sha256()
 with path.open('rb') as file:
  for block in iter(lambda:file.read(8*1024**2),b''):value.update(block)
 return value.hexdigest()
def save_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
def h5_column(group):
 import h5py
 if isinstance(group,h5py.Group):
  if 'categories' in group and 'codes' in group:
   categories=h5_column(group['categories']);codes=group['codes'][:]
   if (codes<0).any():raise ValueError('Missing cell metadata is not expected')
   return categories[codes]
  if 'values' in group and 'mask' in group:
   if group['mask'][:].any():raise ValueError('Missing identifiers are not expected')
   return h5_column(group['values'])
  raise ValueError('Unsupported metadata encoding')
 return group.asstr()[:] if group.dtype.kind in 'OSU' else group[:]
def extract(root,cache):
 import h5py
 from scipy import sparse
 manifest=json.loads((root/'manifest.json').read_text());genes=pd.read_csv(root/'gene_names.csv').gene_name.astype(str).tolist();targets=pd.read_csv(root/'pert_counts.csv').target_gene.astype(str).tolist()
 if len(genes)!=manifest['n_genes'] or len(set(genes))!=len(genes):raise ValueError('Gene specification mismatch')
 if len(targets)!=manifest['n_constructs'] or len(set(targets))!=len(targets) or not set(targets)<=set(genes):raise ValueError('Target specification mismatch')
 cache.mkdir(parents=True,exist_ok=True)
 info={'manifest':manifest,'genes':genes,'targets':targets,'contexts':{},'source_files':{}}
 for name in ['manifest.json','gene_names.csv','pert_counts.csv']:
  info['source_files'][name]={'bytes':(root/name).stat().st_size,'sha256':digest(root/name)}
 for context in manifest['contexts']:
  log(f'Extracting all {context} cells without filtering or normalization')
  path=root/f'context_{context}.h5ad'
  with h5py.File(path,'r') as f:
   axis=h5_column(f['var'][f['var'].attrs['_index']]).astype(str)
   if not np.array_equal(axis,genes):raise ValueError(f'{context} gene order mismatch')
   if set(f['obs'].keys())!={'_index','context','target_gene','ntc_id'}:raise ValueError('Unexpected observation fields; update export to preserve them')
   ids=h5_column(f['obs'][f['obs'].attrs['_index']]).astype(str)
   labels=h5_column(f['obs/context']).astype(str);pert=h5_column(f['obs/target_gene']).astype(str);guides=h5_column(f['obs/ntc_id']).astype(str)
   if set(labels)!={context} or set(pert)!={manifest['control_label']} or len(set(ids))!=len(ids):raise ValueError('Control metadata mismatch')
   group=f['X'];shape=tuple(group.attrs['shape'])
   if group.attrs['encoding-type']!='csr_matrix' or shape!=(manifest['per_context'][context]['control_cells'],len(genes)):raise ValueError('Unexpected count matrix encoding or shape')
   values=group['data'][:];indices=group['indices'][:];indptr=group['indptr'][:]
  for start in range(0,len(values),4_000_000):
   part=values[start:start+4_000_000]
   if not np.isfinite(part).all() or (part<0).any() or (part!=np.floor(part)).any() or (part>np.iinfo(np.int32).max).any():raise ValueError('Expected nonnegative whole raw counts')
  x=sparse.csr_matrix((values.astype(np.int32),indices,indptr),shape=shape)
  x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
  libraries=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel();detected=np.diff(x.indptr)
  if (libraries<=0).any():raise ValueError('Empty control cell')
  guide_names,guide_counts=np.unique(guides,return_counts=True)
  if len(guide_names)!=manifest['per_context'][context]['n_ntc_ids'] or not (guide_counts==manifest['cells_per_pert']).all():raise ValueError('Control guide counts mismatch')
  gene_totals=np.asarray(x.sum(axis=0,dtype=np.int64)).ravel()
  matrix=x.tocsc();matrix.sort_indices()
  for key,value in [('data',matrix.data),('indices',matrix.indices),('indptr',matrix.indptr)]:np.save(cache/f'{context}_{key}.npy',value,allow_pickle=False)
  info['contexts'][context]={'cell_ids':ids.tolist(),'ntc_ids':guides.tolist(),'cell_total_counts':libraries.tolist(),
   'cell_detected_genes':detected.tolist(),'gene_total_counts':gene_totals.tolist(),'gene_detected_cells':np.diff(matrix.indptr).tolist(),
   'shape':[int(value) for value in shape],'nnz':int(matrix.nnz),'source_file':path.name}
  info['source_files'][path.name]={'bytes':path.stat().st_size,'sha256':digest(path)}
  del x,matrix,values,indices,indptr
  gc.collect()
 save_json(cache/'metadata.json',info)
 log('Source extraction complete')
 return info

def specifications(info):
 m=info['manifest'];genes=info['genes'];n=len(genes);target_order={gene:i+1 for i,gene in enumerate(info['targets'])}
 target_flag=np.array([gene in target_order for gene in genes],dtype=np.int32)
 annotations={'gene_name':genes,'gene_index_1based':np.arange(1,n+1),
  'is_requested_knockdown_target':target_flag,'requested_target_index_1based':[target_order.get(gene,0) for gene in genes],
  'requested_cells_per_target_context':target_flag*m['cells_per_pert'],'season':m['season'],'partition':m['partition'],'panel_id':m['panel_id'],
  'source_url':'https://virtualcellchallenge.org/downloads','count_unit':'raw_gene_count','observed_target_gene':m['control_label'],
  'has_observed_knockdown_responses':0,'panel_gene_count':m['n_genes'],'panel_target_count':m['n_constructs'],
  'panel_required_cells_per_target':m['cells_per_pert'],'panel_required_prediction_cells_total':m['n_constructs']*m['cells_per_pert']*len(m['contexts'])}
 headers=[]
 for context in m['contexts']:
  entry=info['contexts'][context];per=m['per_context'][context]
  annotations[f'context_{context}_control_cells']=per['control_cells'];annotations[f'context_{context}_control_guides']=per['n_ntc_ids']
  annotations[f'context_{context}_ground_truth_cells_declared']=per['ground_truth_cells']
  annotations[f'context_{context}_gene_total_raw_count']=entry['gene_total_counts'];annotations[f'context_{context}_gene_detected_control_cells']=entry['gene_detected_cells']
  for ci,(cell,guide,library,detected) in enumerate(zip(entry['cell_ids'],entry['ntc_ids'],entry['cell_total_counts'],entry['cell_detected_genes']),1):
   headers.append(f'context_{context}_col{ci:05d}__cell_id={quote(cell,safe="")}__ntc_id={quote(guide,safe="")}__total_counts={library}__detected_genes={detected}')
 if len(set(headers))!=len(headers):raise ValueError('Duplicate expression column names')
 return pd.DataFrame(annotations),headers

def load_arrays(info,cache):
 return {context:{name:np.load(cache/f'{context}_{name}.npy',mmap_mode='r',allow_pickle=False) for name in ['data','indices','indptr']} for context in info['manifest']['contexts']}

def write_csv(root,cache,output,force,batch_genes):
 if output.exists() and not force:raise FileExistsError('CSV exists; use --force for an intentional replacement')
 info=json.loads((cache/'metadata.json').read_text());annotations,headers=specifications(info);arrays=load_arrays(info,cache)
 total_cells=len(headers);contexts=info['manifest']['contexts'];partial=output.with_suffix('.csv.partial');numeric_hash=hashlib.sha256()
 output.parent.mkdir(parents=True,exist_ok=True)
 with partial.open('w',encoding='utf-8',newline='',buffering=8*1024**2) as file:
  writer=csv.writer(file,lineterminator='\n')
  writer.writerow(annotations.columns.tolist()+headers)
  for start in range(0,len(annotations),batch_genes):
   stop=min(start+batch_genes,len(annotations));block=np.zeros((stop-start,total_cells),dtype=np.int32);offset=0
   for context in contexts:
    arr=arrays[context]
    for gene in range(start,stop):
     lo,hi=arr['indptr'][gene:gene+2];block[gene-start,offset+arr['indices'][lo:hi]]=arr['data'][lo:hi]
    offset+=info['contexts'][context]['shape'][0]
   numeric_hash.update(block.astype('<i4',copy=False).tobytes(order='C'))
   prefix=annotations.iloc[start:stop].astype(str).to_numpy()
   for metadata_row,count_row in zip(prefix,block):
    writer.writerow(metadata_row.tolist()+count_row.tolist())
   if stop%512<batch_genes or stop==len(annotations):log(f'Wrote {stop:,}/{len(annotations):,} genes; {total_cells:,} cell columns')
 partial.replace(output)
 report={'created_utc':datetime.now(timezone.utc).isoformat(),'file':str(output),'bytes':output.stat().st_size,'layout':'genes_as_rows_cells_as_columns',
  'gene_rows':len(annotations),'cell_columns':len(headers),'annotation_columns':annotations.columns.tolist(),'total_columns':len(headers)+len(annotations.columns),
  'requested_target_genes':len(info['targets']),'provided_control_cells':{c:info['contexts'][c]['shape'][0] for c in contexts},
  'numeric_gene_major_int32_sha256':numeric_hash.hexdigest(),'csv_sha256':digest(output),'source_files':info['source_files'],
  'scope':'All provided competition control count matrices and metadata; no observed knockdown-response labels, external references or generated predictions.',
  'count_transformation':'None: raw counts retained exactly, implicit sparse zeros serialized as zeros.',
  'validation':{'completed':False}}
 save_json(root/'analysis/combined_csv_export.json',report)
 log(f'CSV written: {output} ({report["bytes"]/1e9:.3f} GB); verifying next')
 return report

def verify_csv(root,cache,output):
 info=json.loads((cache/'metadata.json').read_text());annotations,headers=specifications(info)
 report_path=root/'analysis/combined_csv_export.json';report=json.loads(report_path.read_text());count_offset=len(annotations.columns)
 expected_annotation=annotations.astype(str).to_numpy();numeric_hash=hashlib.sha256();cell_totals=np.zeros(len(headers),dtype=np.int64);nonzeros=0;gene_count=0
 ranges=[];offset=0
 for context in info['manifest']['contexts']:
  size=info['contexts'][context]['shape'][0];ranges.append((context,offset,offset+size));offset+=size
 with output.open('r',encoding='utf-8',newline='') as file:
  reader=csv.reader(file)
  if next(reader)!=annotations.columns.tolist()+headers:raise ValueError('CSV header/metadata alignment mismatch')
  for gi,row in enumerate(reader):
   if gi>=len(annotations) or len(row)!=count_offset+len(headers) or row[:count_offset]!=expected_annotation[gi].tolist():raise ValueError(f'CSV row/annotation mismatch at gene {gi}')
   counts=np.fromstring(','.join(row[count_offset:]),sep=',',dtype=np.int32)
   if len(counts)!=len(headers) or (counts<0).any():raise ValueError('CSV count parsing mismatch')
   numeric_hash.update(counts.astype('<i4',copy=False).tobytes());cell_totals+=counts;nonzeros+=int(np.count_nonzero(counts))
   for context,lo,hi in ranges:
    if int(counts[lo:hi].sum(dtype=np.int64))!=info['contexts'][context]['gene_total_counts'][gi]:raise ValueError('Per-context gene sum mismatch')
   gene_count=gi+1
   if gene_count%2048==0:log(f'Verified {gene_count:,}/{len(annotations):,} gene rows')
 if gene_count!=len(annotations) or numeric_hash.hexdigest()!=report['numeric_gene_major_int32_sha256']:raise ValueError('Full numeric roundtrip mismatch')
 if nonzeros!=sum(info['contexts'][c]['nnz'] for c in info['manifest']['contexts']):raise ValueError('Nonzero count mismatch')
 for context,lo,hi in ranges:
  if not np.array_equal(cell_totals[lo:hi],info['contexts'][context]['cell_total_counts']):raise ValueError('Per-cell count totals mismatch')
 if digest(output)!=report['csv_sha256']:raise ValueError('CSV checksum changed')
 report['validation']={'completed':True,'every_numeric_value_roundtripped':True,'every_gene_and_cell_identifier_checked':True,
  'all_cell_totals_match':True,'all_context_gene_totals_match':True,'nonzeros':nonzeros,'gene_rows':gene_count,'target_flags':int(annotations.is_requested_knockdown_target.sum())}
 save_json(report_path,report);log('READY: all cells, genes, target flags, counts and identifiers verified')

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=ROOT);parser.add_argument('--output',type=Path)
 parser.add_argument('--phase',choices=['extract','write','verify','all'],default='all');parser.add_argument('--batch-genes',type=int,default=64);parser.add_argument('--force',action='store_true');args=parser.parse_args()
 root=args.root.resolve();output=args.output.resolve() if args.output else root/'competition_controls_combined.csv';cache=root/'analysis/competition_csv_cache'
 if args.batch_genes<1:parser.error('--batch-genes must be positive')
 originals={root/name for name in ['manifest.json','gene_names.csv','pert_counts.csv','context_A.h5ad','context_B.h5ad','context_C.h5ad','prediction.h5ad','prediction.vcc']}
 if output in originals:parser.error('Export cannot replace an original input or submission')
 if args.phase in ['all','write'] and output.exists() and not args.force:parser.error('CSV exists; use --force for an intentional replacement')
 if args.phase in ['all','extract']:extract(root,cache)
 if args.phase in ['all','write']:write_csv(root,cache,output,args.force,args.batch_genes)
 if args.phase in ['all','verify']:verify_csv(root,cache,output)
if __name__=='__main__':main()