"""Compare cached KD with author's construct metadata; do not change frozen models."""
import json
from pathlib import Path
import h5py,numpy as np,pandas as pd
from v13_methods import ROOT,A
from fit_v7_response import read_cache
import build_submission as build

def main():
 data=read_cache();file=ROOT/'reference/K562_gwps_raw_bulk_01.h5ad'
 assert build.sha256(file)=='7cec96b3b76169abbf6b6ab9d10bf00d71d942d89e63292351f745e130b154db'
 with h5py.File(file,'r') as f:
  labels=f['obs/gene_transcript'].asstr()[:];cells=f['obs/num_cells_filtered'][:].astype(float);fold=f['obs/fold_expr'][:].astype(float)
 names=np.array([s.split('_',2)[1] if len(s.split('_',2))==3 and s.split('_',1)[0].isdigit() else '' for s in labels])
 rows=[]
 for i in np.flatnonzero(data['row_families']=='K562'):
  target=data['targets'][data['row_targets'][i]];selected=(names==target)&(cells>=30)&np.isfinite(fold)&(fold<=.5)
  assert selected.any();published=float(np.clip(1-np.average(fold[selected],weights=cells[selected]),0,.99))
  own=int(np.flatnonzero(data['genes']==target)[0]);control=float(data['profiles'][0,own])
  rows.append({'target':target,'cached_KD':float(data['knockdown'][i]),'published_construct_weighted_KD':published,
   'own_control_CPM':control,'target_quantifiable':bool(np.isfinite(control)&(control>5)),
   'source_constructs':int(selected.sum()),'source_cells':float(cells[selected].sum())})
 frame=pd.DataFrame(rows);frame.to_csv(A/'published_k562_knockdown_audit.csv',index=False)
 missing=frame[~frame.target_quantifiable]
 build.write_json(A/'published_knockdown_audit.json',{'source_sha256':build.sha256(file),'rows':len(frame),'unquantifiable_targets':len(missing),
  'published_missing_KD_min':float(missing.published_construct_weighted_KD.min()),'published_missing_KD_median':float(missing.published_construct_weighted_KD.median()),
  'frozen_quality_stage_uses':'Neutral unknown-row KD=.9; does not silently substitute new metadata after plan freezing.',
  'future_test':'Replace neutral unknown-row assumption with matched author construct-weighted KD in a separately frozen comparison.',
  'script_sha256':build.sha256(Path(__file__))})
 print('Author KD available for',len(missing),'unquantifiable targets; median',missing.published_construct_weighted_KD.median(),flush=True)
if __name__=='__main__':main()
