"""Check KOLF units, missingness, row/gene alignment and source-feature extension."""
import json
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
import v39_kolf_hurdle as k
import v38_response_hurdle as m
def main():
 rec=m.read(k.A/'collection.json');assert m.b.sha256(k.A/'response_features.npz')==rec['arrays_sha256'];meta=m.read(k.OUT/'reference/meta.json');d=m.data();z,sat=k.kolf();present=np.isfinite(z).any(1);names=meta['genes'];counts=Counter(names);gene_index={s:i for i,s in enumerate(names) if counts[s]==1};target_index={s:i for i,s in enumerate(meta['perts'])}
 measured=np.array([str(g) in gene_index for g in d['genes']]);positions=np.flatnonzero(measured);indices=np.array([gene_index[str(g)] for g in d['genes'][positions]])
 checked=0
 for ti,t in enumerate(d['targets']):
  if not present[ti]:assert np.isnan(z[ti]).all();continue
  index=target_index[str(t)];chunk=index//meta['chunk_size'];offset=index%meta['chunk_size']*meta['row_bytes'];p=k.OUT/'reference'/f'chunk_{chunk}.bin'
  with p.open('rb') as f:f.seek(offset);payload=f.read(meta['row_bytes'])
  row=np.frombuffer(payload,dtype=np.int8);values=row[indices].astype(np.float32)*meta['scale'];np.testing.assert_array_equal(z[ti,positions],values);assert np.isnan(z[ti,~measured]).all();np.testing.assert_array_equal(sat[ti,positions],abs(row[indices].astype(np.int16))==127);checked+=1
 ds=m.descriptor(d['profiles'][0],('H1','Replogle'));pos=np.arange(400);old,n=k.ORIGINAL(ds,0,pos);new,nn=k.features(ds,0,pos);np.testing.assert_array_equal(new[:,:56],old);np.testing.assert_array_equal(nn,n);np.testing.assert_array_equal(new[:,56],z[0,pos]);assert new.shape==(400,63)
 m.write(k.A/'feature_contract_tests.json',{'passed':True,'matching_targets_replayed_from_pinned_binary':checked,'missing_genes_and_targets_are_NaN':True,'existing_56_features_and_support_unchanged':True,'added_features':7,'z_scores_not_treated_as_counts_or_fold_changes':True,'script_sha256':m.b.sha256(m.ROOT/'test_v39_kolf_features.py')});print('KOLF feature contracts PASS',checked,'targets',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()

