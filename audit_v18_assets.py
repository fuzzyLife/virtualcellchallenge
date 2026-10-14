"""Independently verify compact target sources, pinned author code and all recorded assets."""
from pathlib import Path
import hashlib,json
import numpy as np
import build_submission as b
ROOT=b.ROOT;R=ROOT/'reference/v18_pie';A=ROOT/'experiments/v18_pie_expert/analysis'
def main():
 assets=json.loads((A/'assets.json').read_text())
 for item in assets['files']:
  p=ROOT/item['path'];assert p.stat().st_size==item['bytes'] and b.sha256(p)==item['sha256']
 records=json.loads((A/'ncbi_ranges.json').read_text());meta=json.loads((R/'sources/ncbi_subset/meta.json').read_text());values=np.load(R/'sources/ncbi_subset/embeddings.npy',mmap_mode='r',allow_pickle=False);offset=np.load(R/'sources/ncbi_subset/offsets.npy',allow_pickle=False)
 assert meta['keys']==[v['key'] for v in records['ranges']]
 for i,item in enumerate(records['ranges']):
  raw=values[offset[i]:offset[i+1]].tobytes();assert len(raw)==item['hi']-item['lo']+1 and hashlib.sha256(raw).hexdigest()==item['sha256'];assert int(offset[i+1]-offset[i])==item['tokens']
 tree=json.loads((ROOT/'analysis/post_v17_research/pie_tree.json').read_text())['tree'];code=0
 for item in tree:
  if item['type']=='blob' and item['path'].startswith('src/pie/'):
   raw=(R/'code'/item['path']).read_bytes();assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==item['sha'];code+=1
 assert b.sha256(R/'train.json')=='6143b4d1bdc30f09a323a5af944fe9049c77c724909fb20fdde0039b7faeddfb'
 b.write_json(A/'asset_independent_audit.json',{'passed':True,'recorded_files':len(assets['files']),'original_author_git_blobs':code,'complete_target_token_sequences':len(meta['keys']),'compact_ncbi_payload_bytes':records['total_payload_bytes'],'all_compact_slices_match_original_range_receipts':True,'canonical_split_hash_verified':True,'asset_manifest_sha256':b.sha256(A/'assets.json'),'auditor_sha256':b.sha256(Path(__file__)),'limitation':'The full 43GB NCBI array was not downloaded or globally rehashed. Complete selected ranges fetched from pinned HTTPS asset; range hashes retained and compact reconstruction verified.'});print('Independent PIE source/compact-range/asset audit passed',flush=True)
if __name__=='__main__':main()
