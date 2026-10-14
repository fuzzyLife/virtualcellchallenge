"""V18 blend invariants and author missing-evidence/batch parity tests."""
import argparse,sys,json
from pathlib import Path
import numpy as np

def core():
 from v18_methods import blend,shuffled_indices
 rng=np.random.default_rng(18);base=rng.normal(size=300).astype(np.float32);base[7]=0;pie=rng.normal(size=300);prob=rng.uniform(size=300);mask=np.zeros(300,bool)
 for v in blend(base,pie,prob,mask,7).values():np.testing.assert_array_equal(v,base)
 mask[:100]=True;result=blend(base,pie,prob,mask,7);np.testing.assert_array_equal(result['confidence'][100:],base[100:]);assert all(v[7]==0 for v in result.values());np.testing.assert_allclose(abs(result['amplitude']).sum(),abs(result['confidence']).sum(),rtol=1e-6)
 zero=blend(base,pie,np.zeros(300),mask,7);np.testing.assert_array_equal(zero['confidence'],base)
 assert np.all(shuffled_indices(300)!=np.arange(300)) and len(set(shuffled_indices(300)))==300
 print('V18 fallback, own-target exclusion, amplitude match, confidence bounds and derangement passed.',flush=True)
def upstream():
 from types import SimpleNamespace as S
 root=Path(__file__).resolve().parent;sys.path.insert(0,str(root/'reference/v18_pie/code/src'))
 import torch
 from pie.data.evidence import serve_evidence,Scales,EVIDENCE_KEYS
 from pie.data.dataset import PieDataset,RowRef,SourceLookup,collate
 from pie.sources.contract import read_source
 ev=S(shards={},n_donor_datasets=1,scales=Scales(dp=1.,lfc=1.));idx=S(evidence=ev,gene_cols=np.arange(3),pert_rows={})
 blocks=serve_evidence(idx,dataset='vcc_new',context='unknown',perturbation='VIM',gene_union_ids=np.arange(3));assert set(blocks)==set(EVIDENCE_KEYS)
 for key,value in blocks.items():assert value.shape==(3,3 if key=='evidence_prov' else 4) and np.count_nonzero(value)==0
 source={name:read_source(root/'reference/v18_pie/sources'/('ncbi_subset' if name=='ncbi_text' else name)) for name in ['esm2','ncbi_text','string_space','depmap_gene_effect','context_text','perturbation_text']};lookup=SourceLookup(source,{})
 d=S(dataset='vcc_new',genes=['VIM','ABCD1','ACLY'],ctrl_means=np.array([[1.,0.,.3]],np.float32),meta=S(context_to_id={'unknown':0}))
 dataset=PieDataset([d],[RowRef(0,-1,'unknown','VIM')],[np.arange(3)],lookup,idx,.05,'none');sample=dataset[0];assert 'fold_changes' not in sample and 'lfc_true' not in sample
 expected={name:torch.from_numpy(v.astype(np.float32)) for name in lookup.dims if (v:=lookup.tokens(name,context='unknown',perturbation='VIM')) is not None};assert 'context_text' not in expected
 batch=collate([sample],lookup.dims)
 for name,value in expected.items():np.testing.assert_array_equal(batch.source_tokens[name].numpy()[0],value.numpy());assert batch.source_masks[name].all()
 for key,value in blocks.items():np.testing.assert_array_equal(batch.groups[0].evidence[key].numpy()[0],value)
 np.testing.assert_array_equal(batch.groups[0].ctrl_means.numpy()[0],d.ctrl_means[0]);print('Author serve_evidence missing-record semantics, controls-only PieDataset and collate parity passed.',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--upstream',action='store_true');args=p.parse_args()
 if args.upstream:upstream()
 else:core()
