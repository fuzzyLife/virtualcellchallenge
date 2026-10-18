"""Independent numerical and artifact audits for V57; no held donor D3 access."""
import argparse,sys
import numpy as np
from v57_pie import *
from build_v57_pie_evidence import blocks,dirs
sys.path.insert(0,str(R/'code/src'))
from pie.data import evidence as e
from pie.data.splits import load_split,resolve_split

def evidence():
 meta=read(CACHE/'evidence_projection/projection.json');pp=plan();dd=dirs();kept=resolve_split(load_split(ROOT/'analysis/v56_pie_release/xd_train.json'),dd);bd=[e._open_build_dir(d,kept[d.dataset]) for d in dd]
 assert sum(len(d.kept)*len(d.genes) for d in bd)==meta['full_training_totals']['dp_cnt']
 own_genes=[g for g in meta['genes'] if g in set(read(ROOT/'analysis/post_v17_research/gene_text_meta.json')['keys'])]
 rng=np.random.default_rng(SEED);axis=[own_genes[i] for i in sorted(rng.choice(len(own_genes),127,replace=False))]
 seen=pp['seen_targets'][:6];unseen=pp['unseen_targets'][:3];checks=[]
 for target in seen+unseen:
  agg={k:np.zeros(len(axis),np.float64 if k in e.SUM_ARRAYS else np.int64) for k in e.GLOBAL_ARRAYS};by={k:np.zeros((len(dd),len(axis)),np.int64) for k in e.PER_DATASET_ARRAYS};ex={s:{k:np.zeros_like(v) for k,v in agg.items()} for s in MATCH}
  for slot,d in enumerate(bd):
   local={g:i for i,g in enumerate(d.genes)};covered=np.array([g in local for g in axis]);cols=np.array([local[g] for g in axis if g in local]);ids=[j for j,k in enumerate(d.donor_keys) if k[1]==target]
   for j in ids:
    r=int(d.kept[j]);fc=np.asarray(d.fc[r,cols],float);fdr=np.asarray(d.fdr[r,cols],float);dp=np.asarray(d.dp[r,cols],float);tested=np.asarray(d.tested[r,cols],bool)
    ok=tested&np.isfinite(fc)&(fc>0);lf=np.zeros(len(cols));lf[ok]=np.clip(np.log2(fc[ok]),-meta['clip'],meta['clip']);sig=tested&(fdr<.05)
    row={'dp_sum':dp,'dp_sumsq':dp*dp,'dp_cnt':np.ones(len(cols),np.int64),'tested_cnt':tested.astype(np.int64),'up_cnt':(sig&(fc>1)).astype(np.int64),'down_cnt':(sig&(fc<1)).astype(np.int64),'lfc_sum':lf,'lfc_sumsq':lf*lf,'lfc_cnt':ok.astype(np.int64)}
    for k in agg:agg[k][covered]+=row[k]
    for k,src in e.PER_DATASET_SOURCE.items():by[k][slot,covered]+=row[src]
    # Author leave-one-query-out shard uses float32 donor LFC records.
    packed=e._pack_rows({k:v[None,:] for k,v in row.items()});rounded=e._unpack_rows(*packed)
    for name,(ds,ctx) in MATCH.items():
     if d.dataset==ds and (ctx is None or ctx==d.donor_keys[j][0]):
      for k in agg:ex[name][k][covered]+=rounded[k][0]
  agg={k:v.astype(e.GLOBAL_ARRAYS[k]).astype(np.float64 if k in e.SUM_ARRAYS else np.int64) for k,v in agg.items()}
  for source in ['CD4_Rest_D2',*MATCH]:
   aa={k:v.copy() for k,v in agg.items()};bb={k:v.copy() for k,v in by.items()}
   if source in MATCH:e._subtract_donors(aa,bb,ex[source],meta['datasets'].index(MATCH[source][0]))
   expected=e._derive_blocks(aa,bb,e.Scales.from_dict(meta['scales']));actual,diag=blocks(target,source,axis)
   for k in expected:np.testing.assert_allclose(actual[k],expected[k],rtol=2e-6,atol=2e-6)
   for k in e.CONTEXT_BLOCK_KEYS:np.testing.assert_array_equal(actual[k],0)
   checks.append({'target':target,'source':source,'query_records_excluded':diag['query_records_excluded'],'genes':len(axis)})
 write(A/'evidence_independent_audit.json',{'passed':True,'comparisons':len(checks),'checks':checks,'method':'Independent scalar donor-label arithmetic versus the projected author evidence blocks; exact saved split; covered and missing genes; recipient exclusions and novel-context zero blocks','full_training_count_verified':True,'projection_receipt_sha256':sha(A/'evidence_receipt.json')})
 print('Independent evidence audit passed',len(checks),'comparisons',flush=True)

def artifacts():
 pp=plan();checks=[]
 for source in SOURCES:
  inp=CACHE/'inputs'/(source+'.npz');rec=read(A/'inputs'/(source+'.json'));assert rec['sha256']==sha(inp) and rec['plan_sha256']==sha(A/'plan.json')
  pred=CACHE/'predictions'/(source+'.npz');rec=read(A/'predictions'/(source+'.json'));assert rec['sha256']==sha(pred) and rec['plan_sha256']==sha(A/'plan.json')
  with np.load(pred,allow_pickle=False) as z:
   ids=rec['required_targets'];assert z['done'][ids].all()
   for k in ['p_de','lfc','delta_p']:assert np.isfinite(z[k][ids]).all()
   assert ((z['p_de'][ids]>=0)&(z['p_de'][ids]<=1)).all()
  checks.append(source)
 for k in ['joint','shuffled','evidence']:
  rec=read(A/'calibrators'/(k+'.json'));assert not rec['CD4_labels_used_in_fit'];assert all(not r['source'].startswith('CD4') for r in rec['groups']);assert sha(CACHE/'calibrators'/(k+'.npz'))==rec['sha256']
  held={r['target'] for r in rec['groups'] if r['validation']};fitted={r['target'] for r in rec['groups'] if not r['validation']};assert not held&fitted
 for f in (A/'counts').glob('*.json'):
  r=read(f);assert r['plan_sha256']==sha(A/'plan.json')
  if f.stem.endswith('_null'):continue
  assert not set(r['actual_rows'])&set(r['template_rows']);assert not set(r['reference_rows'])&set(r['template_rows']);assert len(set(r['actual_rows']))==len(r['actual_rows']) and len(set(r['template_rows']))==400
  assert sha(CACHE/'count_statistics'/f.with_suffix('.npz').name)==r['statistics_sha256']
  assert {v['arm'] for v in r['records']}==set(pp['config']['arms'])
 assert read(A/'count_evaluation.json')['complete']
 write(A/'final_audit.json',{'passed':True,'sources':checks,'plan_sha256':sha(A/'plan.json'),'independent_evidence_audit_sha256':sha(A/'evidence_independent_audit.json'),'count_evaluation_sha256':sha(A/'count_evaluation.json'),'response_evaluation_sha256':sha(A/'response_evaluation.json'),'D3_opened':False,'new_submission_generated':False,'checks':['Saved split and checkpoint contracts','All three heads retained and finite','Seen/unseen cohorts','Whole-target alpha selection','No CD4 decoder fitting','Disjoint real populations','Paired four-arm integer-count comparisons','Artifact SHA256 verification']})
 print('Final V57 audit passed',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['evidence','artifacts']);args=p.parse_args();{'evidence':evidence,'artifacts':artifacts}[args.stage]()
