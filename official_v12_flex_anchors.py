"""Compute real official generic/replicate comparators on public development panels."""
import json,argparse,hashlib,gc
from pathlib import Path
from dataclasses import replace
import polars as pl
from cell_eval2 import EvalConfig,build_generic_baseline
from cell_eval2.anchor import compute_replicate_anchor
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v12_distribution_calibrated/analysis/flex_stage'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 manifest=json.loads((A/'official_panel_manifest.json').read_text());out=A/'official_anchors';out.mkdir(exist_ok=True)
 cfg=EvalConfig.from_preset('vcc2026');cfg.device='cpu';cfg.num_threads=2;cfg.de=replace(cfg.de,backend='scanpy')
 records=[]
 for r in manifest['panels']:
  file=ROOT/r['file'];assert sha(file)==r['sha256'];name=r['context'];side=out/(name+'.json')
  if side.exists():
   saved=json.loads(side.read_text());assert saved['real_sha256']==r['sha256'];records.append(saved);continue
  base=build_generic_baseline(str(file),config=cfg,seed=2026100941,allow_degenerate=True)
  base.agg.write_csv(out/(name+'_baseline.csv'));base.results.write_csv(out/(name+'_baseline_targets.csv'))
  splits,anchor=compute_replicate_anchor(str(file),config=cfg,base_seed=2026100947,n_splits=5)
  splits.write_csv(out/(name+'_replicate_splits.csv'));anchor.write_csv(out/(name+'_replicate.csv'))
  record={'context':name,'family':r['family'],'real_sha256':r['sha256'],'n_splits':5,'official_commit':'0ce25cf4c495e76e60dfca4ddc8f822507aa0e18',
   'baseline_emit':'dispersed','backend':'scanpy CPU','files_sha256':{p.name:sha(p) for p in out.glob(name+'_*.csv')},'scope':'Public source context comparator; not private competition anchor.'}
  side.write_text(json.dumps(record,indent=2)+'\n');records.append(record);del base;gc.collect();print('Official anchors completed',name,flush=True)
 (out/'manifest.json').write_text(json.dumps({'completed':True,'panels':records,'script_sha256':sha(Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
