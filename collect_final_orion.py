"""Reuse reviewed Orion extractor in isolation for only newly required final targets."""
from pathlib import Path
import argparse,types,sys,shutil,json
import pandas as pd
import prepare_v7_sources as old
import audit_v7_source_offsets as offsets
import build_submission as b
R=Path(__file__).resolve().parent;O=R/'reference/final_v91/orion';A=R/'experiments/final_v91/analysis'
def main(line):
 O.mkdir(parents=True,exist_ok=True);new=sorted(set(pd.read_csv(R/'final/pert_counts.csv').target_gene)-set(pd.read_csv(R/'pert_counts.csv').target_gene));pd.DataFrame({'target_gene':new}).to_csv(O/'new_targets.csv',index=False)
 source=Path(old.__file__).read_text(encoding='utf-8');text=source.replace("targets=pd.read_csv(ROOT/'pert_counts.csv')","targets=pd.read_csv(OUT/'new_targets.csv')")
 assert text!=source
 module=types.ModuleType('final_orion_isolated');module.__file__=old.__file__;exec(compile(text,old.__file__,'exec'),module.__dict__);module.OUT=O
 offsets.OUT=O
 for c in [line]:
  dest=O/c;dest.mkdir(exist_ok=True);orig=R/'reference/v7_sources'/c/'expression_fragment_audit.json';rec=json.loads(orig.read_text());assert rec['revision']==old.REVISION
  if not (dest/orig.name).exists():shutil.copyfile(orig,dest/orig.name)
  module.extract(c,cap=400,workers=3)
 b.write_json(A/(line+'_new_source.json'),dict(extractor_original_sha256=b.sha256(Path(old.__file__)),offset_auditor_sha256=b.sha256(Path(offsets.__file__)),adapter_sha256=b.sha256(Path(__file__)),changes='Output folder and target-list path only; original pinned release, expression offset audit, library checks and sampler unchanged.',new_targets=new,extraction=json.loads((O/line/'extraction.json').read_text())))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('line',choices=['HCT116','HEK293T']);args=p.parse_args();main(args.line)
