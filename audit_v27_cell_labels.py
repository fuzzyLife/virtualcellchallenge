"""Additional held-out actual-cell labels audit, preserving the V27 input auditor."""
from pathlib import Path
import numpy as np
from v27_methods import ROOT,A,load_plan,write
from v26_methods import raw,family_key
from fit_v7_response import read_cache
import build_submission as b

def main():
 p=load_plan();genes=read_cache()['genes'];checked=0
 for ctx in p['contexts']:
  stages=[s for s in ['null','treatment'] if (A/f'{s}_{family_key(ctx["family"])}_evaluation.json').exists()]
  if not stages:continue
  x,m,measured,_=raw(ctx['family'],genes)
  for stage in stages:
   for cell in (ctx['null_cases'] if stage=='null' else ctx['cells']):
    ids=cell['actual'];assert len(set(ids))==len(ids);assert (m.gene_target.to_numpy(str)[ids]==('Non-Targeting' if stage=='null' else cell['target'])).all();checked+=1
 write(A/'cell_label_audit.json',{'passed':True,'actual_cases_checked':checked,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__))})
 print('V27 independent actual-cell label audit passed',checked,flush=True)
if __name__=='__main__':main()
