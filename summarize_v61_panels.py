"""Limited common-control panel diagnostics; never a 300-target score forecast."""
import numpy as np
import pandas as pd
import v61_pooled_state_response as v
from audit_transfer_validation import pds,corrected_expression

def main():
 v.activate();m=v.m;assert m.read(m.A/'count_assessment.json')['complete'];genes,targets=m.axes();panels={}
 for path in sorted((m.A/'counts').glob('*.json')):
  if path.stem.endswith('_null'):continue
  r=m.read(path);record=r['records'][0];file=m.CACHE/'count_statistics'/path.with_suffix('.npz').name;assert m.sha(file)==r['statistics_sha256'];z=m.load(file);key=(record['source'],record['seed']);panels.setdefault(key,[]).append((record['target'],z))
 rows=[];nonpanel=~np.isin(genes,targets)
 for (source,seed),items in panels.items():
  items=sorted(items,key=lambda a:a[0]);names=[x[0] for x in items];assert len(names)==len(set(names));truth=np.array([z['truth_profile'] for _,z in items]);tj=np.array([z['truth_jackknife'] for _,z in items]);control=items[0][1]['control_profile'];cj=items[0][1]['control_jackknife']
  for _,z in items:np.testing.assert_array_equal(z['control_profile'],control);np.testing.assert_array_equal(z['control_jackknife'],cj)
  own=np.array([np.flatnonzero(genes==t)[0] for t in names]);perfect=pds(truth[:,nonpanel],truth[:,nonpanel],control[nonpanel]);assert abs(perfect-1)<1e-12
  for arm in ['v24','composition50','joint50','shuffled_joint50','pooled50','state_cells50']:
   pred=np.array([z[arm+'_profile'] for _,z in items]);pj=np.array([z[arm+'_jackknife'] for _,z in items]);raw=pds(pred[:,nonpanel],truth[:,nonpanel],control[nonpanel]);constant=pds(np.tile(pred.mean(0),(len(items),1))[:,nonpanel],truth[:,nonpanel],control[nonpanel]);assert abs(constant-.5)<1e-12
   cyclic=np.mean([pds(np.roll(pred,shift,axis=0)[:,nonpanel],truth[:,nonpanel],control[nonpanel]) for shift in range(len(items))]);assert abs(cyclic-.5)<1e-12
   expr=corrected_expression(pred,truth,control,pj,tj,cj,own)
   rows.append(dict(source=source,seed=seed,arm=arm,unique_targets=len(items),raw_panel_PDS=raw,constant_PDS=constant,balanced_label_shift_PDS=cyclic,**expr))
 frame=pd.DataFrame(rows);frame.to_csv(m.A/'panel_diagnostics.csv',index=False);m.write(m.A/'panel_diagnostics_receipt.json',dict(complete=True,panels=len(panels),comparisons=len(frame),script_sha256=m.sha(__file__),metric_helper_sha256=m.sha(m.ROOT/'audit_transfer_validation.py'),data_sha256=m.sha(m.A/'panel_diagnostics.csv'),interpretation='Only 4 or 8 unique target identities per source/seed, not 300-target official discrimination. Common NTC profile verified bitwise. Every panel excludes all 300 target genes for PDS. Expression numerator/denominator are retained separately; unresolved denominators are missing, not replaced. No official anchors, private score forecast or model promotion.',D3_opened=False))
 print(frame.groupby(['source','arm'])[['raw_panel_PDS','common_expression_numerator','common_expression_denominator']].mean().to_string(),flush=True)
if __name__=='__main__':main()
