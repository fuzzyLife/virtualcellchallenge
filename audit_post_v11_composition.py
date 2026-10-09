"""Quantify own-target knockdown's indirect CPM effect using observed controls only."""
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'analysis/post_v11_research'

def main():
 targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene;records=[];inputs={}
 for context in 'ABC':
  path=ROOT/f'analysis/gene_statistics_{context}.csv';d=pd.read_csv(path).set_index('gene_name');inputs[str(path.relative_to(ROOT))]=build.sha256(path)
  np.testing.assert_allclose(d.mean_cpm.sum(),1e6,rtol=1e-10)
  for target in targets:
   share=float(d.loc[target,'mean_cpm']/1e6);record={'context':context,'target_gene':target,'own_mean_CPM':share*1e6}
   for knockdown in [.5,.8,.9]:record[f'other_gene_log2FC_at_kd{int(knockdown*100)}']=float(-np.log2(1-knockdown*share))
   record['other_gene_log2FC_kd90_minus_kd80']=record['other_gene_log2FC_at_kd90']-record['other_gene_log2FC_at_kd80'];records.append(record)
 frame=pd.DataFrame(records);frame.to_csv(OUT/'on_target_composition.csv',index=False)
 build.write_json(OUT/'composition_manifest.json',{'script_sha256':build.sha256(Path(__file__)),'inputs':inputs,
  'scope':'Exact mean-composition algebra for own-target-only reduction and renormalization. Not a raw-cell DE test, capture model or biological total-RNA estimate. No parameter selection.'})
 print(frame.nlargest(8,'other_gene_log2FC_at_kd90').to_string(index=False))
if __name__=='__main__':main()
