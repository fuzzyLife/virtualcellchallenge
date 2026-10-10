"""Audit RPE1 coverage using gene symbols parsed from composite intervention IDs."""
import re,json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v13_response_programs/analysis'

def main():
 path=ROOT/'reference/v13_resources/rpe1_raw_bulk_01.h5ad'
 data=ad.read_h5ad(path,backed='r');_,genes,targets=build.read_spec(ROOT)
 rows=[]
 for i,label in enumerate(data.obs_names.astype(str)):
  match=re.fullmatch(r'\d+_(.+?)_[^_]+_(ENSG\d+|nan)',label)
  symbol=match.group(1) if match else None
  control='non-targeting' in label
  assert symbol is not None or control,label
  n=float(data.obs.num_cells_filtered.iloc[i]);fold=float(data.obs.fold_expr.iloc[i])
  rows.append({'intervention_id':label,'target_symbol':symbol,'target_ensembl':match.group(2) if match and match.group(2)!='nan' else None,
    'non_targeting':control,'core_control':bool(data.obs.core_control.iloc[i]),'num_cells_filtered':n,'fold_expr':fold,
    'competition_target':symbol in targets if symbol else False,'quality_pass':bool(n>=30 and np.isfinite(fold) and fold<=.5)})
 frame=pd.DataFrame(rows);frame.to_csv(A/'rpe1_intervention_coverage.csv',index=False)
 present=sorted(frame.loc[frame.competition_target,'target_symbol'].unique().tolist())
 qc=sorted(frame.loc[frame.competition_target&frame.quality_pass,'target_symbol'].unique().tolist())
 measured=data.var.gene_name.astype(str).to_numpy();overlap=sorted(set(genes)&set(measured))
 old=A/'rpe1_coverage.json';archive=A/'rpe1_coverage_initial.json'
 if old.exists() and not archive.exists():archive.write_bytes(old.read_bytes())
 build.write_json(old,{'targets_present':present,'targets_QC':qc,'core_control_profiles':int(frame.core_control.sum()),
  'total_profiles':len(frame),'parsed_gene_profiles':int(frame.target_symbol.notna().sum()),
  'competition_output_gene_overlap':len(overlap),'measured_genes':len(measured),
  'identifier_policy':'Parse numeric_construct, HGNC symbol, promoter/transcript token, and ENSG identifier (or explicit published nan); assert every unparsed row is a non-targeting control. Do not compare composite IDs directly with HGNC symbols.',
  'scope':'Bulk means, not single-cell counts; quality n>=30 and published fold_expr<=.5. This audit does not change the frozen V13 predictor.',
  'source_sha256':build.sha256(path),'script_sha256':build.sha256(Path(__file__)),
  'coverage_table_sha256':build.sha256(A/'rpe1_intervention_coverage.csv')})
 data.file.close();print('Parsed RPE1 symbols:',len(present),'overlapping interventions;',len(qc),'QC;',len(overlap),'overlapping output genes')
if __name__=='__main__':main()
