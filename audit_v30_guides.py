"""Audit actual guide-resolution limits; no inferred guides or treatment fitting."""
from pathlib import Path
import json
import pandas as pd
import v30_methods as m
import build_submission as b
FILES={'H1':('reference/v5_h1/selected_cells.csv','target_gene','guide_id'),'HCT116':('reference/v7_sources/HCT116/selected_cells.csv','gene_target','guide_target'),'HEK293T':('reference/v7_sources/HEK293T/selected_cells.csv','gene_target','guide_target'),'CD4_Rest_D2':('experiments/v22_state_response/analysis/cd4_selected/cells.csv','gene_target',None),'CD4_Stim48_D2':('experiments/v24_cd4_response/analysis/holdout/cells.csv','gene_target',None)}
records=[]
for name,(rel,target,guide) in FILES.items():
 p=m.ROOT/rel;f=pd.read_csv(p);columns=list(f);guide=guide or next((c for c in ['guide_id','guide_target','guide','sgRNA','guide_ids','guide_identity'] if c in f),None);target=target if target in f else next(c for c in ['target_gene','gene_target','target'] if c in f);f=f[~f[target].astype(str).str.lower().isin(['non-targeting','non_targeting','control'])];rec={'family':name,'metadata_sha256':b.sha256(p),'metadata_columns':columns,'target_column':target,'guide_column':guide,'targets':int(f[target].nunique())}
 if guide:
  n=f.groupby(target)[guide].nunique();rec.update(targets_multiple_guide_labels=int((n>1).sum()),max_labels_per_target=int(n.max()),composite_label_fraction=float(f[guide].astype(str).str.contains('|',regex=False).mean()),interpretation=('Composite construct label cannot identify separate effects of co-delivered guides.' if f[guide].astype(str).str.contains('|',regex=False).any() else 'Separately resolved guide labels available; shared donors/lanes are not independent culture replicates.'))
 else:rec['interpretation']='No separately resolved guide column in the extracted metadata; no guide effects fitted.'
 records.append(rec)
m.write(m.A/'guide_resolution_audit.json',{'records':records,'construct_model_scope':'The V30 hierarchical variant uses original Replogle pseudobulk promoter/transcript constructs with measured fold_expr, not invented individual guide identities.'});print(json.dumps(records,indent=2))
