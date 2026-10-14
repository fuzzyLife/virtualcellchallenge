"""Freeze the same metadata-selected CD4 D2 reserve for the V23 candidate."""
import json
from pathlib import Path
from v23_methods import ROOT,A,write,load_plan
import v22_methods as parent
import build_submission as b

def main():
 load_plan();dest=A/'reserve_plan.json'
 if dest.exists():raise FileExistsError('Preserve reserve identities')
 source=parent.A/'reserve_plan.json';plan=json.loads(source.read_text());amend=json.loads((parent.A/'reserve_source_amendment.json').read_text());assert not amend['treatment_metrics_opened'];assert len(plan['cells'])==8 and 'D2_Rest' in json.loads((parent.A/'cd4_selected/extraction.json').read_text())['source'];assert not (parent.A/'confirmation_evaluation.json').exists()
 plan['parent_metadata_reserve_plan_sha256']=b.sha256(source);plan['parent_plan_sha256']=b.sha256(A/'plan.json');plan['source_amendment_sha256']=b.sha256(parent.A/'reserve_source_amendment.json');plan['source_extraction_sha256']=b.sha256(parent.A/'cd4_selected/extraction.json');plan['selection_basis']='Same eight alphabetical metadata-only powered donor2 targets; no held treatment scores seen. Both stages exclude all CD4 from training.'
 write(dest,plan);print('V23 raw-cell reserve frozen:',[v['target'] for v in plan['cells']],flush=True)
if __name__=='__main__':main()
