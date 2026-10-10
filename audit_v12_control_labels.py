"""Authoritative control-guide identity audit before any Flex model fitting."""
import json,hashlib
from pathlib import Path
import requests,pandas as pd
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v12_distribution_calibrated/analysis';D=ROOT/'reference/v12_resources'
def main():
 url='https://satijalab.org/viperturb-seq/computational-tools/guest_list_protospacer_reference.csv';p=D/'guest_list_protospacer_reference.csv'
 if not p.exists():
  r=requests.get(url,timeout=90);r.raise_for_status();p.write_bytes(r.content)
 f=pd.read_csv(p,header=None,names=['guide_id','protospacer','gene']);assert len(f)==57050
 nt=f[f.gene.eq('NO-TARGET')];ctrl=f[f.gene.eq('CTRL')]
 assert len(nt)==496 and len(ctrl)==3 and set(nt.guide_id)=={f'NO-TARGET_{i}' for i in range(496)}
 record={'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'library_guides':len(f),'true_non_targeting_gene_label':'NO-TARGET','non_targeting_guides':496,'CTRL_is_target_gene_with_three_guides':True,'control_guide_names_in_Seurat':'NO-TARGETg0 through NO-TARGETg495','control_labels_verified_before_first_Flex_response_fit':True}
 (A/'guide_reference.json').write_text(json.dumps(record,indent=2)+'\n')
 parent=A/'flex_source_preprocessing_plan.json'
 correction={'parent_preprocessing_plan_sha256':hashlib.sha256(parent.read_bytes()).hexdigest(),'reason':'Initial whitelist wrongly treated the human gene CTRL as a negative-control label and omitted NO-TARGET. Corrected using the author library, before fitting any Flex effects.',
 'true_control_cells_observed_binB':7949,'supersedes_initial_69_control_count':True,'true_control_policy':'Only NO-TARGET cells whose guide ID is in the verified 496-guide reference. Deduplicate exact repeated cells across bins. Match two lanes as prospectively planned.', 'no_model_fitted_with_wrong_labels':not (D/'flex_source/summaries.npz').exists(),'guide_reference_sha256':hashlib.sha256((A/'guide_reference.json').read_bytes()).hexdigest()}
 assert correction['no_model_fitted_with_wrong_labels']
 (A/'flex_control_identity_correction.json').write_text(json.dumps(correction,indent=2)+'\n')
 print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':main()
