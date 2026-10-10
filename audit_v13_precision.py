"""Quantify incidental floating-point differences on targets with no KD correction."""
from pathlib import Path
import json
import numpy as np
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v13_response_programs/analysis'

def main():
 source=ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'
 with np.load(source,allow_pickle=False) as z:
  genes=z['genes'];targets=z['targets'];own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);rt=z['row_targets']
  control=z['profiles'][z['row_contexts'],own[rt]]
 affected=np.unique(rt[(~np.isfinite(control))|(control<=5)]);unaffected=np.setdiff1d(np.arange(len(targets)),affected)
 olddir=ROOT/'experiments/v11_consensus_coverage/analysis';oldmodel=json.loads((olddir/'reference_model.json').read_text(encoding='utf-8'))
 assert build.sha256(olddir/'reference_model_arrays.npz')==oldmodel['model_arrays_sha256']
 newmodel=json.loads((A/'reference_model.json').read_text(encoding='utf-8'))
 assert build.sha256(A/'reference_model_arrays.npz')==newmodel['model_arrays_sha256']
 with np.load(olddir/'reference_model_arrays.npz',allow_pickle=False) as z:old=z['effects']
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:new=z['effects']
 difference=abs(new[:,unaffected]-old[:,unaffected]);changed=np.any(new!=old,axis=2)
 maximum=float(difference.max()) if difference.size else 0.
 assert maximum<=1e-6,'Unexpected substantive change outside the KD-corrected targets'
 build.write_json(A/'numerical_precision_audit.json',{'passed':True,'KD_corrected_targets':len(affected),'other_targets':len(unaffected),
  'maximum_absolute_LFC_change_on_other_targets':maximum,'mean_absolute_LFC_change_on_other_targets':float(difference.mean()),
  'all_context_target_response_vectors_bitwise_changed':int(changed.sum()),
  'explanation':'The V13 pooling implementation retains float64 intermediates where V11 casts pooled mean/variance to float32 before shrinkage. These low-bit changes are included in the frozen, count-tested candidate. They prevent exact reuse even where KD is unchanged; do not describe V13 as a bitwise single-change ablation.',
  'scope':'Post-freeze arithmetic audit only; no predictor changes or retuning.',
  'old_arrays_sha256':oldmodel['model_arrays_sha256'],'new_arrays_sha256':newmodel['model_arrays_sha256'],'script_sha256':build.sha256(Path(__file__))})
 print('Precision audit:',len(affected),'KD-affected targets;',len(unaffected),'others; max other delta',maximum,'changed groups',int(changed.sum()))
if __name__=='__main__':main()
