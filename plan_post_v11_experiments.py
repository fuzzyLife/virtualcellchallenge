"""Freeze future post-V11 confirmation targets using metadata, before new count tests."""
import hashlib,json
from pathlib import Path
import pandas as pd
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'analysis/post_v11_research'

def main():
 path=OUT/'future_experiment_plan.json'
 if path.exists():raise FileExistsError('Future plan is already frozen; do not silently retune it.')
 used=set();inputs={}
 for directory in ['v8_corrected_transfer','v9_activation_preserving','v10_population_transport','v11_consensus_coverage']:
  source=ROOT/f'experiments/{directory}/analysis/validation_plan.json';p=json.loads(source.read_text());inputs[str(source.relative_to(ROOT))]=build.sha256(source)
  used|={c['target'] for context in p['contexts'] for c in context['cells']}
 for family in ['H1','HCT116','HEK293T']:
  source=ROOT/f'experiments/v7_multisource_adaptive/analysis/cell_split_{family}.json';p=json.loads(source.read_text());used|=set(p['targets']);inputs[str(source.relative_to(ROOT))]=build.sha256(source)
 latest=json.loads((ROOT/'experiments/v11_consensus_coverage/analysis/validation_plan.json').read_text());targets=set(pd.read_csv(ROOT/'pert_counts.csv').target_gene)
 chosen=set();confirm=[]
 for family in ['HCT116','HEK293T']:
  source=ROOT/f'reference/v7_sources/{family}/selected_cells.csv';meta=pd.read_csv(source);inputs[str(source.relative_to(ROOT))]=build.sha256(source)
  for context in [c for c in latest['contexts'] if c['family']==family]:
   counts=meta[meta['sample'].isin(context['batches'])].gene_target.value_counts()
   available=[t for t in targets-used-chosen if counts.get(t,0)>=150]
   available.sort(key=lambda t:hashlib.sha256(f'post-v11-confirm:{context["context"]}:{t}'.encode()).hexdigest())
   selected=available[:3];assert len(selected)==3;chosen|=set(selected)
   confirm.append({'family':family,'context':context['context'],'targets':selected,'available_treated_cells':{t:int(counts[t]) for t in selected},
    'technical_strata_from_v11':True,'control_population_indices_from':str(Path('experiments/v11_consensus_coverage/analysis/validation_plan.json'))})
 assert not chosen&used
 tests=[
  {'id':'D1','hypothesis':'Matched Flex genetic data improve source transfer after controlling for source family.','design':'Extract the 300 targets plus deduplicated NTCs from all VIPerturb bins. Verify raw counts, output gene mask, guide knockdown, cell IDs, batch and library axes. Compare original source pool, Flex replacement within K562, and within-family combined summaries at identical V11 emission. Hold out all K562 datasets together; do not count platforms as independent lineages.','negative_control':'Shuffle target labels within the added source; compare a size-matched downsample of existing K562.'},
  {'id':'D2','hypothesis':'Gene-specific distribution calibration reduces expression error without erasing direction/coverage.','design':'Keep V11 effects fixed; compare exact V11 emission, moment-matched gamma-Poisson, and one 25% gamma/75% compact template. Test null controls first, then actual public perturbations at genuine cell counts. Analyze direction, reach and six-metric panel jointly.','negative_control':'Real control resampling with no response; source target-label permutation; own-target knockdown disabled for nulls.'},
  {'id':'D3','hypothesis':'The strongest response confidence proxies are miscalibrated across contexts.','design':'Fit a bounded transfer-heterogeneity calibration and false-sign reliability curve using inner held-family development folds. Compare scalar empirical Bayes with a multivariate/covariance or pathway prior. Covariances must use training effects only; a novel family has no measured effect and requires an explicit prediction rule.','negative_control':'Permute source identity/uncertainty while preserving marginal weights. Report held-family interval coverage and sign error by abundance/support, not only aggregate MSE.'},
  {'id':'D4','hypothesis':'Cell-state and responder heterogeneity bias source pseudobulk effects.','design':'Compare raw source CPM effects with control-factor-adjusted means and a rank-one perturbation-factor pilot inspired by RNA fingerprinting. Include soft responder mixtures; do not delete weak responders from validation truth. Calibrate factors back to mean/pooled count effects before emission.','negative_control':'Random factors of identical rank; NTC pseudo-perturbations; target-shuffled response factors. Query challenge controls alone cannot identify an unobserved treated fingerprint.'},
  {'id':'D5','hypothesis':'Pathway structure improves residual target-specific response rather than shared stress recovery.','design':'Use pinned Reactome membership for a small residual prior after subtracting training-only shared response. Compare zero effect, training perturbation centroid, V11, and pathway residual. Keep official common-control metrics and report descriptive perturbation-centered diagnostics separately.','negative_control':'Size/degree-matched shuffled pathway sets; exclude each target and its held functional group when testing extrapolation.'}]
 plan={'status':'Planned; no new perturbation candidate fitted or selected.','script_sha256':build.sha256(Path(__file__)),'input_sha256':inputs,
  'baseline_v11_model_sha256':build.sha256(ROOT/'experiments/v11_consensus_coverage/analysis/reference_model.json'),
  'seeds':[2026100911,2026100927],'prediction_cells':400,'reference_modes':['sampled400','full_independent'],
  'new_count_confirmation':confirm,'excluded_previous_count_targets':sorted(used),
  'split_scope':'Selected from target labels and cell-count metadata only. Previously inspected response summaries mean these are new count-test targets, not pristine biological blind data. H1 remains reused development; no fresh 400-cell confirmation is claimed.',
  'selection_rules':['Freeze a candidate using development only before opening new count confirmation. Treat a failed confirmation as failure and do not retune against it.',
   'Exclude the whole evaluated source family, including all platforms and correlated releases. Keep common target-independent templates and independent references.',
   'Require complete six-metric panels; equal-weight source families and technical strata; paired seeds. Report worst-family and actual-cell-count limitations.',
   'Protect V11 fidelity and reach: no aggregate fall greater than 0.01; no family fall greater than 0.03. NMAE must not worsen by more than 0.01. PDS must not fall more than 0.02.',
   'Expression error must improve on development and confirmation. Treat these provisional tolerances as fixed screening rules, not an official anchored score estimate.',
   'Report neutral-call fraction and detection/variance errors in every context; do not claim a null-safety pass from a mean across contexts. Do not inflate observed sample size by duplicating cells.',
   'Limit the initial count-stage to at most four candidates including V11; source, response and emitter changes must be tested separately before combination.'],
  'tests':tests}
 build.write_json(path,plan)
 print('Frozen new count-confirmation targets:',[(c['context'],c['targets']) for c in confirm])
if __name__=='__main__':main()
