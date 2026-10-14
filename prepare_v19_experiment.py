"""Freeze covariance-response tests after the independently failed PIE pilot."""
import json
import build_submission as b
from v19_methods import ROOT,A,NAMES

def main():
 A.mkdir(parents=True,exist_ok=True)
 if (A/'plan.json').exists():raise FileExistsError('Preserve V19 freeze')
 parent=ROOT/'experiments/v18_pie_expert/analysis';old=json.loads((parent/'plan.json').read_text());failed=json.loads((parent/'selection.json').read_text());assert not failed['new_candidate_passes'];assert not (parent/'confirmation_result.json').exists()
 deps=['v19_methods.py','prepare_v19_experiment.py','run_v19_experiment.py','assess_v19_experiment.py','audit_v19_experiment.py','test_v19_covariance.py','v16_methods.py','fit_v7_response.py','v12_cpu_metrics.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py']
 record={'version':19,'contexts':old['contexts'],'seeds':old['seeds'],'candidates':NAMES,'eligible':['covariance'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(parent/'plan.json'),'previous_failed_selection_sha256':b.sha256(parent/'selection.json'),
 'model':'Training-only common-measured top2048 within-batch/depth-adjusted variance genes, all300 panel genes excluded. Shared16-program standardized control SVD (seed19019,n_iter3); source/query covariances .8*C+.2*diag(C)+.05*mean(diag(C))*I. Source shifts: batch- and NTC-depth-regression-adjusted log1p CP10k treatment program means; >=40 original cells, family weight n/(n+100). Original full-library denominator, no treatment-derived query inputs.',
 'correction':'For each source, bounded query_covariance @ inverse(source_covariance) minus identity, applied to observed source program shift. Singular values bounded [.5,2]. Pool supported source families; decode via shared basis and local query mean(C/(1+C))*ln2 derivative floor.1; clip +/-1 log2, weight.25, existing expression gate, total final LFC +/-3. Unsupported targets/genes or query coverage<.8 retain V16. V16 residual10 emitter unchanged.',
 'controls':'Diagonal source/query covariance; frozen permutation of query program labels; V16 with non-own full-axis L1 amplitude matched to covariance response. These share basis, data and correction budget. No per-context whitening to identity and no query-treatment amplitude fitting.',
 'release':'Same available raw and worst-family guards, >=3 metrics improve including DE. Must beat diagonal, shuffled and amplitude controls on >=2 fidelity/reach/Jaccard measures. Same once-only rule on reserved verification, no post-failure retuning. Missing NMAE remains missing.',
 'verification':'Eight corrected V17 pairs, not opened by V17 or failed V18. Previously accessible public biology, not blind biological validation. Same14 development pairs reused explicitly. Covariance transfer is a newly frozen algorithm, not a retuned PIE threshold.',
 'scope':'Adaptation inspired by CIPHER equilibrium response; not reproduction of CIPHER rank1 oracle amplitude model and not causal identification. Controls-only query covariance can change expected response directions; joint count generator unchanged.'}
 b.write_json(A/'plan.json',record);print('V19 frozen',b.sha256(A/'plan.json'),flush=True)
if __name__=='__main__':main()
