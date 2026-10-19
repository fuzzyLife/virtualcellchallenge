"""V67: bounded intervention-network rewiring from measured CRISPR responses."""
from pathlib import Path
import argparse
import gc
from functools import lru_cache
import numpy as np
from threadpoolctl import threadpool_limits
from v57_pie import read, write, sha
import v62_conserved_response as old

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'experiments/v67_intervention_network'
A = OUT / 'analysis'
ARMS = ('network', 'network_shuffled', 'gain_only')


def graph(response, variance, valid):
    """Response rows are interventions, columns are measured target genes."""
    n = len(response)
    z = np.where(valid, np.sign(response) * np.maximum(abs(response) - 2*np.sqrt(variance), 0), 0)
    np.fill_diagonal(z, 0)
    susceptibility = np.eye(n) + np.clip(z, -3, 3).T / np.log2(.1)
    u, s, vt = np.linalg.svd(susceptibility, full_matrices=False)
    inverse = (vt.T / np.maximum(s, .25)) @ u.T
    diagonal = np.diag(inverse)
    normalized = np.divide(inverse, diagonal[:, None], out=np.eye(n), where=abs(diagonal[:, None])>1e-6)
    w = np.eye(n) - normalized
    np.fill_diagonal(w, 0)
    norm = float(np.linalg.norm(w, ord=2))
    w *= min(1., .45/max(norm, 1e-30))
    return w, dict(raw_operator_norm=norm, bounded_operator_norm=float(np.linalg.norm(w, 2)),
                   response_matrix_min_singular=float(s.min()), nonzero_edges=int(np.count_nonzero(w)))


def susceptibility(w, gain):
    h = np.linalg.solve(np.eye(len(w)) - w*np.asarray(gain)[None, :], np.eye(len(w)))
    assert np.isfinite(h).all() and (np.diag(h) > 0).all()
    return h / np.diag(h)[None, :]


def response_map(w, gain):
    if np.array_equal(gain, np.ones(len(w))) or not np.any(w):
        return np.eye(len(w))
    return np.linalg.solve(susceptibility(w, np.ones(len(w))), susceptibility(w, gain))


def availability(query, source):
    q, s = np.asarray(query, float), np.asarray(source, float)
    valid = np.isfinite(q) & np.isfinite(s) & (q >= 0) & (s > 5)
    gain = np.ones(len(s))
    # A fixed saturation convention, not an experimentally identified mechanism.
    gain[valid] = np.sqrt((q[valid]/(q[valid]+20)) / (s[valid]/(s[valid]+20)))
    return np.clip(gain, .5, 1.5)


@lru_cache(3)
def bank(excluded):
    b = old.bank(excluded)
    own = b['own']; models = []
    for fi, family in enumerate(old.FAMILIES):
        if old.GROUP[str(family)] in excluded:
            continue
        present = b['ok'][fi].any(1)
        axis = b['ok'][fi, present].all(0) if present.any() else np.zeros(len(old.data()['genes']), bool)
        # Every output used by the network must be measured for every active source intervention.
        axis[own] = False
        w, trace = graph(b['y'][fi][:, own], b['var'][fi][:, own], b['ok'][fi][:, own])
        models.append(dict(family=str(family), family_index=fi, w=w, axis=axis,
                           own=own, profile=b['profiles'][fi], present=present,
                           responses=b['centered'][fi][:, axis].astype(float),
                           variance=b['var'][fi][:, axis].astype(float), trace=trace))
    return models


def correction(query, models, n_genes, shuffle=False):
    own = models[0]['own']; numerator = np.zeros((len(own), n_genes)); denominator = np.zeros_like(numerator)
    traces = []
    for source in models:
        gain = availability(query[own], source['profile'][own])
        if shuffle:
            gain = gain[np.random.default_rng(67067).permutation(len(own))]
        mapping = response_map(source['w'], gain)
        delta = (mapping.T - np.eye(len(own))) @ source['responses']
        weight = old.WEIGHTS[source['family_index']] * np.clip(.04 / np.maximum(source['variance'], .0025), .1, 4.)
        weight *= source['present'][:, None]
        numerator[:, source['axis']] += weight * delta
        denominator[:, source['axis']] += weight
        traces.append(dict(family=source['family'], genes=int(source['axis'].sum()),
                           active_interventions=int(source['present'].sum()),
                           gain_min=float(gain.min()), gain_max=float(gain.max()),
                           mapping_distance=float(np.linalg.norm(mapping-np.eye(len(own)))), **source['trace']))
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0), traces, denominator > 0


def predict(query, baseline, excluded):
    models = bank(tuple(sorted(excluded)))
    delta, trace, supported = correction(query, models, baseline.shape[1])
    shuffled, shtrace, _ = correction(query, models, baseline.shape[1], True)
    # Matching projection tests whether a graph gain is only an amplitude correction.
    gain = np.divide(np.sum(delta*baseline, 1), np.sum(baseline**2, 1),
                     out=np.zeros(len(baseline)), where=np.sum(baseline**2, 1)>0)
    scalar = gain[:, None]*baseline*supported
    corrections = dict(network=delta, network_shuffled=shuffled, gain_only=scalar)
    own = models[0]['own']
    result = {}
    for arm, v in corrections.items():
        v = .5*np.clip(v, -1, 1)
        v[:, ~np.isfinite(query)] = 0
        # All panel output genes retain the comparator: they serve as network mediators only.
        v[:, own] = 0
        result[arm] = np.clip(baseline+v, -3, 3).astype(np.float32)
        result[arm][np.arange(len(own)), own] = 0
    return result, dict(sources=trace, shuffled_sources=shtrace,
                        maximum_raw_correction=float(abs(delta).max()),
                        mean_absolute_accepted_correction=float(np.mean(abs(result['network']-baseline))))


def freeze():
    A.mkdir(parents=True, exist_ok=True)
    assert not (A/'plan.json').exists()
    deps = ['v67_intervention_network.py','run_v67_experiment.py','test_v67_network.py',
            'audit_v67_boundaries.py','v62_conserved_response.py','v15_ensemble_methods.py','v7_response_methods.py',
            'run_v35_v36.py','run_v37_ensemble.py','v25_methods.py','v12_cpu_metrics.py','v16_methods.py']
    inputs = ['experiments/v7_multisource_adaptive/analysis/source_provenance.json',
              'experiments/v7_multisource_adaptive/analysis/source_cache.npz',
              'experiments/v66_assay_strength/analysis/leaderboard_result.json',
              'experiments/v66_assay_strength/analysis/release_manifest.json',
              'experiments/v34_residual_flow/analysis/plan.json']
    from run_v35_v36 import contexts
    for ctx in contexts():
        name=ctx['family']
        inputs += ['experiments/v66_assay_strength/analysis/predictions/'+name+'.json',
                   'experiments/v37_component_ensemble/analysis/cache/'+name+'/receipt.json']
    write(A/'plan.json', dict(version=67, dependencies={n:sha(ROOT/n) for n in deps},
          inputs={n:sha(ROOT/n) for n in inputs}, arms=ARMS, primary='network',
          hypothesis='Use intervention-to-intervention response structure to model feedback changes with baseline target availability; transfer the resulting source response difference onto V63. No fitted neural residual, STRING graph or observational covariance inverse.',
          model='Threshold cross-target log2 responses by 2 source SE, form normalized susceptibility S=I+thresholded_response.T/log2(.1), invert with singular floor .25, row-normalize inverse(S) by its diagonal, form W=I-normalized_inverse and set diagonal to zero and contract its spectral norm to <=.45. Query outgoing mediator activity is sqrt(saturation(query)/saturation(source)), saturation(p)=p/(p+20), clipped .5..1.5, unavailable source/query pairs unchanged. Clamp-normalize source and query resolvents (I-W*diag(activity))^-1. Mixing map is solve(source_resolvent,query_resolvent). Transport source-common-centered experimental responses, subtract original; precision-pool only outputs present for every active source intervention; apply half of bounded +/-1 correction to V63. All panel output responses remain V63. Assumptions are predictive approximations; 90% knockdown is not an infinitesimal intervention and causal network identification is not claimed.',
          controls='Unchanged V63 and V66, mediator-availability permutation (fixed seed67067), and row-wise projection of the network correction onto baseline responses (gain_only). Primary fixed before outcomes; no promotion of controls or parameter grid.',
          broad_test='All 655 existing target/context pairs with >=30 measured treatment cells, five contexts, three studies; same study exclusion as V63 including HCT/HEK together. Target SHA256(v67-target:+symbol) modulo2 partitions development/confirmation, shared across contexts. Both subsets are reused public biology, not independently blind. Predict all 300 responses before reading any treatment matrix.',
          mean_gate=dict(profile_ratio_vs_v63=.99, profile_ratio_vs_v66=.99, worst_study_profile_ratio=1.02,
                         pds_change_min=-.005, beat_both_controls_profile=True,
                         development_and_confirmation_profile_ratio_max=1.0),
          count_gate='Only a primary passing broad means proceeds to all canonical30 count pairs with seeds0/620000, unchanged residual10 emitter, V63/V66/permuted controls. Require lower expression error than both, no PDS loss>.005, FID loss>.002, REACH loss>.005, JAC loss>.0005 or NMAE error increase>1%; beat shuffle on profile and PDS. If failed, report the negative result, retain prior best and do not generate a renamed duplicate.',
          boundary='All recipient-study responses, errors, valid masks, knockdowns and source profiles excluded from graph and response-bank inference; query NTC profiles remain valid. Poisoning audit required. New independent biological data are not claimed; D3 reserve stays unopened unless all earlier gates pass.',
          references=['https://www.nature.com/articles/s41586-022-05688-9','https://doi.org/10.1006/jtbi.2002.3096'],
          D3_opened=False, automatic_upload=False))
    print('V67 plan frozen', sha(A/'plan.json'), flush=True)


def check():
    p = read(A/'plan.json')
    for name,h in {**p['dependencies'], **p['inputs']}.items():
        assert sha(ROOT/name)==h, name
    return p


def prepare():
    check()
    from run_v35_v36 import contexts, key
    from run_v37_ensemble import checked_folder
    for ctx in contexts():
        name = ctx['family']; path = A/'predictions'/(name+'.npz'); path.parent.mkdir(exist_ok=True)
        if path.exists():
            assert sha(path)==read(path.with_suffix('.json'))['sha256']; continue
        folder, rec = checked_folder(name)
        p66 = ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz')
        assert sha(p66)==read(p66.with_suffix('.json'))['sha256']
        with np.load(p66, allow_pickle=False) as z:
            base=z['v63'].copy(); v66=z['assay75'].copy(); q=z['query'].copy(); bulk=z['bulk'].copy(); measured=z['measured'].copy()
        effects, trace = predict(q, base, (key(name),))
        np.savez_compressed(path, query=q, bulk=bulk, measured=measured, v63=base, v66=v66, **effects)
        write(path.with_suffix('.json'), dict(sha256=sha(path), excluded_study=key(name),
              v66_predictions_sha256=sha(p66), source_receipt_sha256=sha(folder/'receipt.json'), trace=trace))
        bank.cache_clear(); old.bank.cache_clear(); gc.collect()
        print('V67 prepared', name, trace['mean_absolute_accepted_correction'], flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);args=p.parse_args()
    with threadpool_limits(limits=2):globals()[args.stage]()


