"""Gate observational neural templates on disjoint control-guide distributions."""
from neural_runtime import activate_local_runtime
activate_local_runtime()
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import torch
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from train_gene300_neural import load_network,predict_network
from train_gene300_models import split_guides
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments,de,profile_jackknife
ROOT=build.ROOT; A=ROOT/'experiments/v6_fidelity_calibrated/analysis'


def main():
    A.mkdir(parents=True,exist_ok=True)
    _,genes,targets=build.read_spec(ROOT);ipos=np.array([np.flatnonzero(genes==t)[0] for t in targets]);opos=np.flatnonzero(~np.isin(genes,targets))
    atlas,atlas_sha=load_atlas_model();records=[];identities=[]
    for ci,context in enumerate('ABC'):
        source=ad.read_h5ad(ROOT/f'context_{context}.h5ad');x=source.X.astype(np.int32).tocsr()
        roles=split_guides(source.obs.ntc_id.astype(str).to_numpy());train=np.flatnonzero(roles=='train');val=np.flatnonzero(roles=='validation');test=np.flatnonzero(roles=='test')
        model,network=load_network(ROOT/f'analysis/gene300_neural/evaluation_raw_{context}.npz')
        np.testing.assert_array_equal(model['fit_row_indices'],train)
        np.testing.assert_array_equal(model['output_genes'],genes[opos])
        inputs=x[val][:,ipos].toarray();pred=predict_network(model,network,inputs)
        depth,p,b,_=moments(x[train])
        for seed in [0,1]:
            rng=np.random.default_rng(20261007+ci*100+seed)
            local=rng.choice(len(val),400,replace=False);rows=val[local];reference=test[rng.choice(len(test),400,replace=False)]
            actual=x[rows].toarray().astype(float);lib=actual.sum(axis=1);expected=actual.copy();expected[:,opos]=pred[local]
            ref=x[reference];_,rp,rb,rv=moments(ref)
            ref_profile,_=profile_jackknife(ref)
            for name,weight in [('observed',0.),('moment_only',0.),('neural_10',.1),('neural_25',.25)]:
                template=(1-weight)*actual+weight*expected
                if name=='observed':
                    count=actual.astype(np.int32)
                else:
                    count=atlas.dual_moment_counts(template/template.sum(axis=1)[:,None],p,b,depths=lib.astype(np.int64),seed=20261007+ci*100+seed,iterations=200)
                dd,mp,mb,mv=moments(count);table=de(sparse.csr_matrix(count),ref,genes)
                gp,_=profile_jackknife(count)
                records.append({'context':context,'seed':seed,'method':name,'neural_weight':weight,
                    'null_DE_calls':int((table['q']<.05).sum()),'variance_ratio_to_reference':float(mv.sum()/rv.sum()),
                    'mean_detected_genes':float(np.count_nonzero(count,axis=1).mean()),
                    'reference_detected_genes':float(ref.getnnz(axis=1).mean()),
                    'pooled_profile_mse':float(np.mean((gp-ref_profile)**2)),
                    'mean_probability_L1':float(np.abs(mp-rp).sum()),'pooled_probability_L1':float(np.abs(mb-rb).sum())})
            identities.append({'context':context,'seed':seed,'template_ids':source.obs_names[rows].tolist(),'reference_ids':source.obs_names[reference].tolist(),
                'fit_guide_groups':sorted(source.obs.ntc_id.iloc[train].astype(str).unique().tolist()),
                'template_guide_groups':sorted(source.obs.ntc_id.iloc[rows].astype(str).unique().tolist()),
                'reference_guide_groups':sorted(source.obs.ntc_id.iloc[reference].astype(str).unique().tolist())})
        del source,x,pred,network; print(context,'neural distribution gate completed',flush=True)
    frame=pd.DataFrame(records);frame.to_csv(A/'neural_control_gate_scores.csv',index=False)
    summary=frame.groupby(['context','method']).mean(numeric_only=True).reset_index();summary.to_csv(A/'neural_control_gate_summary.csv',index=False)
    gate=[]
    for context in 'ABC':
        group=summary[summary.context==context].set_index('method');base=group.loc['observed']
        for method in ['moment_only','neural_10','neural_25']:
            row=group.loc[method]
            passed=(.85<=row.variance_ratio_to_reference<=1.15 and row.null_DE_calls<=base.null_DE_calls+10 and
                row.pooled_profile_mse<=base.pooled_profile_mse and abs(row.mean_detected_genes/base.mean_detected_genes-1)<=.1)
            gate.append({'context':context,'method':method,'control_gate_passed':bool(passed)})
    build.write_json(A/'neural_template_selection.json',{'production_neural_weight':0.,'control_gates':gate,
        'adoption_policy':'A control gate is necessary but not sufficient; no measured A/B/C perturbation labels support a neural response/template advantage. Keep empirical cells for this candidate.',
        'script_sha256':build.sha256(Path(__file__)),'upstream_sha256':atlas_sha,
        'scope':'Comparative control-guide validation using previously inspected test-guide groups. All neural weights fit train guides only; template moments use train guides; no test-guide fitting.'})
    build.write_json(A/'neural_control_identity.json',identities)
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(4)
    with threadpool_limits(limits=2):main()
