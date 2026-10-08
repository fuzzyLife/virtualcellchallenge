"""V10 decoder control-to-control safety check using only provided challenge NTCs."""
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v10_methods import ROOT,A,SEEDS,emit
GENERATORS=["legacy","compact","selective"]
from v6_methods import moments,de
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build


def main():
    A.mkdir(parents=True,exist_ok=True)
    manifest,genes,targets=build.read_spec(ROOT);atlas,_=load_atlas_model();records=[];inputs=[]
    for ci,context in enumerate(manifest['contexts']):
        x,guides,library,_,_,statistics=build.load_control(ROOT,context,genes,manifest)
        rng=np.random.default_rng(SEEDS[0]);left=[];right=[]
        for guide in sorted(set(guides)):
            ids=rng.permutation(np.flatnonzero(guides==guide));half=len(ids)//2;left.extend(ids[:half]);right.extend(ids[half:])
        assert not set(left)&set(right)
        fit=x[left];ref=x[right];left_guides=guides[left];right_guides=guides[right]
        depth,p,b,var=moments(fit);prepared=(depth,p,b,np.asarray(fit.mean(axis=0)).ravel())
        for seed in SEEDS:
            rows=build.guide_stratified_sample(left_guides,library[left],400,np.random.default_rng(seed+ci))
            reference_rows=build.guide_stratified_sample(right_guides,library[right],400,np.random.default_rng(seed+ci+2000))
            control=ref[reference_rows]
            for generator in GENERATORS:
                count,diag=emit(fit,rows,np.zeros(len(genes)),None,seed+ci,generator,prepared=prepared,atlas=atlas)
                metrics=de(count,control,genes);false=(metrics['q']<.05)&metrics['eligible']
                libs,pi,bi,vi=moments(count)
                records.append({'context':context,'seed':seed,'generator':generator,
                    'tested_genes':int(metrics['eligible'].sum()),'false_de_calls':int(false.sum()),
                    'false_call_fraction':float(false.sum()/max(metrics['eligible'].sum(),1)),
                    'CPM_variance_ratio':float(vi.sum()/var.sum()),
                    'mean_probability_L1_from_population':float(abs(pi-p).sum()),
                    'bulk_probability_L1_from_population':float(abs(bi-b).sum()),**diag})
                print('V10 neutral',context,seed,generator,'false DE',int(false.sum()),flush=True)
            pd.DataFrame(records).to_csv(A/'decoder_neutral_scores.csv',index=False)
        inputs.append({'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
        del x,fit,ref;gc.collect()
    frame=pd.DataFrame(records)
    # An explicit control-preservation gate, fixed before any confirmation scoring.
    summary=frame.groupby('generator').agg(max_false_de_calls=('false_de_calls','max'),
        mean_false_de_calls=('false_de_calls','mean'),max_false_call_fraction=('false_call_fraction','max'),
        mean_variance_ratio=('CPM_variance_ratio','mean')).reset_index()
    summary['neutral_guard']=summary.max_false_call_fraction<=.005
    summary.to_csv(A/'decoder_neutral_summary.csv',index=False)
    build.write_json(A/'decoder_neutral_manifest.json',{'passed_generators':summary.loc[summary.neutral_guard,'generator'].tolist(),
        'criterion':'Diagnostic only for this distribution-model experiment; not a selection gate. Every context/seed neutral test has BH-significant false calls <=0.5% of tested genes. Genuine disjoint provided NTCs, with no target knockdown.',
        'source_inputs':inputs,'scores_sha256':build.sha256(A/'decoder_neutral_scores.csv'),
        'summary_sha256':build.sha256(A/'decoder_neutral_summary.csv'),'script_sha256':build.sha256(ROOT/'benchmark_v10_neutral.py'),
        'methods_sha256':build.sha256(ROOT/'v10_methods.py'),'scope':'Control-only population-model bias diagnostic; does not evaluate hidden perturbation accuracy.'})
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
