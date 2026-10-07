"""Fit the frozen V7 public-source selection on all eligible source responses."""
import json,argparse,time,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import A,ROOT,read_cache
from v7_response_methods import CONFIGS,predict_modes

def main():
    manifest,genes,targets=build.read_spec(ROOT);data=read_cache()
    np.testing.assert_array_equal(genes,data['genes']);np.testing.assert_array_equal(targets,data['targets'])
    cell=json.loads((A/'cell_selection.json').read_text());selection=json.loads((A/'summary_selection.json').read_text())
    assert build.sha256(A/'cell_summary.csv')==cell['summary_sha256']
    assert build.sha256(A/'source_cache.npz')==selection['source_cache_sha256']
    assert build.sha256(ROOT/'v7_response_methods.py')==selection['methods_sha256']
    config=next(c for c in CONFIGS if c['name']==selection['production_candidate']);train=np.ones(len(data['truth']),bool)
    profiles=[];inputs=[]
    for context in manifest['contexts']:
        x,guide,depth,_,_,stats=build.load_control(ROOT,context,genes,manifest)
        profiles.append(np.asarray(x.multiply((1e6/depth)[:,None]).mean(axis=0)).ravel());del x;gc.collect()
        inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
    first=len(data['profiles']);data['profiles']=np.r_[data['profiles'],np.asarray(profiles)]
    components_file=ROOT/'experiments/v4_lineage_response/analysis/source_components.npz'
    v4_record=json.loads(components_file.with_name('reference_model.json').read_text());assert build.sha256(components_file)==v4_record['source_components_sha256']
    with np.load(components_file,allow_pickle=False) as z:components={k:z[k].copy() for k in z.files}
    np.testing.assert_array_equal(components['genes'],genes);np.testing.assert_array_equal(components['targets'],targets)
    effect=[];details=[];candidate=cell['candidate'];neural=None;scale=None
    if candidate in ['neural_response','neural_blend25']:
        import torch
        from train_v7_transition import make_model,features,predict
        neural_rec=json.loads((A/'neural_manifest.json').read_text());assert build.sha256(A/'transition_weights.npz')==neural_rec['weights_sha256']
        sample=features(data,train,np.arange(300),np.full(300,first))
        neural=make_model(sample.shape[1],neural_rec['production_architecture'],len(genes))
        with np.load(A/'transition_weights.npz',allow_pickle=False) as z:
            scale=z['scale'].copy();neural.load_state_dict({k:torch.tensor(z[k]) for k in z.files if k!='scale'})
        neural.eval()
    for ci,context in enumerate(manifest['contexts']):
        modes,detail=predict_modes(data,train,first+ci);pool=modes[config['mode']]*config['amplitude']
        if candidate=='coverage_fill':
            old_mask=components['k562_mask']|components['cd4_mask'][ci]
            if ci==0:old=np.where(components['cd4_mask'][ci],components['cd4'][ci],np.where(components['k562_mask'],components['k562'],0))
            else:
                old_mask|=components['h1_mask'];old=np.where(components['h1_mask'],components['h1'],.5*components['cd4_specific'][ci]+components['k562'])+components['h1_prior']
            value=np.where(old_mask,old,modes['uniform'])
        elif candidate=='coverage_uniform':value=modes['uniform']
        elif candidate=='neural_response':value=predict(neural,scale,data,train,np.arange(300),np.full(300,first+ci))
        elif candidate=='neural_blend25':value=.75*pool+.25*predict(neural,scale,data,train,np.arange(300),np.full(300,first+ci))
        elif candidate=='ridge_blend25':
            ridge_rec=json.loads((A/'ridge_manifest.json').read_text());assert build.sha256(A/'ridge_production_predictions.npz')==ridge_rec['production_predictions_sha256']
            assert build.sha256(A/'source_cache.npz')==ridge_rec['source_cache_sha256']
            for r in ridge_rec['context_original_inputs']:assert build.sha256(ROOT/f'context_{r["context"]}.h5ad')==r['sha256']
            with np.load(A/'ridge_production_predictions.npz',allow_pickle=False) as z:
                supported=z['supported'];value=pool.copy();value[:,supported]=.75*pool[:,supported]+.25*z['predictions'][ci][:,supported]
        elif candidate=='v4_transfer':
            if ci==0:value=np.where(components['cd4_mask'][ci],components['cd4'][ci],np.where(components['k562_mask'],components['k562'],0))
            else:value=np.where(components['h1_mask'],components['h1'],.5*components['cd4_specific'][ci]+components['k562'])+components['h1_prior']
        else:value=pool
        seen=np.isin(np.arange(len(targets)),data['row_targets'][train]);value[~seen]=pool[~seen]
        # Preserve the prior A primary CD4 response when it is directly measured.
        # This source anchoring policy is informed by previous work, not validated A treatment labels.
        if ci==0 and candidate!='v4_transfer':value=np.where(components['cd4_mask'][0],components['cd4'][0],value)
        gate=np.clip(profiles[ci]/5,0,1)
        own=np.array([int(np.flatnonzero(genes==t)[0]) for t in targets]);target_gate=gate[own]
        value=np.clip(value,-3,3)*gate[None,:]*target_gate[:,None];value[np.arange(300),own]=0
        assert value.shape==(300,len(genes)) and np.isfinite(value).all()
        effect.append(value.astype(np.float32));details.append({'context':context,'candidate':candidate,'rank_weight_diagnostics':detail,
            'nonzero_effects':int((np.abs(value)>1e-7).sum()),'mean_absolute_log2_effect':float(np.abs(value).mean())})
    np.savez_compressed(A/'reference_model_arrays.npz',genes=genes,targets=targets,effects=np.asarray(effect))
    dependencies=['train_v7_ridge.py','embedding_response_methods.py','fit_v7_submission.py','v7_response_methods.py','v7_count_methods.py','build_submission.py','v6_methods.py']
    hashes={f:build.sha256(ROOT/f) for f in dependencies}
    build.write_json(A/'reference_model.json',{'model':'V7 expanded measured source coverage and held-lineage selected response model',
        'candidate':candidate,'generator':cell['generator'],'production_summary_config':config,'context_original_inputs':inputs,
        'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'source_cache_sha256':build.sha256(A/'source_cache.npz'),
        'summary_selection_sha256':build.sha256(A/'summary_selection.json'),'cell_selection_sha256':build.sha256(A/'cell_selection.json'),
        'source_components_sha256':build.sha256(components_file),'dependencies':hashes,'contexts':details,
        'A_policy':'Preserve directly measured V4 CD4 effects; use selected expanded-source response for unsupported coordinates. This is a biological transfer assumption, not an identified A lineage.',
        'unseen_target_guard':'Residual heads do not alter targets absent from every eligible training contrast; unseen-target transfer has not been established.',
        'gating':'Control mean CPM/5 clipped to [0,1] for output gene and own target; own response coordinate zero; count emitter applies 90% own-target knockdown.',
        'scope':'No competition treated labels; public-source development evidence does not guarantee a hidden leaderboard gain.'})
    print('V7 production model frozen:',candidate,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-for-selection',action='store_true');args=p.parse_args()
    if args.wait_for_selection:
        while not (A/'cell_selection.json').exists():time.sleep(20)
    with threadpool_limits(limits=2):main()
