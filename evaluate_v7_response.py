"""Compare source-held response models using 400 generated cells and actual Wilcoxon/BH."""
from pathlib import Path
import argparse,json,time,hashlib,gc
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import A,ROOT,raw_inputs,read_cache
from v7_response_methods import CONFIGS
from v7_count_methods import emit
from v6_methods import moments,de,direction_metrics,profile_jackknife,summarize_profiles
from evaluate_atlas_transfer import load_atlas_model
# Previous committed numerical implementation; only missing-metric JSON encoding changed.
SERIALIZATION_COMPATIBLE_SHAS={'c320ed8e8d106864f0a420113fab0888ac0402e5c94549ede7106c026e1a922a'}
SEEDS=[20261008,20261019]
CANDIDATES={'v4_transfer':'v4','coverage_fill':'v4','coverage_uniform':'v4','calibrated_pool':'v4','pool_low_noise':'low_noise','pool_cal10':'cal10','neural_response':'v4','neural_blend25':'v4','ridge_blend25':'v4'}

def baseline(line):
    directory=ROOT/'experiments/v4_lineage_response/analysis'
    record=json.loads((directory/'reference_model.json').read_text());assert build.sha256(directory/'source_components.npz')==record['source_components_sha256']
    _,genes,targets=build.read_spec(ROOT)
    with np.load(directory/'source_components.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
        setting=record['cross_lineage_settings'];transfer=setting['cd4_coefficient']*z['cd4_specific'][2]+setting['k562_coefficient']*z['k562']
        return transfer if line=='H1' else np.where(z['h1_mask'],z['h1'],transfer)+z['h1_prior']

def matched_rows(meta,available,treated,rng,n=400):
    """Sample distinct controls in the same batch proportions as treated cells."""
    batches=meta.iloc[treated]['sample'].value_counts();desired=n*batches/batches.sum()
    quotas=np.floor(desired).astype(int)
    order=(desired-quotas).sort_values(ascending=False).index
    for b in order[:n-int(quotas.sum())]:quotas[b]+=1
    result=[]
    for batch,count in quotas.items():
        pool=available[meta.iloc[available]['sample'].to_numpy()==batch]
        if len(pool)<count:raise ValueError('Insufficient distinct matched controls')
        result.extend(rng.choice(pool,int(count),replace=False).tolist())
    return np.asarray(result,dtype=int)

def population_statistics(x,meta,available):
    """Control-population moments per batch, separate from the 400 sampled templates."""
    labels=meta.iloc[available]['sample'];batches=sorted(labels.unique());lookup={b:i for i,b in enumerate(batches)}
    group=np.array([lookup[b] for b in labels]);count=np.bincount(group,minlength=len(batches))
    block=x[available].astype(float);depth=np.asarray(block.sum(axis=1)).ravel()
    assignment=sparse.csr_matrix((1/count[group],(group,np.arange(len(available)))),shape=(len(batches),len(available)))
    raw=(assignment@block).toarray();mean=(assignment@block.multiply((1/depth)[:,None])).toarray()
    return batches,raw,mean

def matched_population(population,meta,actual):
    batches,raw,mean=population;counts=meta.iloc[actual]['sample'].value_counts()
    weight=np.array([counts.get(b,0)/len(actual) for b in batches]);assert abs(weight.sum()-1)<1e-10
    mean_raw=weight@raw;p=weight@mean;library=float(mean_raw.sum());b=mean_raw/library
    return np.array([library]),p,b,mean_raw

def models(line,ci):
    selection=json.loads((A/'summary_selection.json').read_text())
    config=next(c for c in CONFIGS if c['name']==selection['outer_selections'][line]['candidate'])
    with np.load(A/'outer_predictions.npz',allow_pickle=False) as z:
        pool=z[f'{line}_{ci}_{config["mode"]}']*config['amplitude'];uniform=z[f'{line}_{ci}_uniform']
    with np.load(A/'neural_outer_predictions.npz',allow_pickle=False) as z:neural=z[f'{line}_{ci}_nn_nested'].copy()
    data=read_cache();seen=np.isin(np.arange(len(data['targets'])),data['row_targets'][data['row_families']!=line]);neural[~seen]=0
    with np.load(A/'ridge_outer_predictions.npz',allow_pickle=False) as z:
        ridge=z[f'{line}_{ci}_ridge'];supported=z[f'{line}_{ci}_supported'];ridge_blend=pool.copy();ridge_blend[:,supported]=.75*pool[:,supported]+.25*ridge[:,supported];ridge_blend[~seen]=pool[~seen]
    base=baseline(line)
    with np.load(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz',allow_pickle=False) as z:
        old_mask=z['k562_mask']|z['cd4_mask'][2]
        if line!='H1':old_mask|=z['h1_mask']
    coverage=np.where(old_mask,base,uniform)
    return {'v4_transfer':base,'coverage_fill':coverage,'coverage_uniform':uniform,'calibrated_pool':pool,
        'pool_low_noise':pool,'pool_cal10':pool,'neural_response':neural,'neural_blend25':.75*pool+.25*neural,'ridge_blend25':ridge_blend}

def evaluate_line(line,ci,max_targets):
    d=read_cache();genes=d['genes'];targets=d['targets'];x,meta,measured,source=raw_inputs(line,genes)
    response=models(line,ci);atlas,_=load_atlas_model();controls=np.flatnonzero(meta.gene_target.to_numpy()=='Non-Targeting')
    rng=np.random.default_rng(SEEDS[0]);left=[];right=[]
    for batch,group in meta.iloc[controls].groupby('sample',sort=True):
        ids=rng.permutation(group.index.to_numpy());half=len(ids)//2;left.extend(ids[:half]);right.extend(ids[half:])
    left=np.asarray(left);right=np.asarray(right);assert not set(left)&set(right)
    library=np.asarray(x.sum(axis=1)).ravel();population=population_statistics(x,meta,left)
    query=np.asarray(x[left].multiply((1e6/library[left])[:,None]).mean(axis=0)).ravel()
    gate=np.clip(query/5,0,1);eligible=[];ids_record=[];records=[]
    directory=A/'cell_evaluation';directory.mkdir(exist_ok=True)
    files=['source_cache.npz','outer_predictions.npz','neural_outer_predictions.npz','summary_selection.json','ridge_outer_predictions.npz','ridge_manifest.json']
    input_hash={f:build.sha256(A/f) for f in files}
    input_hash.update({f:build.sha256(ROOT/f) for f in ['evaluate_v7_response.py','v7_count_methods.py','v6_methods.py']})
    order=sorted(range(len(targets)),key=lambda ti:hashlib.sha256(str(targets[ti]).encode()).hexdigest())
    for ti in order:
        treated=np.flatnonzero(meta.gene_target.to_numpy()==targets[ti]);own=int(np.flatnonzero(genes==targets[ti])[0])
        if len(treated)<200 or query[own]<=5:continue
        # Verify batch-matched disjoint control support before selecting the target.
        check=np.random.default_rng(SEEDS[0]+ti)
        actual=check.choice(treated,min(400,len(treated)),replace=False)
        try:matched_rows(meta,left,actual,check);matched_rows(meta,right,actual,check)
        except ValueError:continue
        eligible.append(ti)
        if len(eligible)>=max_targets:break
    if len(eligible)<4:raise ValueError(f'{line}: fewer than four supported targets for 400-cell prediction ({len(eligible)})')
    for index,ti in enumerate(eligible):
        target=targets[ti];own=int(np.flatnonzero(genes==target)[0]);split='calibration' if index%2==0 else 'confirmation'
        for seed in SEEDS:
            file=directory/f'{line}_{ti:03d}_{seed}.npz';side=file.with_suffix('.json')
            if file.exists() and side.exists():
                rec=json.loads(side.read_text())
                compatible=rec['inputs']==input_hash or (rec['inputs'].get('evaluate_v7_response.py') in SERIALIZATION_COMPATIBLE_SHAS and {k:v for k,v in rec['inputs'].items() if k!='evaluate_v7_response.py'}=={k:v for k,v in input_hash.items() if k!='evaluate_v7_response.py'})
                if compatible and rec['sha256']==build.sha256(file):
                    if rec['inputs']!=input_hash:
                        rec['serialization_migration']={'previous_script_sha256':rec['inputs']['evaluate_v7_response.py'],'change':'NaN metrics serialized as null only; numerical formulas, models, candidates, cells and seeds unchanged.'};rec['inputs']=input_hash;build.write_json(side,rec)
                    records.extend(rec['records']);ids_record.append(rec['cells']);continue
                raise ValueError('Existing evaluation differs from frozen inputs')
            rng=np.random.default_rng(seed+ti);treated=np.flatnonzero(meta.gene_target.to_numpy()==target)
            actual=rng.choice(treated,min(400,len(treated)),replace=False)
            try:template=matched_rows(meta,left,actual,rng);reference=matched_rows(meta,right,actual,rng)
            except ValueError:raise ValueError(f'{line} {target}: seed lacks distinct matched controls; preserve failed selection')
            truth=x[actual];control=x[reference];templates=x[template]
            realde=de(truth,control,genes);tp,tj=profile_jackknife(truth);cp,cj=profile_jackknife(control)
            arrays={'truth':tp,'truth_jack':tj,'ctrl':cp,'ctrl_jack':cj};group=[]
            prepared_m=moments(templates);prepared=matched_population(population,meta,actual)
            truth_var=moments(truth)[3].sum();original_var=prepared_m[3].sum()
            for name,generator in CANDIDATES.items():
                effect=response[name][ti]*gate*np.clip(query[own]/5,0,1);effect=effect.copy();effect[own]=0
                prediction=emit(templates,np.arange(400),effect,own,seed+ti,generator,atlas,prepared)
                metrics=direction_metrics(de(prediction,control,genes),realde,genes,own)
                pp,pj=profile_jackknife(prediction);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj
                pred_var=moments(prediction)[3].sum()
                row={'family':line,'target_gene':str(target),'target_index':ti,'split':split,'seed':seed,'candidate':name,
                    'generator':generator,'variance_ratio':float(pred_var/max(truth_var,1e-30)),
                    'template_variance_ratio':float(pred_var/max(original_var,1e-30)),
                    'observed_cells':len(actual),'generated_cells':400,'reference_cells':400,'mean_library':float(np.asarray(prediction.sum(axis=1)).mean()),**metrics}
                group.append(row);del prediction
            cells={'target_gene':str(target),'seed':seed,'split':split,'actual':meta.iloc[actual].cell_id.tolist(),
                'template':meta.iloc[template].cell_id.tolist(),'reference':meta.iloc[reference].cell_id.tolist()}
            assert not set(cells['template'])&set(cells['reference'])
            np.savez_compressed(file,**arrays)
            build.write_json(side,{'inputs':input_hash,'sha256':build.sha256(file),'records':[{k:None if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in r.items()} for r in group],'cells':cells})
            records.extend(group);ids_record.append(cells)
        pd.DataFrame(records).to_csv(A/f'cell_scores_{line}.csv',index=False)
        print('Cell evaluation',line,index+1,'/',len(eligible),target,flush=True)
    frame=pd.DataFrame(records);summaries=[]
    for split in ['calibration','confirmation','all']:
        selected=eligible if split=='all' else eligible[0::2] if split=='calibration' else eligible[1::2]
        for seed in SEEDS:
            arrays=[dict(np.load(directory/f'{line}_{ti:03d}_{seed}.npz',allow_pickle=False)) for ti in selected]
            for name in CANDIDATES:
                r=frame[(frame.seed==seed)&(frame.candidate==name)&frame.target_index.isin(selected)]
                stats=r.mean(numeric_only=True).to_dict()
                profile=summarize_profiles(np.array([a[name+'_pred'] for a in arrays]),np.array([a['truth'] for a in arrays]),
                    np.array([a['ctrl'] for a in arrays]),np.array([a[name+'_jack'] for a in arrays]),
                    np.array([a['truth_jack'] for a in arrays]),np.array([a['ctrl_jack'] for a in arrays]),genes,targets[selected],targets)
                summaries.append({'family':line,'candidate':name,'split':split,'seed':seed,'targets':len(selected),**stats,**profile})
    pd.DataFrame(summaries).to_csv(A/f'cell_summary_{line}.csv',index=False)
    build.write_json(A/f'cell_split_{line}.json',{'source':source,'inputs':input_hash,'excluded_training_family':line,
        'targets':targets[eligible].tolist(),'cells':ids_record,'selection':'Target-ID hash order; >=200 distinct treated cells (up to 400 observed) and expressed own target; exactly 400 distinct generated templates and reference cells. Batch-matched, distinct template/reference controls; no response-based target selection.',
        'scope':'Public development; calibration and confirmation target sets disjoint; previously inspected data, not a pristine blind evaluation.'})
    del x;gc.collect()

def select():
    f=pd.concat([pd.read_csv(A/f'cell_summary_{line}.csv') for line in ['H1','HCT116','HEK293T']],ignore_index=True)
    f['utility']=(f.pds+f.fidelity+f.reach+f.jaccard+1-np.minimum(f.nmae,2)+1-np.minimum(f.expression_ratio,2))/6
    f.to_csv(A/'cell_summary.csv',index=False)
    means=f.groupby(['split','family','candidate']).mean(numeric_only=True).groupby(['split','candidate']).mean()
    rows=[]
    for name in CANDIDATES:
        row={'candidate':name}
        for split in ['calibration','confirmation']:
            a=means.loc[(split,name)];b=means.loc[(split,'v4_transfer')]
            row[split+'_utility']=float(a.utility)
            row[split+'_guard']=bool(a.fidelity>=b.fidelity-.01 and a.expression_ratio<=1.05*b.expression_ratio and a.template_variance_ratio>=.9*b.template_variance_ratio)
        rows.append(row)
    tab=pd.DataFrame(rows);cal=tab[tab.calibration_guard].sort_values('calibration_utility',ascending=False);winner=cal.iloc[0].candidate
    chosen=tab[tab.candidate==winner].iloc[0];base=tab[tab.candidate=='v4_transfer'].iloc[0]
    confirmed=bool(chosen.confirmation_guard and chosen.confirmation_utility>base.confirmation_utility)
    selected=winner if confirmed else 'v4_transfer';tab.to_csv(A/'cell_selection_scores.csv',index=False)
    build.write_json(A/'cell_selection.json',{'candidate':selected,'calibration_winner':winner,'confirmed_improvement':confirmed,
        'generator':CANDIDATES[selected],'summary_sha256':build.sha256(A/'cell_summary.csv'),
        'utility':'Equal mean of raw PDS, fidelity, reach, Jaccard, clipped 1-NMAE and 1-expression-ratio; equal weight to source families. Not leaderboard anchor normalization.',
        'guards':'FID within .01 of baseline, expression ratio <=1.05 baseline, template variance >=.9 baseline. Freeze calibration winner before inspecting confirmation; fall back if not confirmed.',
        'scope':'Previously inspected public development, source-held models, disjoint target calibration/confirmation, two seeds. No hidden score estimate.'})
    print('Selected generated-cell candidate',selected,'confirmed',confirmed,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-for-models',action='store_true');p.add_argument('--select-only',action='store_true');p.add_argument('--line',choices=['H1','HCT116','HEK293T','all'],default='all');p.add_argument('--targets',type=int,default=18);a=p.parse_args()
    if a.wait_for_models or not a.select_only:
        while not all((A/f).exists() for f in ['summary_selection.json','neural_manifest.json','ridge_manifest.json']):time.sleep(20)
    with threadpool_limits(limits=2):
        if not a.select_only:
            for line,ci in [('H1',4),('HCT116',5),('HEK293T',6)]:
                if a.line in ['all',line]:evaluate_line(line,ci,a.targets)
        if a.line=='all' or a.select_only:select()
