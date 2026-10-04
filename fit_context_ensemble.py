"""Fit a context-specific response ensemble with independent-source validation."""
from pathlib import Path
import argparse
import hashlib
import json
import gc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls
import build_submission as build
from build_cd4_model import load_subset, eligible, align_disjoint

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v3_context_ensemble/analysis'
CONDITIONS=['Rest','Stim8hr','Stim48hr']
REFERENCES=['K562','HCT116','HEK293T','H1','CD4']


def target_folds(targets, n=5):
    unique=sorted(set(targets))
    order=sorted(unique,key=lambda t:hashlib.sha256(t.encode()).hexdigest())
    mapping={target:i%n for i,target in enumerate(order)}
    return np.array([mapping[t] for t in targets])


def prior_mean(values, errors, valid, conditions, train):
    result={}
    for condition in CONDITIONS:
        rows=train & (conditions==condition)
        weights=valid[rows]/(errors[rows]**2+.04)
        denominator=weights.sum(axis=0)
        result[condition]=np.clip(np.divide((weights*values[rows]).sum(axis=0),denominator,
            out=np.zeros(values.shape[1]),where=denominator>0),-.5,.5)
    return result


def estimate(values, errors, q, conditions, settings, prior):
    if settings['kind']=='zero':
        return np.zeros_like(values)
    if settings['kind']=='hard':
        return settings['alpha']*np.where(q<settings['q'],values,0)
    mu=np.stack([prior[c] for c in conditions]) if settings['prior'] else np.zeros_like(values)
    weight=settings['tau2']/(settings['tau2']+errors**2)
    return np.clip(settings['alpha']*(mu+weight*(values-mu)),-3,3)


def loss(pred_left,pred_right,left,right,valid,rows):
    denom=np.maximum(valid[rows].sum(axis=1),1)
    values=(((pred_left[rows]-right[rows])**2+(pred_right[rows]-left[rows])**2)*valid[rows]).sum(axis=1)/denom
    return float(values.mean())


def candidates():
    result=[{'kind':'zero','alpha':0.}]
    result.extend({'kind':'hard','q':q,'alpha':a} for q in [.05,.2] for a in [.3,.7,1.])
    result.extend({'kind':'eb','tau2':tau,'alpha':a,'prior':prior} for tau in [.0025,.01,.04] for a in [.5,1.] for prior in [False,True])
    return result


def validate_cd4(analysis):
    directory=ROOT/'reference/cd4_probe'
    source_report=json.loads((ROOT/'analysis/model_audit/cd4_extraction.json').read_text())
    for record in source_report:
        if build.sha256(directory/record['subset'])!=record['sha256']:
            raise ValueError('CD4 source subset hash mismatch')
    names=[Path(r['subset']).stem for r in source_report if Path(r['subset']).stem!='merged']
    if len(names)!=2 or set(names[0].split('_')) & set(names[1].split('_')):
        raise ValueError('Donor groups must be disjoint')
    lmeta,ld=load_subset(directory,names[0]); rmeta,rd=load_subset(directory,names[1])
    rows,left,right,lq,rq,valid=align_disjoint((lmeta,ld),(rmeta,rd))
    symbols=np.array(sorted(set(ld['genes'])&set(rd['genes'])))
    lm={g:i for i,g in enumerate(ld['genes'])}; rm={g:i for i,g in enumerate(rd['genes'])}
    le=ld['lfcSE'][rows.position_left.to_numpy()][:,[lm[g] for g in symbols]].astype(float)
    re=rd['lfcSE'][rows.position_right.to_numpy()][:,[rm[g] for g in symbols]].astype(float)
    valid &= np.isfinite(le)&np.isfinite(re)&(le>=0)&(re>=0)
    le=np.nan_to_num(le,nan=1e6,posinf=1e6,neginf=1e6)
    re=np.nan_to_num(re,nan=1e6,posinf=1e6,neginf=1e6)
    conditions=rows.culture_condition.to_numpy()
    targets=rows.target_contrast_gene_name_left.to_numpy()
    folds=target_folds(targets)
    grid=candidates(); candidate_scores=[]; nested=[]
    baseline=loss(np.zeros_like(left),np.zeros_like(right),left,right,valid,np.ones(len(rows),bool))
    for fold in range(5):
        test=folds==fold; train=~test
        inner=target_folds(targets[train],3)
        internal=[[] for _ in grid]
        train_indices=np.flatnonzero(train)
        for part in range(3):
            tune=np.zeros(len(rows),bool); tune[train_indices[inner==part]]=True
            fit=train & ~tune
            lp=prior_mean(left,le,valid,conditions,fit); rp=prior_mean(right,re,valid,conditions,fit)
            for ci,settings in enumerate(grid):
                pl=estimate(left,le,lq,conditions,settings,lp)
                pr=estimate(right,re,rq,conditions,settings,rp)
                internal[ci].append((int(tune.sum()),loss(pl,pr,left,right,valid,tune)))
        internal_scores=[sum(n*v for n,v in group)/sum(n for n,_ in group) for group in internal]
        selected=int(np.argmin(internal_scores))
        lp=prior_mean(left,le,valid,conditions,train); rp=prior_mean(right,re,valid,conditions,train)
        for ci,settings in enumerate(grid):
            pl=estimate(left,le,lq,conditions,settings,lp); pr=estimate(right,re,rq,conditions,settings,rp)
            error=loss(pl,pr,left,right,valid,test)
            candidate_scores.append({'fold':fold,'candidate':ci,'settings':json.dumps(settings,sort_keys=True),'contrasts':int(test.sum()),'loss':error,'relative_zero_loss':error/baseline})
            if ci==selected:
                nested.append({'fold':fold,'candidate':ci,'settings':settings,'contrasts':int(test.sum()),'loss':error,
                    'test_targets':sorted(set(targets[test])),'train_targets':sorted(set(targets[train]))})
        print('CD4 nested fold',fold,'selected',grid[selected],flush=True)
    score=pd.DataFrame(candidate_scores)
    score.to_csv(analysis/'cd4_candidate_outer_scores.csv',index=False)
    means=score.groupby('candidate').apply(lambda x:np.average(x.loss,weights=x.contrasts),include_groups=False)
    selected=int(means.idxmin())
    # Never accept a method worse than the existing hard-threshold candidate in outer folds.
    existing=next(i for i,s in enumerate(grid) if s=={'kind':'hard','q':.2,'alpha':.3})
    result={'donor_groups':names,'contrasts':len(rows),'targets':len(set(targets)),'genes':len(symbols),'selected':grid[selected],
        'selected_candidate':selected,'selected_outer_cv_relative_zero_loss':float(means[selected]/baseline),
        'v2_outer_cv_relative_zero_loss':float(means[existing]/baseline),
        'nested_selection_relative_zero_loss':float(sum(x['contrasts']*x['loss'] for x in nested)/len(rows)/baseline),
        'criterion':'Equal target-condition, bidirectional clipped log2FC squared error; own target excluded; 5 target-grouped outer folds, 3 inner folds',
        'nested_folds':nested,'limitations':['CD4 donor validation only; no hidden VCC labels.', 'Source donor subsets are enriched for larger responses.',
        'Final settings use all outer-fold results; their selected minimum is optimistic. Nested selection error evaluates the selection procedure.']}
    build.write_json(analysis/'cd4_model_selection.json',result)
    print('CD4 validation',json.dumps({k:v for k,v in result.items() if k!='nested_folds'}),flush=True)
    return grid[selected]


def fit_cd4(settings, genes,targets,analysis):
    meta,data=load_subset(ROOT/'reference/cd4_probe','merged')
    source=data['genes']; good=eligible(meta)
    conditions=meta.culture_condition.to_numpy()
    values=np.clip(np.nan_to_num(data['log_fc'].astype(float),nan=0,posinf=0,neginf=0),-3,3)
    errors=np.nan_to_num(data['lfcSE'].astype(float),nan=1e6,posinf=1e6,neginf=1e6)
    valid=np.isfinite(data['log_fc'])&np.isfinite(data['lfcSE'])&(data['lfcSE']>=0)
    valid &= source[None,:]!=meta.target_contrast_gene_name.to_numpy()[:,None]
    prior=prior_mean(values,errors,valid,conditions,good)
    fitted=estimate(values,errors,data['adj_p_value'],conditions,settings,prior)
    fitted[~np.isfinite(data['log_fc'])]=0
    position={g:i for i,g in enumerate(genes)}
    si=np.array([i for i,g in enumerate(source) if g in position]); di=np.array([position[source[i]] for i in si])
    effects=np.zeros((3,len(targets),len(genes)),dtype=np.float32)
    available=np.zeros((3,len(targets)),bool)
    for ci,condition in enumerate(CONDITIONS):
        for ti,target in enumerate(targets):
            rows=np.flatnonzero(good & (conditions==condition)&(meta.target_contrast_gene_name.to_numpy()==target))
            if len(rows):
                weights=np.sqrt(meta.n_cells_target.to_numpy()[rows])
                effects[ci,ti,di]=np.average(fitted[rows][:,si],axis=0,weights=weights)
                available[ci,ti]=True
    build.write_json(analysis/'cd4_condition_coverage.json',{'conditions':dict(zip(CONDITIONS,available.sum(axis=1).astype(int).tolist())),
        'measured_genes':len(di),'targets_any_condition':int(available.any(axis=0).sum()),'settings':settings})
    measured=np.zeros(len(genes),bool); measured[di]=True
    return effects,available,measured


def mixture_fit(profiles, query):
    scale=1/np.sqrt(np.maximum(profiles.mean(axis=1),1e-6))
    coefficients=nnls(profiles*scale[:,None],query*scale)[0]
    weights=coefficients/coefficients.sum() if coefficients.sum()>0 else np.full(profiles.shape[1],1/profiles.shape[1])
    return coefficients,weights,scale


def context_mixtures(genes,analysis):
    directory=ROOT/'reference/control_profiles'
    references={}
    for name in REFERENCES+['CD4_'+c for c in CONDITIONS]:
        path=directory/(name+'.csv')
        if build.sha256(path)!=json.loads(path.with_suffix('.json').read_text())['profile_sha256']:
            raise ValueError('Basal reference hash mismatch')
        references[name]=pd.read_csv(path).set_index('gene_name').control_count_sum
    targets=set(pd.read_csv(ROOT/'pert_counts.csv').target_gene)
    records=[]; bootstraps=[]; holdouts=[]; results={}
    rng=np.random.default_rng(20261005)
    for context in ['A','B','C']:
        stats=pd.read_csv(ROOT/f'analysis/gene_statistics_{context}.csv').set_index('gene_name').mean_raw_count
        results[context]={}
        for kind,names in [('cell_profiles',REFERENCES),('cd4_conditions',['CD4_'+c for c in CONDITIONS])]:
            axis=sorted((set(stats.index)&set.intersection(*(set(references[n].index) for n in names)))-targets)
            matrix=np.stack([references[n].reindex(axis).to_numpy() for n in names],axis=1)
            matrix=matrix/matrix.sum(axis=0)[None,:]*1e6
            query=stats.reindex(axis).to_numpy(); query=query/query.sum()*1e6
            # A gene absent from every source basis cannot inform mixture coefficients.
            expressed=matrix.mean(axis=1)>5
            matrix=matrix[expressed]; query=query[expressed]; selected=np.array(axis)[expressed]
            coefficient,weights,scale=mixture_fit(matrix,query)
            for fold in range(5):
                test=np.array([int(hashlib.sha256(g.encode()).hexdigest(),16)%5==fold for g in selected])
                coeff,_,_=mixture_fit(matrix[~test],query[~test])
                errors=[np.mean(((query[test]-matrix[test,j])*scale[test])**2) for j in range(len(names))]
                mixed=np.mean(((query[test]-matrix[test]@coeff)*scale[test])**2)
                holdouts.append({'context':context,'basis':kind,'fold':fold,'genes_tested':int(test.sum()),'mixture_weighted_mse':mixed,'best_single_profile_weighted_mse':min(errors),'relative_single_error':mixed/min(errors)})
            for replicate in range(100):
                rows=rng.integers(0,len(matrix),len(matrix))
                _,bw,_=mixture_fit(matrix[rows],query[rows])
                bootstraps.append({'context':context,'basis':kind,'replicate':replicate,**dict(zip(names,bw))})
            results[context][kind]=dict(zip(names,weights.tolist()))
            for name,w in zip(names,weights):
                records.append({'context':context,'basis':kind,'profile':name,'weight':w,'fitted_genes':len(matrix),
                    'residual_weighted_mse':float(np.mean(((query-matrix@coefficient)*scale)**2)),
                    'query_cpm_fraction_on_fitted_genes':float(query.sum()/1e6),
                    'interpretation':'Expression-fit weights, not biological cell fractions'})
    pd.DataFrame(records).to_csv(analysis/'context_mixture_weights.csv',index=False)
    pd.DataFrame(bootstraps).to_csv(analysis/'mixture_gene_bootstraps.csv',index=False)
    pd.DataFrame(holdouts).to_csv(analysis/'mixture_heldout_gene_errors.csv',index=False)
    build.write_json(analysis/'context_mixture_weights.json',results)
    print('Context mixture weights',json.dumps(results),flush=True)
    return results


def h1_response_model(genes,targets,analysis):
    directory=ROOT/'reference/h1_response'
    record=json.loads((directory/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:
        if build.sha256(directory/name)!=record[key]:
            raise ValueError('H1 extracted subset hash mismatch')
    matrix=sparse.load_npz(directory/'selected_counts.npz').astype(float)
    meta=pd.read_csv(directory/'selected_cells.csv')
    source=pd.read_csv(directory/'gene_names.csv').gene_name.to_numpy()
    lib=np.asarray(matrix.sum(axis=1)).ravel()
    cpm=matrix.multiply((1e6/lib)[:,None]).tocsr()
    batches=sorted(meta.batch.unique()); partitions={batch:i%2 for i,batch in enumerate(batches)}
    split=meta.batch.map(partitions).to_numpy()
    controls=meta.target_gene.to_numpy()=='non-targeting'
    target_meta=meta.target_gene.to_numpy()
    control_groups={}
    for batch in batches:
        rows=np.flatnonzero(controls & (meta.batch.to_numpy()==batch))
        sub=cpm[rows]
        mean=np.asarray(sub.mean(axis=0)).ravel()
        variance=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-mean**2,0)
        control_groups[batch]=(mean,variance/max(len(rows)-1,1),len(rows))
    def response(target,part):
        mask=(target_meta==target)
        if part is not None:
            mask &= split==part
        rows=np.flatnonzero(mask); n=len(rows)
        if n<20:
            return None
        sub=cpm[rows]; mean=np.asarray(sub.mean(axis=0)).ravel()
        varmean=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-mean**2,0)/max(n-1,1)
        counts=meta.iloc[rows].batch.value_counts()
        baseline=np.zeros(len(source)); basevar=np.zeros(len(source))
        for batch,count in counts.items():
            bm,bv,_=control_groups[batch]
            baseline+=count/n*bm; basevar+=(count/n)**2*bv
        lfc=np.clip(np.log2((mean+1)/(baseline+1)),-3,3)
        se=np.sqrt(varmean/(mean+1)**2+basevar/(baseline+1)**2)/np.log(2)
        return lfc,se,baseline,n
    shared_targets=sorted(set(target_meta)-{'non-targeting'})
    validation=[]; full={}; pairs=[]
    for target in shared_targets:
        a=response(target,0); b=response(target,1); full[target]=response(target,None)
        if a is None or b is None:
            continue
        valid=(a[2]>5)&(b[2]>5)&(source!=target)
        pairs.append((target,a,b,valid))
    grid=[{'alpha':a,'noise_z':z} for a in [0.,.3,.5,.7,1.] for z in [0.,1.,1.96]]
    def predict(value,setting):
        lfc,se,_,_=value
        return setting['alpha']*np.sign(lfc)*np.maximum(np.abs(lfc)-setting['noise_z']*se,0)
    def error(subset,setting):
        values=[]; anchors=[]
        for _,a,b,mask in subset:
            pa=predict(a,setting); pb=predict(b,setting)
            values.append(np.mean(((pa-b[0])**2+(pb-a[0])**2)[mask]))
            anchors.append(np.mean((a[0]**2+b[0]**2)[mask]))
        return float(np.mean(values)),float(np.mean(anchors))
    if len(pairs)<3:
        raise ValueError('Insufficient independent H1 target/batch contrasts')
    for setting in grid:
        e,anchor=error(pairs,setting)
        validation.append({**setting,'relative_zero_loss':e/anchor,'targets':len(pairs)})
    best=min(validation,key=lambda x:(x['relative_zero_loss'],x['alpha'],x['noise_z']))
    nested=[]
    for held in pairs:
        train=[p for p in pairs if p[0]!=held[0]]
        chosen=min(grid,key=lambda s:error(train,s)[0])
        e,anchor=error([held],chosen)
        nested.append({'target_gene':held[0],**chosen,'loss':e,'zero_loss':anchor})
    pd.DataFrame(validation).to_csv(analysis/'h1_batch_validation_grid.csv',index=False)
    pd.DataFrame(nested).to_csv(analysis/'h1_target_holdout_validation.csv',index=False)
    effects=np.zeros((len(targets),len(genes)),dtype=np.float32); available=np.zeros(len(targets),bool)
    positions={g:i for i,g in enumerate(genes)}
    si=np.array([i for i,g in enumerate(source) if g in positions]); di=np.array([positions[source[i]] for i in si])
    target_pos={t:i for i,t in enumerate(targets)}
    coverage=[]
    gene_available=np.zeros((len(targets),len(genes)),bool)
    nested_ratio=sum(row['loss'] for row in nested)/sum(row['zero_loss'] for row in nested)
    use=best['alpha']>0 and nested_ratio<1
    for target,value in full.items():
        lfc,se,base,n=value
        own=np.flatnonzero(source==target)
        own_lfc=float(lfc[own[0]]) if len(own) else None
        accept=use and n>=40 and own_lfc is not None and own_lfc<=np.log2(.7)
        ti=target_pos[target]
        if accept:
            effect=predict(value,best)
            # Unsupported low-expression genes are not treated as measured response estimates.
            effect[base<=5]=0
            effects[ti,di]=effect[si]; available[ti]=True
            gene_available[ti,di]=base[si]>5
        coverage.append({'target_gene':target,'source_cells':n,'own_log2FC':own_lfc,'accepted':bool(accept)})
    pd.DataFrame(coverage).to_csv(analysis/'h1_response_coverage.csv',index=False)
    build.write_json(analysis/'h1_model_selection.json',{'selected':best,'nested_target_holdout_relative_zero_loss':nested_ratio,
        'accepted_targets':int(available.sum()),'source_targets':len(full),'validation_targets':len(pairs),
        'scope':'Disjoint public 2025 source batches, not challenge or cross-cell-type validation','source':record,
        'gene_response_definition':'Batch-matched mean CPM log2 ratio with 1-CPM continuity value; sampling SE soft threshold; genes with control CPM <=5 excluded',
        'quality':'At least 40 selected cells and >=30% measured source on-target reduction; requires nonzero settings and held-target selection error below zero response'})
    print('H1 settings',best,'nested error',nested_ratio,'accepted',int(available.sum()),flush=True)
    return effects,available,gene_available


def combine_sources(cd4,cd4_available,k562,k562_available,h1,h1_available,mixtures,cd4_gene_mask=None,k562_gene_mask=None,h1_gene_mask=None):
    """Mix available response sources; never invent responses for unavailable types."""
    output=np.zeros((3,cd4.shape[1],cd4.shape[2]),dtype=np.float32)
    details=[]
    for ci,context in enumerate(['A','B','C']):
        cell_weights=mixtures[context]['cell_profiles']
        condition_weights=np.array([mixtures[context]['cd4_conditions']['CD4_'+c] for c in CONDITIONS])
        # Basal expression does not identify culture condition; regularize toward consensus.
        condition_weights=.5*condition_weights+.5/len(CONDITIONS)
        for ti in range(cd4.shape[1]):
            w=condition_weights*cd4_available[:,ti]
            # Missing conditions receive the available condition consensus rather than zeros.
            if w.sum()==0 and cd4_available[:,ti].any():
                w=cd4_available[:,ti].astype(float)
            c_effect=np.einsum('c,cg->g',w,cd4[:,ti])/max(w.sum(),1e-12)
            available=np.array([cd4_available[:,ti].any(),k562_available[ti],h1_available[ti]])
            weights=np.array([cell_weights['CD4'],cell_weights['K562'],cell_weights['H1']])*available
            # Unsupported HCT116/HEK293T basal components cannot supply perturbation responses.
            # A source with a vanishing basal fit still supplies the measured fallback if alone.
            if weights.sum()<1e-8:
                weights=available.astype(float)
            if weights.sum()>0:
                weights/=weights.sum()
            masks=np.stack([np.ones(cd4.shape[2],bool) if cd4_gene_mask is None else cd4_gene_mask,
                np.ones(cd4.shape[2],bool) if k562_gene_mask is None else k562_gene_mask,
                np.ones(cd4.shape[2],bool) if h1_gene_mask is None else h1_gene_mask[ti]]) & available[:,None]
            gene_weights=weights[:,None]*masks
            denominator=gene_weights.sum(axis=0)
            empty=(denominator<1e-12) & masks.any(axis=0)
            gene_weights[:,empty]=masks[:,empty]
            denominator=gene_weights.sum(axis=0)
            output[ci,ti]=(gene_weights[0]*c_effect+gene_weights[1]*k562[ti]+gene_weights[2]*h1[ti])/np.maximum(denominator,1e-12)
            details.append({'context':context,'target_index':ti,'cd4_weight':weights[0],'k562_weight':weights[1],
                'h1_weight':weights[2],'cd4_available':bool(available[0]),'k562_available':bool(available[1]),'h1_available':bool(available[2])})
    return output,pd.DataFrame(details)


def select_templates(genes,targets,manifest,analysis):
    summaries=[]
    for ci,context in enumerate(['A','B','C']):
        x,guides,libs,_,_,_=build.load_control(ROOT,context,genes,manifest)
        baseline=np.asarray(x.sum(axis=0)).ravel().astype(float)
        baseline=np.log1p(50000*baseline/baseline.sum())
        downstream=~np.isin(genes,targets)
        scores=[]; chosen=None
        for candidate in range(64):
            rng=np.random.default_rng(np.random.SeedSequence([20261005,ci,candidate]))
            rows=build.guide_stratified_sample(guides,libs,400,rng)
            counts=np.asarray(x[rows].sum(axis=0)).ravel().astype(float)
            profile=np.log1p(50000*counts/counts.sum())
            error=float(np.mean((profile[downstream]-baseline[downstream])**2))
            scores.append({'context':context,'candidate':candidate,'control_profile_mse':error,'mean_library':float(libs[rows].mean()),'median_library':float(np.median(libs[rows]))})
            if chosen is None or error<chosen[0]:
                chosen=(error,rows,candidate)
        error,rows,candidate=chosen
        pd.DataFrame({'control_row':rows,'ntc_id':guides[rows],'total_counts':libs[rows]}).to_csv(analysis/f'control_template_{context}.csv',index=False)
        summaries.extend(scores)
        print(context,'shared balanced template selected',candidate,'error',error,flush=True)
        del x; gc.collect()
    pd.DataFrame(summaries).to_csv(analysis/'control_template_candidates.csv',index=False)


def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest,genes,targets=build.read_spec(ROOT)
    selection=OUT/'cd4_model_selection.json'
    settings=json.loads(selection.read_text())['selected'] if selection.exists() else validate_cd4(OUT)
    cd4,ca,cg=fit_cd4(settings,genes,targets,OUT)
    h1,ha,hg=h1_response_model(genes,targets,OUT)
    kdir=OUT/'k562_source';kdir.mkdir(exist_ok=True)
    km=build.reference_model(ROOT,genes,targets,kdir,shrinkage=.3,noise_z=1.96)
    k562=np.zeros((len(targets),len(genes)),dtype=np.float32)
    k562[:,km[3]]=.3*km[0][0]
    mixtures=context_mixtures(genes,OUT)
    kg=np.zeros(len(genes),bool); kg[km[3]]=True
    effects,details=combine_sources(cd4,ca,k562,km[1][0],h1,ha,mixtures,cg,kg,hg)
    details['target_gene']=np.tile(targets,3)
    details.to_csv(OUT/'target_source_weights.csv',index=False)
    for ti,target in enumerate(targets):
        effects[:,ti,np.flatnonzero(genes==target)[0]]=0
    np.savez_compressed(OUT/'reference_model_arrays.npz',effects=effects,available=np.any(effects!=0,axis=2),
        basal=np.zeros((3,len(genes)),dtype=np.float32),dest_positions=np.arange(len(genes)),genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),cell_lines=np.array(['A','B','C']))
    source_summary={'model':'Context-specific uncertainty-calibrated CD4 condition ensemble with measured H1 and conservative K562 fallback',
        'context_specific':True,'basal_profile_available':False,'context_mixture_weights':mixtures,
        'condition_weight_regularization':'Half fitted CD4 condition weights plus half equal-condition consensus; fixed conservative assumption',
        'model_arrays_sha256':build.sha256(OUT/'reference_model_arrays.npz'),
        'source_data_sha256':{name:build.sha256(ROOT/'reference/cd4_probe'/name) for name in ['merged.npz','merged.csv','CE0006864_CE0008162.npz','CE0006864_CE0008162.csv','CE0008678_CE0010866.npz','CE0008678_CE0010866.csv']},
        'cd4_model_selection':json.loads(selection.read_text()),'h1_model_selection':json.loads((OUT/'h1_model_selection.json').read_text()),
        'k562_shrinkage_assumption':.3,'noise_z_assumption':1.96,'source_weights':'Normalized per target AND measured gene over available CD4/K562/H1 responses; HCT116 and HEK293T response estimates unavailable',
        'simulation':'400 common guide/depth-balanced control-template cells per context for all targets; zero-aware upward counts',
        'limitations':['No challenge perturbation ground truth and no estimate of leaderboard score.',
        'Basal mixture coefficients are expression-fit weights, not identified biological cell fractions.',
        'Basal fit is not validation of response-transfer weights.',
        'K562 shrinkage is conservative and unvalidated across contexts.',
        'CD4 donor subsets are enriched and source H1 subset is small.'],
        'script_sha256':build.sha256(Path(__file__)),'extraction_script_sha256':build.sha256(ROOT/'extract_h1_responses.py')}
    build.write_json(OUT/'reference_model.json',source_summary)
    select_templates(genes,targets,manifest,OUT)
    print('V3 source fitting and template selection complete',flush=True)


def fit_model(root,genes,targets,analysis,*unused):
    if root.resolve()!=ROOT.resolve():
        raise ValueError('Prepare model using this project root')
    with np.load(analysis/'reference_model_arrays.npz',allow_pickle=False) as z:
        if not np.array_equal(z['genes'],genes) or not np.array_equal(z['targets'],targets):
            raise ValueError('Prepared V3 axes differ')
        effects=z['effects'].copy(); available=z['available'].copy()
    summary=json.loads((analysis/'reference_model.json').read_text())
    if summary['model_arrays_sha256']!=build.sha256(analysis/'reference_model_arrays.npz'):
        raise ValueError('Prepared model array checksum mismatch')
    if summary['script_sha256']!=build.sha256(Path(__file__)):
        raise ValueError('Model fitter changed; rerun fit_context_ensemble.py')
    return effects,available,np.zeros((3,len(genes))),np.arange(len(genes)),{'shrinkage':1.,'temperature':0.},['A','B','C'],summary


if __name__=='__main__':
    prepare()
