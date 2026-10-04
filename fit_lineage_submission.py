"""Fit V4 lineage-guided responses with target-held-out shared response calibration."""
from pathlib import Path
import argparse
import json
import gc
import shutil
import numpy as np
import pandas as pd
from scipy import sparse, stats
import build_submission as build
import fit_context_ensemble as v3
from build_cd4_model import load_subset, eligible, align_disjoint

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v4_lineage_response/analysis'
CONDITIONS=v3.CONDITIONS

def bh_valid(p,valid):
    """BH within each row's finite, eligible gene universe."""
    order=np.argsort(np.where(valid,p,2.),axis=1,kind='stable')
    values=np.take_along_axis(p,order,axis=1)
    present=np.take_along_axis(valid,order,axis=1)
    number=valid.sum(axis=1)[:,None]
    adjusted=np.where(present,values*number/np.arange(1,p.shape[1]+1)[None,:],1.)
    adjusted=np.minimum.accumulate(adjusted[:,::-1],axis=1)[:,::-1]
    result=np.ones_like(p)
    np.put_along_axis(result,order,np.minimum(adjusted,1.),axis=1)
    return result


def response_prior(values,errors,valid,conditions,train):
    """Other-target condition mean; training masks also exclude each row's own target."""
    result={}
    for condition in CONDITIONS:
        use=train & (conditions==condition)
        weight=valid[use]/(errors[use]**2+.04)
        total=weight.sum(axis=0)
        result[condition]=np.clip(np.divide((weight*values[use]).sum(axis=0),total,
            out=np.zeros(values.shape[1]),where=total>0),-.5,.5)
    return result


def predict_cd4(values,errors,conditions,setting,prior):
    specific=setting['specific_gain']*values*.04/(.04+errors**2)
    background=np.stack([prior[c] for c in conditions])
    return np.clip(specific+setting['prior_gain']*background,-3,3)


def evaluate(left_prediction,right_prediction,left,right,lq,rq,valid,left_se400,right_se400,rows):
    indices=np.flatnonzero(rows);v=valid[indices]
    pl=left_prediction[indices];pr=right_prediction[indices];l=left[indices];r=right[indices]
    denom=np.maximum(v.sum(axis=1),1)
    loss=(((pl-r)**2+(pr-l)**2)*v).sum(axis=1)/denom
    zero=((l*l+r*r)*v).sum(axis=1)/denom
    summaries=[]
    for prediction,truth,truthq,se in [(pl,r,rq[indices],left_se400[indices]),(pr,l,lq[indices],right_se400[indices])]:
        p=2*stats.norm.sf(np.abs(prediction)/np.maximum(se,1e-12))
        called=(bh_valid(p,v)<.05)&v&(truth!=0)
        real=(truthq<.05)&v
        correct=called&(np.sign(prediction)==np.sign(truth))
        npred=called.sum(axis=1);nreal=real.sum(axis=1)
        fidelity=np.divide(correct.sum(axis=1),np.maximum(npred,nreal),out=np.full(len(indices),np.nan),where=np.maximum(npred,nreal)>0)
        union=(called|real).sum(axis=1)
        jac=np.divide((called&real).sum(axis=1),union,out=np.ones(len(indices)),where=union>0)
        trusted=called&real
        summaries.append({'direction_fidelity_proxy':float(np.nanmean(fidelity)),
           'jaccard_proxy':float(jac.mean()),'mean_predicted_calls':float(npred.mean()),
           'mean_reference_calls':float(nreal.mean()),'sign_accuracy_on_reference_DE':float(correct[trusted].mean()) if trusted.any() else None})
    result={'contrasts':len(indices),'loss':float(loss.mean()),'zero_loss':float(zero.mean()),
         'relative_zero_loss':float(loss.sum()/zero.sum())}
    for key in ['direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']:
        result[key]=float(np.mean([x[key] for x in summaries]))
    accuracy=[x['sign_accuracy_on_reference_DE'] for x in summaries if x['sign_accuracy_on_reference_DE'] is not None]
    result['sign_accuracy_on_reference_DE']=float(np.mean(accuracy)) if accuracy else None
    result['selection_score']=result['direction_fidelity_proxy']+result['jaccard_proxy']+.5*(1-result['relative_zero_loss'])
    return result


def choose(records):
    acceptable=[x for x in records if x['relative_zero_loss']<=1.02]
    if not acceptable:raise ValueError('No source-validation candidate passes the expression-error guardrail')
    return max(acceptable,key=lambda x:(x['selection_score'],-x['setting']['prior_gain'],-x['setting'].get('specific_gain',x['setting'].get('cd4_coefficient',0)+x['setting'].get('k562_coefficient',0))))


def donor_inputs():
    directory=ROOT/'reference/cd4_probe'
    provenance=json.loads((ROOT/'analysis/model_audit/cd4_extraction.json').read_text())
    for record in provenance:
        if build.sha256(directory/record['subset'])!=record['sha256']:raise ValueError('CD4 source changed')
    names=[Path(x['subset']).stem for x in provenance if Path(x['subset']).stem!='merged']
    if len(names)!=2 or set(names[0].split('_'))&set(names[1].split('_')):raise ValueError('Donor groups overlap')
    lm,ld=load_subset(directory,names[0]);rm,rd=load_subset(directory,names[1])
    rows,left,right,lq,rq,valid=align_disjoint((lm,ld),(rm,rd))
    symbols=np.array(sorted(set(ld['genes'])&set(rd['genes'])))
    le=ld['lfcSE'][rows.position_left.to_numpy()][:,[list(ld['genes']).index(g) for g in symbols]].astype(float)
    re=rd['lfcSE'][rows.position_right.to_numpy()][:,[list(rd['genes']).index(g) for g in symbols]].astype(float)
    valid &= np.isfinite(le)&np.isfinite(re)&(le>0)&(re>0)
    le=np.nan_to_num(le,nan=1e6,posinf=1e6,neginf=1e6);re=np.nan_to_num(re,nan=1e6,posinf=1e6,neginf=1e6)
    conditions=rows.culture_condition.to_numpy();targets=rows.target_contrast_gene_name_left.to_numpy()
    se400left=le*np.sqrt(rows.n_cells_target_left.to_numpy()[:,None]/400)
    se400right=re*np.sqrt(rows.n_cells_target_right.to_numpy()[:,None]/400)
    return rows,left,right,lq,rq,valid,le,re,conditions,targets,se400left,se400right,names


def select_cd4():
    OUT.mkdir(parents=True,exist_ok=True)
    rows,left,right,lq,rq,valid,le,re,conditions,targets,se400l,se400r,names=donor_inputs()
    folds=v3.target_folds(targets)
    grid=[{'specific_gain':a,'prior_gain':b} for a in [.5,1.,1.5] for b in [0.,.5,1.,2.]]
    nested=[];outer=[];inner_records=[]
    for fold in range(5):
        test=folds==fold;train=~test;inner=v3.target_folds(targets[train],3);indices=np.flatnonzero(train)
        candidate_rows=[[] for _ in grid]
        for part in range(3):
            tune=np.zeros(len(rows),bool);tune[indices[inner==part]]=True;fit=train&~tune
            lp=response_prior(left,le,valid,conditions,fit);rp=response_prior(right,re,valid,conditions,fit)
            for ci,setting in enumerate(grid):
                pl=predict_cd4(left,le,conditions,setting,lp);pr=predict_cd4(right,re,conditions,setting,rp)
                result=evaluate(pl,pr,left,right,lq,rq,valid,se400l,se400r,tune)
                candidate_rows[ci].append(result)
        choices=[]
        for ci,setting in enumerate(grid):
            values=candidate_rows[ci];n=sum(x['contrasts'] for x in values)
            avg={key:sum(x[key]*x['contrasts'] for x in values)/n for key in ['loss','zero_loss','direction_fidelity_proxy','jaccard_proxy','selection_score']}
            avg['relative_zero_loss']=avg['loss']/avg['zero_loss'];avg['setting']=setting;avg['candidate']=ci
            choices.append(avg);inner_records.append({'outer_fold':fold,**avg})
        chosen=choose(choices)
        lp=response_prior(left,le,valid,conditions,train);rp=response_prior(right,re,valid,conditions,train)
        for ci,setting in enumerate(grid):
            pl=predict_cd4(left,le,conditions,setting,lp);pr=predict_cd4(right,re,conditions,setting,rp)
            result=evaluate(pl,pr,left,right,lq,rq,valid,se400l,se400r,test)
            record={'fold':fold,'candidate':ci,'setting':setting,**result};outer.append(record)
            if ci==chosen['candidate']:
                nested.append({**record,'train_targets':sorted(set(targets[train])),'test_targets':sorted(set(targets[test]))})
        print(f'CD4 outer fold {fold}: {chosen["setting"]}; held-target proxy score {nested[-1]["selection_score"]:.5f}',flush=True)
    combined=[]
    for ci,setting in enumerate(grid):
        group=[r for r in outer if r['candidate']==ci];n=sum(x['contrasts'] for x in group)
        result={key:sum(x[key]*x['contrasts'] for x in group)/n for key in ['loss','zero_loss','direction_fidelity_proxy','jaccard_proxy','selection_score','mean_predicted_calls','mean_reference_calls']}
        result['relative_zero_loss']=result['loss']/result['zero_loss'];result.update({'candidate':ci,'setting':setting});combined.append(result)
    selected=choose(combined)
    n=sum(x['contrasts'] for x in nested)
    nested_summary={key:sum(x[key]*x['contrasts'] for x in nested)/n for key in ['loss','zero_loss','direction_fidelity_proxy','jaccard_proxy','selection_score','mean_predicted_calls','mean_reference_calls']}
    nested_summary['relative_zero_loss']=nested_summary['loss']/nested_summary['zero_loss']
    baseline=combined[0]
    # A source proxy is insufficient to force a change if held-out performance does not improve.
    use=nested_summary['selection_score']>baseline['selection_score'] and nested_summary['relative_zero_loss']<=1.02
    final=selected['setting'] if use else grid[0]
    pd.DataFrame([{**r,'setting':json.dumps(r['setting'],sort_keys=True)} for r in outer]).to_csv(OUT/'cd4_outer_candidate_scores.csv',index=False)
    pd.DataFrame([{**r,'setting':json.dumps(r['setting'],sort_keys=True)} for r in combined]).to_csv(OUT/'cd4_candidate_summary.csv',index=False)
    build.write_json(OUT/'cd4_lineage_selection.json',{'selected':final,'selected_all_fold_candidate':selected,'nested_procedure':nested_summary,
        'v3_fixed_baseline':baseline,'changed_estimator_supported':bool(use),'nested_folds':nested,'inner_selections':inner_records,
        'disjoint_donor_groups':names,'targets':len(set(targets)),'contrasts':len(rows),'genes':left.shape[1],
        'prediction_call_proxy':'Two-sided normal p-values from predicted log2FC and source SE scaled to 400 cells; within-row eligible-gene BH; not a Wilcoxon or competition score',
        'selection_objective':'Direction-fidelity proxy + Jaccard proxy + 0.5*(1-relative bidirectional log2FC MSE), subject to MSE <=1.02 of zero',
        'scope':'Nested target-held-out within-CD4 donor validation, not anonymous-context or leaderboard validation; enriched donor subsets',
        'source_hashes':{x:build.sha256(ROOT/'reference/cd4_probe'/x) for x in ['merged.npz','merged.csv',names[0]+'.npz',names[0]+'.csv',names[1]+'.npz',names[1]+'.csv']}})
    print('CD4 selected',final,'nested',nested_summary,'V3 baseline',baseline,flush=True)
    return final


def h1_inputs():
    directory=ROOT/'reference/h1_response';record=json.loads((directory/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:
        if build.sha256(directory/name)!=record[key]:raise ValueError('H1 source changed')
    matrix=sparse.load_npz(directory/'selected_counts.npz').astype(float);meta=pd.read_csv(directory/'selected_cells.csv')
    source=pd.read_csv(directory/'gene_names.csv').gene_name.to_numpy();libraries=np.asarray(matrix.sum(axis=1)).ravel()
    cpm=matrix.multiply((1e6/libraries)[:,None]).tocsr();batches=sorted(meta.batch.unique())
    partition=meta.batch.map({batch:i%2 for i,batch in enumerate(batches)}).to_numpy()
    controls=meta.target_gene.to_numpy()=='non-targeting';batch_labels=meta.batch.to_numpy()
    baselines={}
    for batch in batches:
        x=cpm[controls&(batch_labels==batch)];mean=np.asarray(x.mean(axis=0)).ravel()
        variance=np.maximum(np.asarray(x.power(2).mean(axis=0)).ravel()-mean**2,0)
        baselines[batch]=(mean,variance/max(x.shape[0]-1,1))
    targets=sorted(set(meta.target_gene)-{'non-targeting'});all_values=[]
    for part in [0,1,None]:
        values=[];errors=[];bases=[];qvalues=[];nvalues=[]
        for target in targets:
            selected=(meta.target_gene.to_numpy()==target)
            if part is not None:selected &= partition==part
            idx=np.flatnonzero(selected);n=len(idx)
            if n<20:raise ValueError('H1 held-batch response too sparse')
            x=cpm[idx];mean=np.asarray(x.mean(axis=0)).ravel()
            varmean=np.maximum(np.asarray(x.power(2).mean(axis=0)).ravel()-mean**2,0)/max(n-1,1)
            counts=meta.iloc[idx].batch.value_counts();base=np.zeros(len(source));basevar=np.zeros(len(source))
            for batch,count in counts.items():
                bm,bv=baselines[batch];base+=count/n*bm;basevar+=(count/n)**2*bv
            lfc=np.clip(np.log2((mean+1)/(base+1)),-3,3)
            se=np.sqrt(varmean/(mean+1)**2+basevar/(base+1)**2)/np.log(2)
            valid=(base>5)&(source!=target)&(se>0)
            # Actual held-batch Wilcoxon calls, not the estimator's normal approximation.
            control_idx=np.flatnonzero(controls if part is None else controls&(partition==part))
            q=np.ones(len(source));columns=np.flatnonzero(valid)
            for start in range(0,len(columns),256):
                cols=columns[start:start+256]
                p=stats.mannwhitneyu(cpm[idx][:,cols].toarray(),cpm[control_idx][:,cols].toarray(),axis=0,method='asymptotic').pvalue
                q[cols]=p
            q[valid]=stats.false_discovery_control(q[valid],method='bh')
            values.append(lfc);errors.append(se);bases.append(base);qvalues.append(q);nvalues.append(n)
        all_values.append({'values':np.array(values),'errors':np.array(errors),'baselines':np.array(bases),
                           'q':np.array(qvalues),'n':np.array(nvalues)})
    return np.array(targets),source,all_values,record


def predict_h1(values,errors,conditions,setting,prior):
    specific=setting['specific_gain']*np.sign(values)*np.maximum(np.abs(values)-errors,0)
    background=np.stack([prior[c] for c in conditions])
    return np.clip(specific+setting['prior_gain']*background,-3,3)


def select_h1(targets,genes,parts):
    left,right=parts[:2];l=left['values'];r=right['values'];le=left['errors'];re=right['errors']
    valid=(left['baselines']>5)&(right['baselines']>5)&(genes[None,:]!=targets[:,None])&(le>0)&(re>0)
    conditions=np.repeat('Rest',len(targets));sl=le*np.sqrt(left['n'][:,None]/400);sr=re*np.sqrt(right['n'][:,None]/400)
    grid=[{'specific_gain':a,'prior_gain':b} for a in [.7,1.,1.5] for b in [0.,.5,1.,2.]]
    folds=v3.target_folds(targets,3);outer=[];nested=[]
    for fold in range(3):
        test=folds==fold;train=~test;inner=v3.target_folds(targets[train],2);indices=np.flatnonzero(train)
        scores=[[] for _ in grid]
        for part in range(2):
            tune=np.zeros(len(targets),bool);tune[indices[inner==part]]=True;fit=train&~tune
            lp=response_prior(l,le,valid,conditions,fit);rp=response_prior(r,re,valid,conditions,fit)
            for ci,setting in enumerate(grid):
                pl=predict_h1(l,le,conditions,setting,lp);pr=predict_h1(r,re,conditions,setting,rp)
                scores[ci].append(evaluate(pl,pr,l,r,left['q'],right['q'],valid,sl,sr,tune))
        choices=[]
        for ci,setting in enumerate(grid):
            total=sum(x['contrasts'] for x in scores[ci]);record={key:sum(x[key]*x['contrasts'] for x in scores[ci])/total for key in ['loss','zero_loss','selection_score']}
            record.update({'relative_zero_loss':record['loss']/record['zero_loss'],'setting':setting,'candidate':ci});choices.append(record)
        chosen=choose(choices);lp=response_prior(l,le,valid,conditions,train);rp=response_prior(r,re,valid,conditions,train)
        for ci,setting in enumerate(grid):
            pl=predict_h1(l,le,conditions,setting,lp);pr=predict_h1(r,re,conditions,setting,rp)
            result=evaluate(pl,pr,l,r,left['q'],right['q'],valid,sl,sr,test)
            record={'fold':fold,'candidate':ci,'setting':setting,**result};outer.append(record)
            if ci==chosen['candidate']:nested.append({**record,'train_targets':targets[train].tolist(),'test_targets':targets[test].tolist()})
        print('H1 outer fold',fold,chosen['setting'],'proxy score',nested[-1]['selection_score'],flush=True)
    summary=[]
    for ci,setting in enumerate(grid):
        group=[x for x in outer if x['candidate']==ci];n=sum(x['contrasts'] for x in group)
        result={key:sum(x[key]*x['contrasts'] for x in group)/n for key in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']}
        result.update({'setting':setting,'candidate':ci,'relative_zero_loss':result['loss']/result['zero_loss']});summary.append(result)
    selected=choose(summary);n=sum(x['contrasts'] for x in nested)
    procedure={key:sum(x[key]*x['contrasts'] for x in nested)/n for key in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']}
    procedure['relative_zero_loss']=procedure['loss']/procedure['zero_loss'];baseline=summary[0]
    supported=procedure['selection_score']>baseline['selection_score'] and procedure['relative_zero_loss']<=1.02
    final=selected['setting'] if supported else grid[0]
    pd.DataFrame([{**x,'setting':json.dumps(x['setting'],sort_keys=True)} for x in outer]).to_csv(OUT/'h1_outer_candidate_scores.csv',index=False)
    pd.DataFrame([{**x,'setting':json.dumps(x['setting'],sort_keys=True)} for x in summary]).to_csv(OUT/'h1_candidate_summary.csv',index=False)
    build.write_json(OUT/'h1_lineage_selection.json',{'selected':final,'nested_procedure':procedure,'v3_fixed_baseline':baseline,
        'changed_estimator_supported':bool(supported),'nested_folds':nested,'targets':len(targets),
        'scope':'Target-held-out, disjoint public H1 batches; 400-cell normal call proxy compared with actual source Wilcoxon calls; not challenge accuracy'})
    print('H1 selected',final,'nested',procedure,'V3 baseline',baseline,flush=True)
    return final


def fit_cd4_arrays(setting,genes,targets,conditions_weights):
    meta,data=load_subset(ROOT/'reference/cd4_probe','merged');source=data['genes'];good=eligible(meta)
    condition=meta.culture_condition.to_numpy();row_target=meta.target_contrast_gene_name.to_numpy()
    raw=np.clip(np.nan_to_num(data['log_fc'].astype(float),nan=0,posinf=0,neginf=0),-3,3)
    se=np.nan_to_num(data['lfcSE'].astype(float),nan=1e6,posinf=1e6,neginf=1e6)
    valid=np.isfinite(data['log_fc'])&np.isfinite(data['lfcSE'])&(data['lfcSE']>0)&(source[None,:]!=row_target[:,None])
    priors=response_prior(raw,se,valid,condition,good)
    position={g:i for i,g in enumerate(genes)};si=np.array([i for i,g in enumerate(source) if g in position]);di=np.array([position[source[i]] for i in si])
    specific=np.zeros((3,len(targets),len(genes)),np.float32);supported=np.zeros_like(specific,dtype=bool);target_available=np.zeros((3,len(targets)),bool)
    # Do not put a target's own measurements into its fitted shared component.
    backgrounds=np.empty_like(specific)
    for ti,target in enumerate(targets):
        own=row_target==target;prior=response_prior(raw,se,valid,condition,good&~own)
        for ci,name in enumerate(CONDITIONS):
            backgrounds[ci,ti]=0;backgrounds[ci,ti,di]=setting['prior_gain']*prior[name][si]
            idx=np.flatnonzero(good&own&(condition==name))
            if len(idx):
                weight=np.sqrt(meta.n_cells_target.to_numpy()[idx]);weight=weight[:,None]*valid[idx]
                total=weight.sum(axis=0)
                fitted=setting['specific_gain']*raw[idx]*.04/(.04+se[idx]**2)
                mean=np.divide((weight*fitted).sum(axis=0),total,out=np.zeros(len(source)),where=total>0)
                specific[ci,ti,di]=mean[si];supported[ci,ti,di]=total[si]>0;target_available[ci,ti]=True
    effects=[];target_effects=[];masks=[]
    for context in ['A','B','C']:
        weights=np.array([conditions_weights[context]['cd4_conditions']['CD4_'+c] for c in CONDITIONS]);weights=.5*weights+.5/3
        actual=weights[:,None,None]*supported;den=actual.sum(axis=0)
        target=np.divide((actual*specific).sum(axis=0),den,out=np.zeros_like(specific[0]),where=den>0)
        background=np.einsum('c,ctg->tg',weights,backgrounds)
        effects.append(np.clip(target+background,-3,3));target_effects.append(target);masks.append(supported.any(axis=0))
    return np.array(effects,np.float32),np.array(target_effects,np.float32),np.array(masks),target_available.any(axis=0)


def select_cross_lineage(h1_targets,h1_genes,parts,cd4_effect,k562_effect,targets,genes,h1_setting):
    positions={g:i for i,g in enumerate(genes)};si=np.array([i for i,g in enumerate(h1_genes) if g in positions]);di=np.array([positions[h1_genes[i]] for i in si])
    ti=np.array([list(targets).index(t) for t in h1_targets]);conditions=np.repeat('Rest',len(ti))
    l,r=parts[:2];left=l['values'][:,si];right=r['values'][:,si];le=l['errors'][:,si];re=r['errors'][:,si]
    valid=(l['baselines'][:,si]>5)&(r['baselines'][:,si]>5)&(h1_genes[si][None,:]!=h1_targets[:,None])&(le>0)&(re>0)
    sl=le*np.sqrt(l['n'][:,None]/400);sr=re*np.sqrt(r['n'][:,None]/400)
    # These fixed predictions use only CD4/K562 target data, never an H1 target response.
    c=cd4_effect[ti][:,di];k=k562_effect[ti][:,di]
    grid=[{'cd4_coefficient':a,'k562_coefficient':b,'prior_gain':p} for a in [0.,.5,1.] for b in [0.,.5,1.] for p in [0.,.5,1.]]
    records=[];nested=[]
    for held in range(len(ti)):
        train=np.arange(len(ti))!=held;test=~train
        lp=response_prior(left,le,valid,conditions,train);rp=response_prior(right,re,valid,conditions,train)
        candidates=[]
        for ci,setting in enumerate(grid):
            # Same-transfer prediction evaluated against both independent H1 batch partitions.
            base=setting['cd4_coefficient']*c+setting['k562_coefficient']*k
            pl=np.clip(base+setting['prior_gain']*np.stack([lp[x] for x in conditions]),-3,3)
            pr=np.clip(base+setting['prior_gain']*np.stack([rp[x] for x in conditions]),-3,3)
            fit=evaluate(pl,pr,left,right,l['q'][:,si],r['q'][:,si],valid,sl,sr,train)
            fit.update({'setting':setting,'candidate':ci});candidates.append(fit)
            result=evaluate(pl,pr,left,right,l['q'][:,si],r['q'][:,si],valid,sl,sr,test)
            records.append({'held_target':h1_targets[held],'candidate':ci,'setting':setting,**result})
        chosen=choose(candidates)
        nested.append({**next(x for x in records if x['held_target']==h1_targets[held] and x['candidate']==chosen['candidate']),
            'train_targets':h1_targets[train].tolist(),'test_targets':[h1_targets[held]]})
    summary=[]
    for ci,setting in enumerate(grid):
        group=[x for x in records if x['candidate']==ci]
        avg={key:float(np.mean([x[key] for x in group])) for key in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']}
        avg.update({'relative_zero_loss':avg['loss']/avg['zero_loss'],'candidate':ci,'setting':setting});summary.append(avg)
    selected=choose(summary);procedure={key:float(np.mean([x[key] for x in nested])) for key in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']}
    procedure['relative_zero_loss']=procedure['loss']/procedure['zero_loss']
    pd.DataFrame([{**x,'setting':json.dumps(x['setting'],sort_keys=True)} for x in records]).to_csv(OUT/'cross_lineage_held_target_scores.csv',index=False)
    pd.DataFrame([{**x,'setting':json.dumps(x['setting'],sort_keys=True)} for x in summary]).to_csv(OUT/'cross_lineage_candidate_summary.csv',index=False)
    # The 2% MSE guardrail constrains training/final-candidate selection, not a guaranteed
    # held-out error bound. Accept only an improved held-target composite; report the
    # nested squared-error tradeoff separately rather than claiming MSE improvement.
    final=selected['setting'] if procedure['selection_score']>summary[0]['selection_score'] else grid[0]
    build.write_json(OUT/'cross_lineage_selection.json',{'selected':final,'nested_procedure':procedure,'zero_baseline':summary[0],
       'selected_all_fold_candidate':selected,'nested_folds':nested,'targets':len(ti),
       'acceptance_policy':'Candidate selection obeys 1.02 zero-relative MSE guardrail; held-target procedure must improve composite score, with its MSE tradeoff reported rather than claimed as an improvement',
       'scope':'Source-context transfer into 13 public H1 targets; H1 response for each scored target excluded from all predictions; small, selected target panel and not B/C validation'})
    print('Cross-lineage selected',final,'nested',procedure,flush=True)
    return final

def h1_final_arrays(h1_targets,source,parts,genes,targets,setting,transfer_setting):
    full=parts[2];position={g:i for i,g in enumerate(genes)}
    si=np.array([i for i,g in enumerate(source) if g in position]);di=np.array([position[source[i]] for i in si])
    condition=np.repeat('Rest',len(h1_targets));valid=(full['baselines']>5)&(full['errors']>0)&(source[None,:]!=h1_targets[:,None])
    prior_all=response_prior(full['values'],full['errors'],valid,condition,np.ones(len(h1_targets),bool))['Rest']
    background=np.zeros((len(targets),len(genes)),np.float32)
    specific=np.zeros_like(background);direct_mask=np.zeros_like(background,dtype=bool);available=np.zeros(len(targets),bool)
    for ti,target in enumerate(targets):
        use=h1_targets!=target
        prior=response_prior(full['values'],full['errors'],valid,condition,use)['Rest'] if not use.all() else prior_all
        background[ti,di]=transfer_setting['prior_gain']*prior[si]
        idx=np.flatnonzero(h1_targets==target)
        if len(idx):
            j=idx[0];own=np.flatnonzero(source==target)
            if full['n'][j]>=40 and len(own) and full['values'][j,own[0]]<=np.log2(.7):
                specific[ti,di]=setting['specific_gain']*np.sign(full['values'][j,si])*np.maximum(np.abs(full['values'][j,si])-full['errors'][j,si],0)
                direct_mask[ti,di]=valid[j,si];available[ti]=True
                # Direct-response fits use their own independently selected background coefficient.
                background[ti,di]=setting['prior_gain']*prior[si]
    measured=np.zeros(len(genes),bool);measured[di]=np.any(valid[:,si],axis=0)
    background[:,~measured]=0;specific[~direct_mask]=0
    return specific,background,direct_mask,available


def assemble_effects(cd4,cd4_specific,cd4_mask,k562,k_mask,h1,h1_prior,h1_mask,h1_available,transfer,means,genes,targets):
    """Matched CD4 response for A; independently calibrated nonimmune transfer for B/C."""
    effects=np.zeros_like(cd4,dtype=np.float32);records=[]
    for ci,context in enumerate(['A','B','C']):
        if context=='A':
            # Prefer measured, matched-lineage CD4; K562 only fills genuinely unsupported coordinates.
            effect=np.where(cd4_mask[ci],cd4[ci],np.where(k_mask,k562,0))
        else:
            transferred=transfer['cd4_coefficient']*cd4_specific[ci]+transfer['k562_coefficient']*k562
            effect=np.where(h1_mask,h1,transferred)+h1_prior
        expression_gate=np.clip(means[ci]/5,0,1)
        target_gate=np.clip(means[ci,[list(genes).index(t) for t in targets]]/5,0,1)
        effect=np.clip(effect*expression_gate[None,:]*target_gate[:,None],-3,3)
        for ti,target in enumerate(targets):
            effect[ti,list(genes).index(target)]=0
            records.append({'context':context,'target_gene':target,'target_index':ti,'target_expression_gate':float(target_gate[ti]),
               'cd4_measured_genes':int(cd4_mask[ci,ti].sum()),'k562_measured_genes':int(k_mask[ti].sum()),
               'h1_direct_response_available':bool(h1_available[ti]),'nonzero_effect_genes':int(np.count_nonzero(effect[ti])),
               'effects_above_lfc_0_1':int((np.abs(effect[ti])>.1).sum()),'shared_h1_prior_nonzero_genes':0 if ci==0 else int(np.count_nonzero(h1_prior[ti]))})
        effects[ci]=effect
    return effects,pd.DataFrame(records)


def prepare():
    OUT.mkdir(parents=True,exist_ok=True);manifest,genes,targets=build.read_spec(ROOT)
    # Frozen source checksums and original-matrix marker checks are prerequisites, not identity labels.
    audit=json.loads((ROOT/'analysis/target_archive_audit/archive_validation.json').read_text())
    if audit['status']!='passed' or audit['mismatched_values']!=0:raise ValueError('Lineage input audit did not pass')
    for record in audit['original_inputs']:
        if build.sha256(ROOT/record['file'])!=record['sha256']:raise ValueError('Challenge controls changed')
    selection=OUT/'cd4_lineage_selection.json'
    csetting=json.loads(selection.read_text())['selected'] if selection.exists() else select_cd4()
    mixtures=json.loads((ROOT/'experiments/v3_context_ensemble/analysis/context_mixture_weights.json').read_text())
    cd4,cs,cm,ca=fit_cd4_arrays(csetting,genes,targets,mixtures)
    old=OUT/'k562_source';old.mkdir(exist_ok=True)
    expected=json.loads((ROOT/'experiments/v3_context_ensemble/analysis/k562_source/reference_model.json').read_text())
    if build.sha256(ROOT/'reference/K562_gwps_raw_bulk_01.h5ad')!=expected['source_file_sha256']:raise ValueError('K562 input changed')
    kfit=build.reference_model(ROOT,genes,targets,old,shrinkage=.3,noise_z=1.96)
    ke=np.zeros((len(targets),len(genes)),np.float32);km=np.zeros_like(ke,dtype=bool)
    ke[:,kfit[3]]=.3*kfit[0][0];km[:,kfit[3]]=kfit[1][0,:,None]
    del kfit;gc.collect()
    print('Preparing H1 disjoint-batch contrasts and actual source rank tests',flush=True)
    ht,hg,hparts,hrecord=h1_inputs();hsetting=select_h1(ht,hg,hparts)
    transfer=select_cross_lineage(ht,hg,hparts,cs[2],ke,targets,genes,hsetting)
    he,hprior,hm,ha=h1_final_arrays(ht,hg,hparts,genes,targets,hsetting,transfer)
    means=[];baselines=[]
    for context in ['A','B','C']:
        data=pd.read_csv(ROOT/f'analysis/gene_statistics_{context}.csv')
        if not np.array_equal(data.gene_name,genes):raise ValueError('Context mean axis differs')
        means.append(data.mean_cpm.to_numpy());baselines.append(build.sha256(ROOT/f'analysis/gene_statistics_{context}.csv'))
    means=np.array(means)
    effects,records=assemble_effects(cd4,cs,cm,ke,km,he,hprior,hm,ha,transfer,means,genes,targets)
    records.to_csv(OUT/'target_response_sources.csv',index=False)
    np.savez_compressed(OUT/'source_components.npz',cd4=cd4,cd4_specific=cs,cd4_mask=cm,cd4_available=ca,
        k562=ke,k562_mask=km,h1=he,h1_prior=hprior,h1_mask=hm,h1_available=ha,context_mean_cpm=means,
        genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str))
    np.savez_compressed(OUT/'reference_model_arrays.npz',effects=effects,available=np.any(effects!=0,axis=2),
        basal=np.zeros((3,len(genes)),np.float32),dest_positions=np.arange(len(genes)),genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),cell_lines=np.array(['A','B','C']))
    source_checks={str(path.relative_to(ROOT)):build.sha256(path) for path in [
        ROOT/'reference/cd4_probe/merged.npz',ROOT/'reference/cd4_probe/merged.csv',
        ROOT/'reference/cd4_probe/CE0006864_CE0008162.npz',ROOT/'reference/cd4_probe/CE0006864_CE0008162.csv',
        ROOT/'reference/cd4_probe/CE0008678_CE0010866.npz',ROOT/'reference/cd4_probe/CE0008678_CE0010866.csv',
        ROOT/'reference/K562_gwps_raw_bulk_01.h5ad',old/'reference_model_arrays.npz',
        ROOT/'reference/h1_response/selected_counts.npz',ROOT/'reference/h1_response/selected_cells.csv',ROOT/'reference/h1_response/gene_names.csv',
        ROOT/'experiments/v3_context_ensemble/analysis/context_mixture_weights.json',ROOT/'analysis/target_archive_audit/archive_validation.json',
        ROOT/'analysis/target_archive_audit/broader_lineage_expression_summary.csv']}
    summary={'model':'Lineage-guided, response-coverage-calibrated perturbation model v4','context_specific':True,'basal_profile_available':False,
        'context_mixture_weights':{c:{'policy':'CD4 primary with unsupported-gene K562 fallback' if c=='A' else 'H1 direct estimates when measured; H1 shared prior plus held-target-calibrated CD4/K562 transfer otherwise',
            'cd4_conditions':mixtures[c]['cd4_conditions']} for c in ['A','B','C']},
        'cd4_selection':json.loads(selection.read_text()),'h1_selection':json.loads((OUT/'h1_lineage_selection.json').read_text()),
        'cross_lineage_selection':json.loads((OUT/'cross_lineage_selection.json').read_text()),
        'cd4_settings':csetting,'h1_settings':hsetting,'cross_lineage_settings':transfer,
        'model_arrays_sha256':build.sha256(OUT/'reference_model_arrays.npz'),'source_components_sha256':build.sha256(OUT/'source_components.npz'),
        'source_data_sha256':source_checks,'context_original_inputs':audit['original_inputs'],'context_statistics_sha256':dict(zip(['A','B','C'],baselines)),
        'script_sha256':build.sha256(Path(__file__)),'reuse_fitter_sha256':build.sha256(ROOT/'fit_context_ensemble.py'),
        'expression_gates':'Multiply by min(1, context gene mean CPM/5) and min(1, context target mean CPM/5); fixed modeling assumption, not fitted biology',
        'simulation':'400 separately guide/depth-stratified original control cells per target/context; binomial thinning and zero-aware depth-scaled Poisson increases',
        'source_limitations':['No hidden competition responses or leaderboard accuracy estimate.',
           'Lineage assignments are marker-based hypotheses, not known cell identities.',
           'Shared H1 responses are based on 13 selected public targets, not an unbiased challenge response distribution.',
           'Normal 400-cell call proxies are not challenge Wilcoxon scoring; nested source validation cannot guarantee B/C transfer.',
           'Selection trades modest squared-error degradation for downstream response coverage and direction yield.',
           'No new HCT116/HEK293T response model is fit from underpowered pilot cells.'],
        'effect_counts':{c:int((np.abs(effects[i])>.1).sum()) for i,c in enumerate(['A','B','C'])}}
    build.write_json(OUT/'reference_model.json',summary)
    print('V4 prepared',json.dumps(summary['effect_counts']),flush=True)


def fit_model(root,genes,targets,analysis):
    summary=json.loads((analysis/'reference_model.json').read_text())
    if summary['script_sha256']!=build.sha256(Path(__file__)):raise ValueError('V4 fitter changed; rerun preparation')
    if summary['model_arrays_sha256']!=build.sha256(analysis/'reference_model_arrays.npz'):raise ValueError('V4 effect checksum mismatch')
    with np.load(analysis/'reference_model_arrays.npz',allow_pickle=False) as z:
        if not np.array_equal(z['genes'],genes) or not np.array_equal(z['targets'],targets):raise ValueError('V4 axes differ')
        effects=z['effects'].copy();available=z['available'].copy()
    return effects,available,np.zeros((3,len(genes))),np.arange(len(genes)),{'shrinkage':1.,'temperature':0.},['A','B','C'],summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--select-only',action='store_true')
    args=parser.parse_args()
    if args.select_only:select_cd4()
    else:prepare()
