"""Calibrate a conservative promoter-neighbor prior on disjoint CD4 donors and held targets."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import build_submission as build
import fit_context_ensemble as v3
from build_cd4_model import load_subset
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v5_atlas_calibrated/analysis'
DONORS=['CE0006864_CE0008162','CE0008678_CE0010866']


def quality(meta):
    return ((meta.n_cells_target>=100)&meta.ontarget_significant.astype(bool)&~meta.distal_offtarget_flag.astype(bool)&
            ~meta.single_guide_estimate.astype(bool)&~meta.low_target_gex.astype(bool)).to_numpy()


def ramp(distance):
    return 1-np.clip(np.log(np.maximum(distance,500)/500)/np.log(10),0,1)


def estimate(values,weights,train):
    valid=train&np.isfinite(values)&np.isfinite(weights)&(weights>0)
    return float(np.clip(np.sum(weights[valid]*values[valid])/np.sum(weights[valid]),-3,0)) if valid.any() else 0.


def predictions(left,right,scale,train,alpha,base_left=None,base_right=None,support_left=None,support_right=None):
    # Each side fits its coefficient using only its own training-target donor observations.
    mu_l=estimate(left/np.maximum(scale,.05),scale**2,train)
    mu_r=estimate(right/np.maximum(scale,.05),scale**2,train)
    pl=np.clip(alpha*mu_l*scale,-3,0);pr=np.clip(alpha*mu_r*scale,-3,0)
    if base_left is not None:
        pl=base_left.copy() if alpha==0 else np.minimum(base_left,pl)
        pr=base_right.copy() if alpha==0 else np.minimum(base_right,pr)
    if support_left is not None:
        pl=np.where(support_left,pl,base_left);pr=np.where(support_right,pr,base_right)
    return pl,pr


def loss(pl,pr,left,right,mask):
    return float(np.mean(((pl[mask]-right[mask])**2+(pr[mask]-left[mask])**2)/2))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    pairs=pd.read_csv(ROOT/'analysis/upstream_review/promoter_neighbor_pairs.csv')
    lm,ld=load_subset(ROOT/'reference/cd4_probe',DONORS[0]);rm,rd=load_subset(ROOT/'reference/cd4_probe',DONORS[1])
    lm=lm.assign(position=np.arange(len(lm)));rm=rm.assign(position=np.arange(len(rm)))
    joined=lm.loc[quality(lm)].merge(rm.loc[quality(rm)],on=['target_contrast','culture_condition'],suffixes=('_l','_r'))
    joined=joined.merge(pairs,left_on='target_contrast_gene_name_l',right_on='target')
    official=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy();target_names=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    with np.load(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz') as z:k=z['k562'].copy()
    gi={g:i for i,g in enumerate(official)};ti={t:i for i,t in enumerate(target_names)}
    lpos={g:i for i,g in enumerate(ld['genes'])};rpos={g:i for i,g in enumerate(rd['genes'])};rows=[]
    for row in joined.itertuples():
        if row.neighbor not in lpos or row.neighbor not in rpos:continue
        a=ld['log_fc'][row.position_l,lpos[row.neighbor]];b=rd['log_fc'][row.position_r,rpos[row.neighbor]]
        if not np.isfinite(a) or not np.isfinite(b):continue
        le=float(ld['lfcSE'][row.position_l,lpos[row.neighbor]]);re=float(rd['lfcSE'][row.position_r,rpos[row.neighbor]])
        fallback=float(k[ti[row.target],gi[row.neighbor]])
        bl=float(a)*.04/(.04+le**2) if not row.neighboring_gene_KD_l and np.isfinite(le) and le>0 else fallback
        br=float(b)*.04/(.04+re**2) if not row.neighboring_gene_KD_r and np.isfinite(re) and re>0 else fallback
        rows.append({'target':row.target,'neighbor':row.neighbor,'condition':row.culture_condition,'distance':row.distance,
                     'base_left':bl,'base_right':br,'left_log2fc':float(np.clip(a,-3,3)),'right_log2fc':float(np.clip(b,-3,3)),
                     'left_q':float(ld['adj_p_value'][row.position_l,lpos[row.neighbor]]),
                     'right_q':float(rd['adj_p_value'][row.position_r,rpos[row.neighbor]])})
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'promoter_donor_observations.csv',index=False)
    target=frame.target.to_numpy();left=frame.left_log2fc.to_numpy();right=frame.right_log2fc.to_numpy();s=ramp(frame.distance.to_numpy())
    bl=frame.base_left.to_numpy();br=frame.base_right.to_numpy()
    keys=[frame.target,frame.neighbor]
    support_l=((frame.left_q<.05)&(frame.left_log2fc<0));support_r=((frame.right_q<.05)&(frame.right_log2fc<0))
    support_l=(support_l&support_l.groupby(keys).transform('sum').ge(2)).to_numpy()
    support_r=(support_r&support_r.groupby(keys).transform('sum').ge(2)).to_numpy()
    folds=v3.target_folds(target);grid=[0.,.25,.5,1.];records=[];nested=[]
    for outer in range(5):
        test=folds==outer;train=~test;inner=v3.target_folds(target[train],3);indices=np.flatnonzero(train)
        candidates=[]
        for alpha in grid:
            total=0.;n=0
            for part in range(3):
                tune=np.zeros(len(target),bool);tune[indices[inner==part]]=True;fit=train&~tune
                pl,pr=predictions(left,right,s,fit,alpha,bl,br,support_l,support_r);total+=loss(pl,pr,left,right,tune)*tune.sum();n+=tune.sum()
            candidates.append((total/n,alpha))
        chosen=min(candidates)[1]
        for alpha in grid:
            pl,pr=predictions(left,right,s,train,alpha,bl,br,support_l,support_r)
            record={'fold':outer,'alpha':alpha,'observations':int(test.sum()),'loss':loss(pl,pr,left,right,test)};records.append(record)
            if alpha==chosen:nested.append({**record,'train_targets':sorted(set(target[train])),'test_targets':sorted(set(target[test]))})
    def combine(group):return float(np.average([r['loss'] for r in group],weights=[r['observations'] for r in group]))
    baseline=combine([r for r in records if r['alpha']==0]);heldout=combine(nested)
    summary=[{'alpha':alpha,'held_target_loss':combine([r for r in records if r['alpha']==alpha])} for alpha in grid]
    selected=min(summary,key=lambda r:(r['held_target_loss'],r['alpha']))['alpha']
    agree=float(np.mean((left<0)&(right<0)))
    accepted=len(set(target))>=10 and heldout<.9*baseline and selected>0 and agree>=.7
    mu_l=estimate(left/np.maximum(s,.05),s**2,np.ones(len(target),bool));mu_r=estimate(right/np.maximum(s,.05),s**2,np.ones(len(target),bool))
    magnitude=(mu_l+mu_r)/2*selected if accepted else 0.
    result={'accepted':bool(accepted),'selected_alpha':selected,'fitted_close_neighbor_log2fc':magnitude,
            'held_target_loss':heldout,'v4_donor_or_k562_baseline_loss':baseline,'relative_v4_neighbor_loss':heldout/baseline,
            'both_donors_negative_fraction':agree,'observations':len(frame),'targets':len(set(target)),
            'nested_folds':nested,'candidate_summary':summary,'script_sha256':build.sha256(Path(__file__)),
            'source_sha256':{str(p.relative_to(ROOT)):build.sha256(p) for p in [ROOT/'analysis/upstream_review/promoter_neighbor_pairs.csv',ROOT/'experiments/v4_lineage_response/analysis/source_components.npz']+
                [ROOT/'reference/cd4_probe'/(n+ext) for n in DONORS for ext in ['.csv','.npz']]},
            'policy':'Require significant negative source response in at least two source conditions; Target-held-out, donor-independent distance ramp compared with donor-specific V4 CD4/K562 estimator; minimum 10 targets, >=10% nested MSE improvement, >=70% both-donor negative signs; apply only to CD4-matched context A and leave existing stronger negative effects intact.',
            'scope':'A hypothesis for conserved guide/promoter-neighbor suppression; boundaries are GENCODE gene TSS, not guide positions; no B/C collateral-effect validation.'}
    build.write_json(OUT/'promoter_selection.json',result)
    pd.DataFrame(records).to_csv(OUT/'promoter_outer_scores.csv',index=False)
    print(json.dumps({k:v for k,v in result.items() if k not in ['nested_folds','source_sha256']},indent=2),flush=True)


if __name__=='__main__':main()
