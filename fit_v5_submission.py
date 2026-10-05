"""Fit V5: unique public H1 coverage, separate mean/pooled responses and calibrated A promoter prior."""
from pathlib import Path
import argparse,json,gc
import numpy as np
import pandas as pd
from scipy import sparse,stats
import build_submission as build
import fit_lineage_submission as v4
import fit_context_ensemble as v3
from calibrate_v5_promoters import ramp,quality
from build_cd4_model import load_subset
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v5_atlas_calibrated/analysis'
DATA=ROOT/'reference/v5_h1'


def moments(x):
    x=x.astype(float);n=x.shape[0];library=np.asarray(x.sum(axis=1)).ravel();mean_library=library.mean()
    p=np.asarray(x.sum(axis=0)).ravel()/library.sum()
    # Delta-method variance for the library-weighted ratio, retaining count/depth covariance.
    residual_second=np.asarray(x.power(2).mean(axis=0)).ravel()-2*p*np.asarray(x.T.dot(library)).ravel()/n+p*p*np.mean(library**2)
    var=np.maximum(residual_second,0)/max(n-1,1)/mean_library**2
    return p,var


def source_inputs():
    report=json.loads((DATA/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:assert build.sha256(DATA/name)==report[key]
    raw=sparse.load_npz(DATA/'selected_counts.npz').astype(float);meta=pd.read_csv(DATA/'selected_cells.csv');genes=pd.read_csv(DATA/'gene_names.csv').gene_name.to_numpy()
    assert meta.source_cell_id.is_unique and (raw.sum(axis=1)>0).all()
    library=np.asarray(raw.sum(axis=1)).ravel();cpm=raw.multiply((1e6/library)[:,None]).tocsr()
    batches=sorted(meta.batch.unique());partition=meta.batch.map({b:i%2 for i,b in enumerate(batches)}).to_numpy()
    controls=meta.target_gene.to_numpy()=='non-targeting';labels=meta.batch.to_numpy();control_stats={}
    for batch in batches:
        rows=np.flatnonzero(controls&(labels==batch));assert len(rows)>=20
        x=cpm[rows];mean=np.asarray(x.mean(axis=0)).ravel();var=np.maximum(np.asarray(x.power(2).mean(axis=0)).ravel()-mean**2,0)/max(len(rows)-1,1)
        bp,bv=moments(raw[rows]);control_stats[batch]=(mean,var,bp,bv)
    target_names=np.array(sorted(set(meta.target_gene)-{'non-targeting'}));parts=[]
    for part in [0,1,None]:
        values=[];errors=[];bases=[];qvalues=[];numbers=[];bulk_values=[];bulk_errors=[]
        for target in target_names:
            selected=meta.target_gene.to_numpy()==target
            if part is not None:selected&=partition==part
            rows=np.flatnonzero(selected);n=len(rows)
            if n<20:raise ValueError(f'H1 {target} part {part} has only {n} selected cells')
            x=cpm[rows];mean=np.asarray(x.mean(axis=0)).ravel();vm=np.maximum(np.asarray(x.power(2).mean(axis=0)).ravel()-mean**2,0)/max(n-1,1)
            base=np.zeros(len(genes));bv=np.zeros(len(genes));bbase=np.zeros(len(genes));bbvar=np.zeros(len(genes))
            for batch,num in meta.iloc[rows].batch.value_counts().items():
                cm,cv,bp,bvar=control_stats[batch];base+=num/n*cm;bv+=(num/n)**2*cv
                # Match the perturbation's batch library contribution in pooled control space.
                fraction=library[rows[labels[rows]==batch]].sum()/library[rows].sum()
                bbase+=fraction*bp;bbvar+=fraction**2*bvar
            lfc=np.clip(np.log2((mean+1)/(base+1)),-3,3);se=np.sqrt(vm/(mean+1)**2+bv/(base+1)**2)/np.log(2)
            bp,bpv=moments(raw[rows]);delta=np.clip(np.log1p(50000*bp)-np.log1p(50000*bbase),-3,3)
            bse=50000*np.sqrt(bpv/(1+50000*bp)**2+bbvar/(1+50000*bbase)**2)
            valid=(base>5)&(genes!=target)&(se>0);control_rows=np.flatnonzero(controls if part is None else controls&(partition==part))
            q=np.ones(len(genes));columns=np.flatnonzero(valid)
            for start in range(0,len(columns),256):
                cols=columns[start:start+256]
                q[cols]=stats.mannwhitneyu(cpm[rows][:,cols].toarray(),cpm[control_rows][:,cols].toarray(),axis=0,method='asymptotic').pvalue
            if valid.any():q[valid]=stats.false_discovery_control(q[valid],method='bh')
            values.append(lfc);errors.append(se);bases.append(base);qvalues.append(q);numbers.append(n);bulk_values.append(delta);bulk_errors.append(bse)
        parts.append({'values':np.array(values),'errors':np.array(errors),'baselines':np.array(bases),'q':np.array(qvalues),'n':np.array(numbers),
                      'bulk_values':np.array(bulk_values),'bulk_errors':np.array(bulk_errors)})
        print('H1 contrast part',part,'computed for',len(target_names),'targets',flush=True)
    controls_record={'left_batches':[str(b) for i,b in enumerate(batches) if i%2==0],'right_batches':[str(b) for i,b in enumerate(batches) if i%2==1],
                     'left_control_ids':meta.loc[controls&(partition==0),'source_cell_id'].tolist(),'right_control_ids':meta.loc[controls&(partition==1),'source_cell_id'].tolist()}
    assert not set(controls_record['left_control_ids'])&set(controls_record['right_control_ids'])
    return target_names,genes,parts,report,controls_record


def shrink(values,errors,gain):return np.clip(gain*np.sign(values)*np.maximum(np.abs(values)-errors,0),-3,3)


def select_bulk(targets,genes,parts):
    left,right=parts[:2];l=left['bulk_values'];r=right['bulk_values'];le=left['bulk_errors'];re=right['bulk_errors']
    valid=(left['baselines']>5)&(right['baselines']>5)&(genes[None,:]!=targets[:,None])&(le>0)&(re>0)
    sl=le*np.sqrt(left['n'][:,None]/400);sr=re*np.sqrt(right['n'][:,None]/400)
    gains=[.3,.5,.7,1.];folds=v3.target_folds(targets);outer=[];nested=[]
    def evaluate(gain,mask):return v4.evaluate(shrink(l,le,gain),shrink(r,re,gain),l,r,left['q'],right['q'],valid,sl,sr,mask)
    for fold in range(5):
        test=folds==fold;train=~test
        choices=[{'gain':gain,**evaluate(gain,train)} for gain in gains]
        baseline=next(x for x in choices if x['gain']==.7)
        eligible=[x for x in choices if x['loss']<=baseline['loss']+1e-12 and x['relative_zero_loss']<=1.02]
        chosen=max(eligible,key=lambda x:(x['selection_score'],-x['gain']))
        for gain in gains:
            record={'fold':fold,'gain':gain,**evaluate(gain,test)};outer.append(record)
            if gain==chosen['gain']:nested.append({**record,'train_targets':targets[train].tolist(),'test_targets':targets[test].tolist()})
    def aggregate(rows):
        n=sum(x['contrasts'] for x in rows)
        result={key:sum(x[key]*x['contrasts'] for x in rows)/n for key in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy']}
        result['relative_zero_loss']=result['loss']/result['zero_loss'];return result
    summary=[{'gain':gain,**aggregate([r for r in outer if r['gain']==gain])} for gain in gains]
    baseline=next(x for x in summary if x['gain']==.7);procedure=aggregate(nested)
    supported=procedure['loss']<=baseline['loss']+1e-12 and procedure['selection_score']>=baseline['selection_score']-1e-12
    eligible=[x for x in summary if x['loss']<=baseline['loss']+1e-12 and x['relative_zero_loss']<=1.02]
    chosen=max(eligible,key=lambda x:(x['selection_score'],-x['gain']))['gain'] if supported else .7
    result={'selected_gain':chosen,'nested_procedure':procedure,'baseline_gain_0_7':baseline,'candidate_summary':summary,'nested_folds':nested,
            'changed_amplitude_supported':bool(supported and chosen!=.7),
            'scope':'Five held-target folds using disjoint source batches. Pooled log1p(50000*p) error and source direction/overlap proxy; not a hidden score. The 0.7 baseline amplitude is a conservative modeling reference, not a published separate-bulk V4 estimator.'}
    build.write_json(OUT/'bulk_selection.json',result);pd.DataFrame(outer).to_csv(OUT/'bulk_outer_scores.csv',index=False)
    print('Selected pooled response gain',chosen,'nested',procedure,flush=True);return chosen


def h1_arrays(ht,hg,parts,genes,targets,hsetting,transfer,bulk_gain):
    full=parts[2];gp={g:i for i,g in enumerate(genes)};si=np.array([i for i,g in enumerate(hg) if g in gp]);di=np.array([gp[hg[i]] for i in si])
    own=np.array([int(np.flatnonzero(hg==t)[0]) for t in ht]);knockdown=full['values'][np.arange(len(ht)),own]
    own_se=full['errors'][np.arange(len(ht)),own];own_base=full['baselines'][np.arange(len(ht)),own]
    accepted=(full['n']>=40)&(knockdown<=np.log2(.7))&(own_base>5)&(knockdown<-1.96*own_se)
    valid=(full['baselines']>5)&(full['errors']>0)&(hg[None,:]!=ht[:,None])&accepted[:,None]
    bvalid=valid&(full['bulk_errors']>0)
    h=np.zeros((len(targets),len(genes)),np.float32);hp=np.zeros_like(h);bulk=np.zeros_like(h);bulkp=np.zeros_like(h);mask=np.zeros_like(h,bool);available=np.zeros(len(targets),bool)
    conditions=np.repeat('Rest',len(ht));qc=[]
    for i,target in enumerate(ht):qc.append({'target_gene':target,'selected_cells':int(full['n'][i]),'source_control_mean_CPM':float(own_base[i]),
        'own_log2fc':float(knockdown[i]),'own_log2fc_SE':float(own_se[i]),'direct_response_accepted':bool(accepted[i])})
    for ti,target in enumerate(targets):
        train=ht!=target
        prior=v4.response_prior(full['values'],full['errors'],valid,conditions,train)['Rest']
        bp=v4.response_prior(full['bulk_values'],full['bulk_errors'],bvalid,conditions,train)['Rest']
        idx=np.flatnonzero(ht==target);ok=bool(len(idx) and accepted[idx[0]])
        gain=hsetting['prior_gain'] if ok else transfer['prior_gain']
        hp[ti,di]=gain*prior[si];bulkp[ti,di]=gain*bp[si]
        if ok:
            j=idx[0];h[ti,di]=shrink(full['values'][j,si],full['errors'][j,si],hsetting['specific_gain'])
            bulk[ti,di]=shrink(full['bulk_values'][j,si],full['bulk_errors'][j,si],bulk_gain)
            mask[ti,di]=valid[j,si];available[ti]=True
    h[~mask]=0;bulk[~mask]=0
    pd.DataFrame(qc).to_csv(OUT/'h1_target_QC.csv',index=False)
    return h,hp,bulk,bulkp,mask,available


def prepare():
    OUT.mkdir(parents=True,exist_ok=True);manifest,genes,targets=build.read_spec(ROOT)
    previous=ROOT/'experiments/v4_lineage_response/analysis';model=json.loads((previous/'reference_model.json').read_text())
    assert build.sha256(previous/'reference_model_arrays.npz')==model['model_arrays_sha256']
    assert build.sha256(previous/'source_components.npz')==model['source_components_sha256']
    for path,digest in model['source_data_sha256'].items():assert build.sha256(ROOT/path)==digest,path
    with np.load(previous/'source_components.npz') as z:components={key:z[key].copy() for key in z.files}
    assert np.array_equal(components['genes'],genes) and np.array_equal(components['targets'],targets)
    saved=OUT/'h1_contrasts.npz';metadata=OUT/'h1_contrasts.json'
    source_report=json.loads((DATA/'extraction.json').read_text())
    if metadata.exists():
        record=json.loads(metadata.read_text());assert record['source_extraction_sha256']==build.sha256(DATA/'extraction.json') and record['script_sha256']==build.sha256(Path(__file__))
        assert record['arrays_sha256']==build.sha256(saved)
        with np.load(saved) as z:
            ht=z['targets'].copy();hg=z['genes'].copy();parts=[{key:z[str(i)+'_'+key].copy() for key in ['values','errors','baselines','q','n','bulk_values','bulk_errors']} for i in range(3)]
    else:
        ht,hg,parts,source_report,control_record=source_inputs()
        np.savez_compressed(saved,targets=np.asarray(ht,dtype=str),genes=np.asarray(hg,dtype=str),**{str(i)+'_'+key:value for i,p in enumerate(parts) for key,value in p.items()})
        record={'source_extraction_sha256':build.sha256(DATA/'extraction.json'),'arrays_sha256':build.sha256(saved),'script_sha256':build.sha256(Path(__file__)),
                'targets':len(ht),'genes':len(hg),'globally_disjoint_batch_controls':control_record}
        build.write_json(metadata,record)
    old_out=v4.OUT
    try:v4.OUT=OUT;hsetting=v4.select_h1(ht,hg,parts)
    finally:v4.OUT=old_out
    transfer=model['cross_lineage_settings'];bulk_gain=select_bulk(ht,hg,parts)
    h,hp,hb,hbp,hm,ha=h1_arrays(ht,hg,parts,genes,targets,hsetting,transfer,bulk_gain)
    effects,records=v4.assemble_effects(components['cd4'],components['cd4_specific'],components['cd4_mask'],components['k562'],components['k562_mask'],
                h,hp,hm,ha,transfer,components['context_mean_cpm'],genes,targets)
    promoter=json.loads((OUT/'promoter_selection.json').read_text());pairs=pd.read_csv(ROOT/'analysis/upstream_review/promoter_neighbor_pairs.csv')
    gate=np.clip(components['context_mean_cpm'][0]/5,0,1);gi={g:i for i,g in enumerate(genes)};ti={t:i for i,t in enumerate(targets)};changes=[]
    pm,pdta=load_subset(ROOT/'reference/cd4_probe','merged');pgi={g:i for i,g in enumerate(pdta['genes'])}
    for row in pairs.itertuples():
        if row.neighbor not in pgi:continue
        ix=np.flatnonzero(quality(pm)&(pm.target_contrast_gene_name.to_numpy()==row.target))
        j=pgi[row.neighbor];valid_source=np.isfinite(pdta['log_fc'][ix,j])&np.isfinite(pdta['adj_p_value'][ix,j])
        negative=valid_source&(pdta['log_fc'][ix,j]<0)&(pdta['adj_p_value'][ix,j]<.05)
        if len(set(pm.iloc[ix[negative]].culture_condition))<2:continue
        if promoter['accepted']:
            i,j=ti[row.target],gi[row.neighbor];prior=promoter['fitted_close_neighbor_log2fc']*ramp(row.distance)*gate[j]*gate[gi[row.target]]
            before=float(effects[0,i,j]);effects[0,i,j]=min(before,prior)
            if effects[0,i,j]!=before:changes.append({'context':'A','target_gene':row.target,'neighbor':row.neighbor,'distance':row.distance,'before':before,'after':float(effects[0,i,j])})
    pd.DataFrame(changes,columns=['context','target_gene','neighbor','distance','before','after']).to_csv(OUT/'promoter_effect_changes.csv',index=False)
    # Direct H1 pooled responses are independently measured; other-source responses use a conservative transformed prior.
    bulk_delta=np.zeros_like(effects)
    for ci,context in enumerate('ABC'):
        original=sparse.csr_matrix(__import__('anndata').read_h5ad(ROOT/f'context_{context}.h5ad').X)
        pooled=np.asarray(original.sum(axis=0)).ravel().astype(float);pooled/=pooled.sum()
        baseline=50000*pooled
        if ci==0:
            bulk_delta[ci]=np.log1p(baseline[None,:]*np.exp2(effects[ci]))-np.log1p(baseline[None,:])
        else:
            transferred=transfer['cd4_coefficient']*components['cd4_specific'][ci]+transfer['k562_coefficient']*components['k562']
            target_gate=np.clip(components['context_mean_cpm'][ci,[gi[t] for t in targets]]/5,0,1)
            gene_gate=np.clip(components['context_mean_cpm'][ci]/5,0,1)
            translated=np.log1p(baseline[None,:]*np.exp2(np.clip(transferred,-3,3)))-np.log1p(baseline[None,:])
            bulk_delta[ci]=(np.where(hm,hb,bulk_gain*translated)+hbp)*gene_gate[None,:]*target_gate[:,None]
        for i,target in enumerate(targets):bulk_delta[ci,i,gi[target]]=0
        del original;gc.collect()
    np.savez_compressed(OUT/'reference_model_arrays.npz',effects=effects,bulk_delta=bulk_delta,genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str))
    records.to_csv(OUT/'target_response_sources.csv',index=False)
    np.savez_compressed(OUT/'h1_components.npz',h1=h,h1_prior=hp,h1_bulk=hb,h1_bulk_prior=hbp,h1_mask=hm,h1_available=ha)
    result={'model':'V5 atlas-calibrated lineage response with individual-cell dual moments','v4_parent_model_sha256':build.sha256(previous/'reference_model.json'),
        'v4_parent_arrays_sha256':model['model_arrays_sha256'],'v4_source_components_sha256':model['source_components_sha256'],
        'h1_source_extraction':source_report,'h1_source_extraction_sha256':build.sha256(DATA/'extraction.json'),
        'h1_settings':hsetting,'bulk_gain':bulk_gain,'cd4_settings_unchanged':model['cd4_settings'],'cross_lineage_settings':transfer,
        'public_H1_targets':len(ht),'usable_direct_H1_targets':int(ha.sum()),'promoter_changes_A_only':len(changes),
        'model_arrays_sha256':build.sha256(OUT/'reference_model_arrays.npz'),'h1_components_sha256':build.sha256(OUT/'h1_components.npz'),
        'script_sha256':build.sha256(Path(__file__)),'reused_v4_fitter_sha256':build.sha256(ROOT/'fit_lineage_submission.py'),
        'source_original_inputs':model['context_original_inputs'],'promoter_selection_sha256':build.sha256(OUT/'promoter_selection.json'),
        'h1_contrasts_sha256':build.sha256(saved),'h1_contrasts_metadata_sha256':build.sha256(metadata),
        'limitations':['No hidden response validation or promised score improvement.','H1 source sampling is bounded and source ordered, not all released cells.',
          'Source batch/target holdouts do not guarantee anonymous-context transfer.','Pooled-response amplitude is calibrated in H1 and used as a transfer assumption for B/C.',
          'Promoter prior uses 17 selected CD4 targets; applied only to A.','Cell library sizes are held to original controls; treatment-dependent depth changes are not modeled.']}
    build.write_json(OUT/'reference_model.json',result);print(json.dumps({k:result[k] for k in ['public_H1_targets','usable_direct_H1_targets','h1_settings','bulk_gain','promoter_changes_A_only']},indent=2),flush=True)


if __name__=='__main__':prepare()
