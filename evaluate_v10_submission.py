"""V10 public development: fixed V4 responses, freeze calibration winner, then confirmation.

No treated batch mixture enters model inputs. Each control-defined stratum has
one common reference and a fixed target-independent template draw per seed.
"""
import argparse,json,gc,time
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from fit_v7_response import read_cache,raw_inputs
from v10_methods import ROOT,A,SEEDS,SETTINGS,emit
from v8_methods import gated,panel_metrics
from v10_methods import response_bank
from v6_methods import moments,de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build


def load_plan():
    plan=json.loads((A/'validation_plan.json').read_text())
    if plan['script_sha256']!=build.sha256(ROOT/'prepare_v10_validation.py'):raise ValueError('Preparation code changed after plan freeze')
    if plan['source_cache_sha256']!=build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'):raise ValueError('Source labels changed')
    return plan


def balanced(frame,fields):
    return frame.groupby(['family','context','candidate'])[fields].mean().groupby(['family','candidate']).mean().groupby('candidate').mean()


def utility(row):
    # Development proxy only. All six metrics remain visible; no hidden anchors are fabricated.
    if row['common_expression_ratio'] is None:return np.nan
    values=[2*(row['pds']-.5),row['fidelity'],row['reach'],row['jaccard'],
            1-np.clip(row['nmae'],0,3),1-np.clip(row['common_expression_ratio'],0,4)]
    return float(np.mean(values)) if np.isfinite(values).all() else np.nan


def candidates():
    return SETTINGS


def count_stage(split,workers=2):
    plan=load_plan()
    if split=='confirmation':
        frozen=json.loads((A/'count_selection.json').read_text())
        assert frozen['confirmation_read']==False and frozen['calibration_scores_sha256']==build.sha256(A/'count_calibration_scores.csv')
    data=read_cache();genes=data['genes'];targets=data['targets'];settings=candidates();atlas,_=load_atlas_model()
    cache=A/'count_cache';cache.mkdir(exist_ok=True);records=[];summaries=[];diagnostics=[]
    inputs={'plan':build.sha256(A/'validation_plan.json'),
            'script':build.sha256(ROOT/'evaluate_v10_submission.py'),'methods':build.sha256(ROOT/'v10_methods.py')}
    for family in ['H1','HCT116','HEK293T']:
        bank,_=response_bank(data,excluded=family)
        x,meta,measured,source=raw_inputs(family,genes)
        for context in [c for c in plan['contexts'] if c['family']==family]:
            pool=x[context['left']];population=moments(pool);query=population[1]*1e6
            prepared=(population[0],population[1],population[2],np.asarray(pool.mean(axis=0)).ravel())
            for seed in SEEDS:
                cells=[c for c in context['cells'] if c['seed']==seed and c['split']==split]
                control=x[cells[0]['reference']];cp,cj=profile_jackknife_chunked(control)
                templates=x[cells[0]['template']]
                assert all(c['template']==cells[0]['template'] and c['reference']==cells[0]['reference'] for c in cells)
                template_profile,template_jack=profile_jackknife_chunked(templates)
                np.testing.assert_allclose(pds(np.broadcast_to(template_profile,(len(cells),len(genes)))[:,~np.isin(genes,targets)],
                                               np.array([profile_jackknife_chunked(x[c['actual']])[0] for c in cells])[:,~np.isin(genes,targets)],cp[~np.isin(genes,targets)]),.5,atol=1e-12)
                null_de=de(templates,control,genes)
                diagnostics.append({'family':family,'context':context['context'],'split':split,'seed':seed,
                                    'constant_profile_pds':.5,'neutral_false_de_calls':int((null_de['q']<.05).sum()),
                                    'neutral_tested_genes':int(null_de['eligible'].sum()),
                                    'control_pool_overlap':len(set(cells[0]['template'])&set(cells[0]['reference']))})
                def evaluate(cell):
                    target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0])
                    file=cache/f'{context["context"]}_{ti:03d}_{seed}_{split}.npz';side=file.with_suffix('.json')
                    if file.exists() and side.exists():
                        rec=json.loads(side.read_text())
                        if rec['inputs']!=inputs or rec['sha256']!=build.sha256(file):raise ValueError('Count evaluation cache differs from frozen inputs')
                        with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
                        return rec['records'],arrays
                    real=x[cell['actual']];real_de=de(real,control,genes);tp,tj=profile_jackknife_chunked(real)
                    arrays={'truth':tp,'truth_jack':tj,'control_null_pred':template_profile,'control_null_jack':template_jack}
                    group=[]
                    def store(name,count,diag):
                        pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj
                        metrics=direction_metrics(de(count,control,genes),real_de,genes,own)
                        libs,p,b,var=moments(count)
                        row={'family':family,'context':context['context'],'target':target,'seed':seed,'split':split,
                             'candidate':name,'actual_cells':len(real.indptr)-1,'generated_cells':count.shape[0],
                             'mean_library':float(libs.mean()),'mean_detected_genes':float(count.getnnz(axis=1).mean()),
                             'total_CPM_variance':float(var.sum()*1e12),**metrics,**diag}
                        group.append(row)
                    null_metrics=direction_metrics(null_de,real_de,genes,own)
                    group.append({'family':family,'context':context['context'],'target':target,'seed':seed,'split':split,
                                  'candidate':'control_null','actual_cells':real.shape[0],'generated_cells':400,**null_metrics})
                    for name,(response,generator) in settings.items():
                        effect=gated(bank[response][ti].copy(),query,own)
                        block,diag=emit(templates,np.arange(400),effect,own,seed+ti,generator,prepared=prepared,atlas=atlas)
                        store(name,block,diag);del block
                    np.savez_compressed(file,**arrays)
                    build.write_json(side,{'inputs':inputs,'sha256':build.sha256(file),'records':[{k:None if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()} for row in group]})
                    print('V10',split,context['context'],seed,target,'counts evaluated',flush=True)
                    return group,arrays
                with ThreadPoolExecutor(max_workers=workers) as pool_executor:
                    evaluated=list(pool_executor.map(evaluate,cells))
                for group,arrays in evaluated:records.extend(group)
                truth=np.array([a['truth'] for _,a in evaluated]);tj=np.array([a['truth_jack'] for _,a in evaluated]);names=[c['target'] for c in cells]
                for name in list(settings)+['control_null']:
                    pp=np.array([a[name+'_pred'] for _,a in evaluated]);pj=np.array([a[name+'_jack'] for _,a in evaluated])
                    score=panel_metrics(pp,truth,cp,pj,tj,cj,genes,names,targets)
                    subset=pd.DataFrame([r for group,_ in evaluated for r in group if r['candidate']==name])
                    means=subset.mean(numeric_only=True).to_dict()
                    row={'family':family,'context':context['context'],'seed':seed,'split':split,'candidate':name,**means,**score}
                    row['utility']=utility(row);summaries.append(row)
                pd.DataFrame(records).to_csv(A/f'count_{split}_scores.csv',index=False)
                pd.DataFrame(summaries).to_csv(A/f'count_{split}_panels.csv',index=False)
                pd.DataFrame(diagnostics).to_csv(A/f'count_{split}_diagnostics.csv',index=False)
        del x;gc.collect()
    summary=balanced(pd.DataFrame(summaries),['utility','pds','fidelity','reach','jaccard','nmae','common_expression_ratio'])
    summary.reset_index().to_csv(A/f'count_{split}_summary.csv',index=False)
    print(summary.to_string(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['calibration','confirmation'],required=True);p.add_argument('--workers',type=int,default=2);args=p.parse_args()
    with threadpool_limits(limits=2):
        count_stage(args.stage,args.workers)
