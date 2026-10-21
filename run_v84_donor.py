"""V84 project-unseen donor count test; controls-only prediction and paired scoring."""
import argparse,gc,json,hashlib,time
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import v84_intervention_operator as m
from v57_pie import read,write,sha
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
from audit_v76_counts import spot_ranks

def setup():
 folder=m.ROOT/'experiments/v83_context_expansion/analysis/donor';receipt=read(folder/'extraction.json');assert sha(folder/'counts.npz')==receipt['counts_sha256'];assert sha(folder/'cells.csv')==receipt['cells_sha256'];cohort=read(folder/'cohort.json');assert sha(folder/'cohort.json')==receipt['cohort_sha256'];p=m.A/'donor_predictions.npz';assert sha(p)==read(p.with_suffix('.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:prediction={k:z[k].copy() for k in z.files}
 return sparse.load_npz(folder/'counts.npz'),pd.read_csv(folder/'cells.csv'),cohort,prediction

def counts(offset,chunk,chunks):
 x,meta,cohort,pred=setup();genes=pred['genes'];targets=pred['targets'];reference=x[cohort['right']];assert (meta.iloc[cohort['left']].gene_target=='Non-Targeting').all() and (meta.iloc[cohort['right']].gene_target=='Non-Targeting').all();cc,cj=profile_jackknife_chunked(reference);atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();selection=read(m.A/'selection.json');arms=['v66','program25','program25_shuffled']
 if selection['selected'] and selection['selected']!='program25':arms.extend([selection['selected'],selection['selected']+'_shuffled'])
 dest=m.A/'donor_counts'/str(offset);dest.mkdir(parents=True,exist_ok=True)
 for ci,case in enumerate(cohort['cases']):
  if ci%chunks!=chunk:continue
  view=next(v for v in case['views'] if v['offset']==offset);target=case['target'];ti=list(targets).index(target);own=list(genes).index(target);actual=x[view['actual']];template=x[view['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert (meta.iloc[view['actual']].gene_target==target).all();truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
  for arm in arms:
   path=dest/f'{ci}__{arm}.json';arr=path.with_suffix('.npz')
   if path.exists():assert sha(arr)==read(path)['arrays_sha256'];continue
   start=time.time();count,diag=emit(template,np.arange(400),pred[arm][ti],own,view['seed'],'residual10',(libs,pred['cp'],pred['bulk']),atlas);assert count.shape==(400,len(genes)) and np.issubdtype(count.dtype,np.integer) and (count.data>=0).all();np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
   pr,pj=profile_jackknife_chunked(count);dd=de(count,reference,genes);scores=direction_metrics(dd,truth,genes,own);coordinates=spot_ranks(count,reference,dd);axis=np.arange(len(genes))!=own;err=float(np.mean((pr[axis]-tp[axis])**2));h=hashlib.sha256()
   for z in [count.data,count.indices,count.indptr]:h.update(np.ascontiguousarray(z).tobytes())
   np.savez_compressed(arr,pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
   record=dict(target=target,case_index=ci,arm=arm,seed_offset=offset,seed=view['seed'],prediction_cells=400,actual_cells=400,reference_cells=reference.shape[0],unique_prediction_control_cells=len(set(view['template'])),original_control_library_totals_exact=True,whole_CD4_response_family_excluded=True,query_treated_cells_used_for_inference=False,count_sha256=h.hexdigest(),profile_MSE=err,manual_rank_coordinates_checked=coordinates,**scores,emitter_diagnostics=diag,arrays_sha256=sha(arr),seconds=time.time()-start);write(path,record);print('V84 donor',offset,target,arm,'profile',round(err,6),'FID',record['fidelity'],flush=True);del count;gc.collect()

def assess(offset):
 _,_,cohort,pred=setup();dest=m.A/'donor_counts'/str(offset);selection=read(m.A/'selection.json');arms=['v66','program25','program25_shuffled']
 if selection['selected'] and selection['selected']!='program25':arms.extend([selection['selected'],selection['selected']+'_shuffled'])
 summary={};records=[]
 for arm in arms:
  rows=[];values=[]
  for ci,case in enumerate(cohort['cases']):
   p=dest/f'{ci}__{arm}.json';r=read(p);assert sha(p.with_suffix('.npz'))==r['arrays_sha256'];rows.append(r);records.append(r)
   with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:values.append({k:z[k].copy() for k in ['pred','jack','truth','truth_jack','control','control_jack']})
  panel=panel_metrics(np.array([v['pred'] for v in values]),np.array([v['truth'] for v in values]),values[0]['control'],np.array([v['jack'] for v in values]),np.array([v['truth_jack'] for v in values]),values[0]['control_jack'],pred['genes'],[r['target'] for r in rows],pred['targets'])
  for k in ['fidelity','reach','jaccard','nmae','profile_MSE']:
   vv=[r[k] for r in rows if r.get(k) is not None and np.isfinite(r[k])];panel[k]=float(np.mean(vv)) if vv else None
  panel['nmae_available_targets']=[r['target'] for r in rows if r['nmae'] is not None];summary[arm]=panel
 base=summary['v66'];checks={}
 for arm in [a for a in arms if a!='v66' and not a.endswith('_shuffled')]:
  v=summary[arm];fake=summary[arm+'_shuffled'];checks[arm]=dict(profile=v['profile_MSE']<=base['profile_MSE'],pds=v['pds']>=base['pds']-.005,fidelity=v['fidelity']>=base['fidelity']-.002,reach=v['reach']>=base['reach']-.005,jaccard=v['jaccard']>=base['jaccard']-.0005,nmae_available=v['nmae'] is not None and base['nmae'] is not None and v['nmae']<=1.01*base['nmae'],beats_shuffle=v['profile_MSE']<fake['profile_MSE'] and (v['fidelity']>fake['fidelity'] or v['jaccard']>fake['jaccard']))
  gains=[v['profile_MSE']<base['profile_MSE'],v['pds']>base['pds'],v['fidelity']>base['fidelity'],v['reach']>base['reach'],v['jaccard']>base['jaccard']]
  if v['nmae'] is not None and base['nmae'] is not None:gains.append(v['nmae']<base['nmae'])
  checks[arm]['three_gains']=sum(gains)>=3
 write(m.A/f'donor_assessment_{offset}.json',dict(summary=summary,checks=checks,passes={a:all(c.values()) for a,c in checks.items()},cases=records,seed_offset=offset,primary='program25',selected_on_old_public_development=selection['selected'],new_donor_not_used_for_selection=True,scope='Project-unseen donor in previously used CD4 study. Whole CD4 responses excluded. Technical seeds reuse donor biology. Public raw metrics; no hidden overall forecast.'))
 print(json.dumps(dict(offset=offset,summary=summary,checks=checks)),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['counts','assess']);p.add_argument('--offset',type=int,default=0);p.add_argument('--chunk',type=int,default=0);p.add_argument('--chunks',type=int,default=1);a=p.parse_args()
 with threadpool_limits(limits=2):(counts(a.offset,a.chunk,a.chunks) if a.stage=='counts' else assess(a.offset))
