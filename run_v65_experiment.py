"""V65 fixed response screen and paired integer-count diagnostics against V63."""
import argparse,gc,json,time,shutil,subprocess,sys
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v65_closed_response as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from audit_transfer_validation import profile,profile_jackknife_chunked,pds
from v8_methods import panel_metrics
from v12_cpu_metrics import de,direction_metrics
from v16_methods import emit

def cases(ctx):
 seen=set();out=[]
 for i,c in enumerate(ctx['cells']):
  if c['target'] not in seen:out.append(i);seen.add(c['target'])
 return out

def load_predictions(name):
 p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def means():
 m.check();d=m.old.data();genes=d['genes'];targets=d['targets'];own=np.array([list(genes).index(t) for t in targets]);records=[];panels=[]
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name);z=load_predictions(name);ref=profile(sparse.load_npz(folder/'reference.npz'));pr={a:[] for a in ('v24','v63')+m.ARMS};truth=[];names=[]
  from v6_methods import moments
  ref_cpm=moments(sparse.load_npz(folder/'reference.npz'))[1]*1e6
  for i in cases(ctx):
   cell=ctx['cells'][i];ti=list(targets).index(cell['target']);actual=sparse.load_npz(folder/f'{i}_actual.npz');tp=profile(actual);truth.append(tp);names.append(cell['target']);axis=np.arange(len(genes))!=own[ti]
   am=moments(actual)[1]*1e6;valid=z['measured']&(ref_cpm>5)&axis;true_lfc=np.clip(np.log2((am+.1)/(ref_cpm+.1)),-3,3)
   for arm in pr:
    bp=z['bulk']*np.exp2(z[arm][ti]);bp[own[ti]]*=.1;bp/=bp.sum();pp=np.log1p(50000*bp);pr[arm].append(pp)
    cp=np.nan_to_num(z['query'])*np.exp2(z[arm][ti]);cp[own[ti]]*=.1;cp*=1e6/cp.sum();lfc=np.clip(np.log2((cp+.1)/(ref_cpm+.1)),-3,3);delta=lfc[valid]-true_lfc[valid]
    records.append(dict(source=name,study=key(name),target=cell['target'],case_index=i,arm=arm,profile_error=float(np.sum((pp[axis]-tp[axis])**2)),lfc_square=float(delta@delta),lfc_absolute=float(abs(delta).sum())))
  for arm,pp in pr.items():
   axis=z['measured']&~np.isin(genes,targets);panels.append(dict(source=name,study=key(name),arm=arm,pds=pds(np.asarray(pp)[:,axis],np.asarray(truth)[:,axis],ref[axis])))
  print('V65 mean',name,len(names),flush=True)
 summaries={}
 for arm in m.ARMS:
  studies={}
  for study in ('H1','Orion','CD4'):
   sr={}
   for metric in ('profile_error','lfc_square','lfc_absolute'):
    ratios=[]
    for name in sorted(set(r['source'] for r in records if r['study']==study)):
     base=sum(r[metric] for r in records if r['source']==name and r['arm']=='v63');new=sum(r[metric] for r in records if r['source']==name and r['arm']==arm);ratios.append(new/base)
    sr[metric+'_ratio']=float(np.mean(ratios))
   pp=[p['pds'] for p in panels if p['study']==study and p['arm']==arm];bb=[p['pds'] for p in panels if p['study']==study and p['arm']=='v63'];sr['pds_change']=float(np.mean(pp)-np.mean(bb));studies[study]=sr
  summary={k:float(np.mean([r[k] for r in studies.values()])) for k in next(iter(studies.values()))};summary['studies']=studies
  summary['eligible']=summary['profile_error_ratio']<1 and summary['lfc_square_ratio']<=1.02 and summary['lfc_absolute_ratio']<=1.02 and summary['pds_change']>=-.01 and max(s['profile_error_ratio'] for s in studies.values())<=1.03
  summaries[arm]=summary
 eligible=[a for a in m.ARMS if summaries[a]['eligible']];selected=min(eligible,key=lambda a:summaries[a]['profile_error_ratio']) if eligible else None
 write(m.A/'mean_results.json',dict(records=records,panels=panels,summary=summaries,selected=selected,selection_frozen_before_count_results=True,plan_sha256=sha(m.A/'plan.json')))
 print(json.dumps(dict(summary=summaries,selected=selected),indent=2),flush=True)

def counts(source):
 m.check();selected=read(m.A/'mean_results.json')['selected'];assert selected,'No eligible response candidate'
 atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();d=m.old.data();genes=d['genes'];targets=d['targets'];dest=m.A/'counts';dest.mkdir(exist_ok=True);arms=('v24','v63',selected,selected+'_shuffled')
 for ctx in contexts():
  name=ctx['family']
  if name!=source:continue
  folder,_=checked_folder(name);z=load_predictions(name);reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
  with np.load(folder/'effects.npz',allow_pickle=False) as raw:cp=raw['cp'].copy();bulk=raw['bulk'].copy()
  for offset in (0,620000):
   outputs={a:[] for a in arms}
   for i in cases(ctx):
    cell=ctx['cells'][i];ti=list(targets).index(cell['target']);own=list(genes).index(cell['target']);template=sparse.load_npz(folder/f'{i}_template.npz');actual=sparse.load_npz(folder/f'{i}_actual.npz');libs=template.sum(1).A1.astype(np.int64);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual);axis=np.arange(len(genes))!=own
    for arm in arms:
     path=dest/f'{name}__{i}__{offset}__{arm}.json';arr=path.with_suffix('.npz');start=time.time()
     if not path.exists():
      reused=None
      if arm in ('v24','v63'):
       stage='development' if i<4 and offset==0 else 'additional';oldarm='guarded_assay' if arm=='v63' else 'v24';old=m.ROOT/'experiments/v63_identity_guard/analysis'/stage/f'{name}__{i}__{offset}__{oldarm}.json'
       if old.exists():
        r=read(old);assert sha(old.with_suffix('.npz'))==r['arrays_sha256'];shutil.copyfile(old.with_suffix('.npz'),arr);r=dict(r,arm=arm,reused_source=str(old.relative_to(m.ROOT)),reused_receipt_sha256=sha(old));reused=True
      if not reused:
       c,diag=emit(template,np.arange(len(libs)),z[arm][ti],own,cell['seed']+ti+offset,'residual10',(libs,cp,bulk),atlas);np.testing.assert_array_equal(c.sum(1).A1,libs);assert c.dtype.kind in 'iu' and c.data.min(initial=0)>=0
       dd=de(c,reference,genes);pr,pj=profile_jackknife_chunked(c);metrics=direction_metrics(dd,truth,genes,own)
       np.savez_compressed(arr,pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
       r=dict(family=name,study=key(name),target=cell['target'],case_index=i,seed_offset=offset,arm=arm,cells=len(libs),row_libraries_exact=True,**metrics,**diag,arrays_sha256=sha(arr),seconds=time.time()-start);del c;gc.collect()
      write(path,r)
     r=read(path);assert sha(arr)==r['arrays_sha256']
     with np.load(arr,allow_pickle=False) as rr:pr=rr['pred'].copy();pj=rr['jack'].copy();np.testing.assert_array_equal(rr['truth'],tp)
     r['profile_square_error']=float(np.sum((pr[axis]-tp[axis])**2));outputs[arm].append((pr,pj,tp,tj,r));print('V65 counts',name,i,offset,arm,round(time.time()-start,1),flush=True)
   for arm,vals in outputs.items():
    pp,jp,tt,jt=(np.asarray([v[k] for v in vals]) for k in range(4));rows=[v[4] for v in vals];panel=panel_metrics(pp,tt,cc,jp,jt,cj,genes,[r['target'] for r in rows],targets)
    for metric in ('fidelity','reach','jaccard','nmae','profile_square_error'):
     vv=[r[metric] for r in rows if r.get(metric) is not None and np.isfinite(r[metric])];panel[metric]=float(np.mean(vv)) if vv else None
    write(dest/f'{name}__{offset}__{arm}__panel.json',dict(family=name,study=key(name),seed_offset=offset,arm=arm,**panel))

def parallel():
 selected=read(m.A/'mean_results.json')['selected'];assert selected
 jobs=[c['family'] for c in contexts()];active={}
 while jobs or active:
  while jobs and len(active)<2:
   name=jobs.pop(0);log=(m.OUT/(name+'_counts.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'run_v65_experiment.py'),'counts','--source',name],stdout=log,stderr=subprocess.STDOUT);active[name]=(p,log);print('Started V65',name,flush=True)
  for name,(p,log) in list(active.items()):
   if p.poll() is not None:log.close();assert p.returncode==0,(name,p.returncode);del active[name];print('Completed V65',name,flush=True)
  if jobs or active:time.sleep(3)
 assess()

def assess():
 selected=read(m.A/'mean_results.json')['selected'];panels=[read(p) for p in (m.A/'counts').glob('*__panel.json')];arms=('v24','v63',selected,selected+'_shuffled');assert len(panels)==40,len(panels)
 metrics=('pds','fidelity','reach','jaccard','nmae','profile_square_error');summaries={};studies={}
 for arm in arms:
  studies[arm]={}
  for study in ('H1','Orion','CD4'):
   studies[arm][study]={}
   for k in metrics:
    vals=[p[k] for p in panels if p['arm']==arm and p['study']==study and p.get(k) is not None and np.isfinite(p[k])];studies[arm][study][k]=float(np.mean(vals)) if vals else None
  summaries[arm]={k:float(np.mean([s[k] for s in studies[arm].values() if s[k] is not None])) for k in metrics}
 b=summaries['v63'];v=summaries[selected];sh=summaries[selected+'_shuffled'];g=read(m.A/'plan.json')['count_guards']
 tests=dict(profile=v['profile_square_error']<b['profile_square_error'],pds=v['pds']>=b['pds']-g['pds_loss'],fidelity=v['fidelity']>=b['fidelity']-g['fidelity_loss'],reach=v['reach']>=b['reach']-g['reach_loss'],jaccard=v['jaccard']>=b['jaccard']-g['jaccard_loss'],nmae=v['nmae']<=b['nmae']*g['nmae_ratio'],beats_shuffle_profile=v['profile_square_error']<sh['profile_square_error'],beats_shuffle_pds=v['pds']>=sh['pds'])
 ratio=float(np.mean([v1['profile_square_error']/studies['v63'][k]['profile_square_error'] for k,v1 in studies[selected].items()]))
 write(m.A/'count_assessment.json',dict(complete=True,release_supported=all(tests.values()),checks=tests,summary=summaries,studies=studies,paired_equal_study_profile_error_ratio=ratio,selected=selected,unique_target_context_pairs=30,seeds=[0,620000],comparisons=240,D3_opened=False,private_score=None,scope='Reused public cohorts; two seeds are technical repeats. Equal context/seed within study then equal studies; NMAE undefined studies omitted and explicitly visible in studies. No official anchors.'))
 print(json.dumps(dict(checks=tests,summary=summaries,profile_ratio=ratio),indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['means','counts','parallel','assess']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='counts':counts(a.source)
  else:globals()[a.stage]()
