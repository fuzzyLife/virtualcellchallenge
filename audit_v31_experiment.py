"""Independent raw-guide reconstruction and development arithmetic audit for V31."""
from pathlib import Path
import json
import numpy as np
from scipy import special
from threadpoolctl import threadpool_limits
import v31_methods as m
import build_submission as b

def moments(x,ids):
 a=x[ids].toarray().astype(float);a*=1e6/a.sum(1,keepdims=True);return a.mean(0),a.var(0)/max(len(ids)-1,1)

def source_audit(condition,ci,meta_record,t):
 genes=m.data()['genes'];x,meta,measured,source=m.raw(condition,genes);labels=meta.gene_target.to_numpy(str);lanes=meta['sample'].to_numpy(str);gg=meta.guide_id.to_numpy(str);allguide=meta_record['guides'];targets=sorted(set(v['target'] for v in allguide));selected=targets[:2]+targets[-1:];output=[]
 for target in selected:
  records=[r for r in allguide if r['target']==target];rows=[];eff=[];var=[];valid=[];n=[];kds=[];own=int(np.flatnonzero(genes==target)[0]);seen=set()
  for rec in records:
   ids=np.array(rec['original_rows'],int);assert len(ids)==rec['n'] and len(set(ids))==len(ids) and not seen&set(ids);seen.update(ids);assert (labels[ids]==target).all() and (gg[ids]==rec['guide']).all();mu,mv=moments(x,ids);cp=np.zeros(len(genes));cv=cp.copy()
   for lane in sorted(set(lanes[ids])):
    control=np.flatnonzero((labels=='Non-Targeting')&(lanes==lane));assert len(control)>=20;c,v=moments(x,control);w=float(np.mean(lanes[ids]==lane));cp+=w*c;cv+=w*w*v
   lfc=np.clip(np.log2((mu+.1)/(cp+.1)),-3,3);variance=np.maximum((mv/(mu+.1)**2+cv/(cp+.1)**2)/np.log(2)**2,.0025);kd=float(np.clip(1-(mu[own]+.1)/(cp[own]+.1),0,.99)) if cp[own]>5 else None
   if kd is not None:assert abs(kd-rec['knockdown'])<1e-10
   mask=measured&(cp>5);mask[own]=False;gain=float(np.clip(.9/max(.35,.9 if kd is None else kd),.8,1.8));eff.append(lfc*gain);var.append(variance*gain*gain);valid.append(mask);n.append(len(ids));kds.append(kd)
  e=np.array(eff);v=np.array(var);mask=np.array(valid);nn=np.array(n,float);initial=nn[:,None]*mask;total=initial.sum(0);w=np.divide(initial,total,out=np.zeros_like(initial),where=total>0);start=(w*e).sum(0);tau=np.maximum((w*(e-start)**2).sum(0)-(w*v).sum(0),0);prec=mask/(v+tau+.0025);floor=np.min(np.where(mask,prec,np.inf),axis=0);floor[~np.isfinite(floor)]=0;prec=np.minimum(prec,10*floor);den=prec.sum(0);wm=np.divide(prec,den,out=np.zeros_like(prec),where=den>0);mu=(wm*e).sum(0);vv=np.divide(1,den,out=np.full_like(den,1e6),where=den>0)+tau;confidence=np.clip(2*special.ndtr(abs(mu)/np.sqrt(np.maximum(vv,.0025)))-1,0,1)*nn.sum()/(nn.sum()+200);mu[den==0]=0;confidence[den==0]=0;ti=int(np.flatnonzero(m.data()['targets']==target)[0]);er=float(abs(mu.astype(np.float32)-t[f'guide_effect_{ci}'][ti]).max());ec=float(abs(confidence.astype(np.float32)-t[f'guide_confidence_{ci}'][ti]).max());assert er<2e-6 and ec<2e-6,(condition,target,er,ec);output.append({'target':target,'guides':len(records),'treated_rows':len(seen),'max_response_error':er,'max_confidence_error':ec})
 return {'condition':condition,'source_checks':output}

def main():
 plan=m.load_plan()
 for name,h in plan['inputs'].items():assert b.sha256(m.ROOT/name)==h
 t,rec=m.load_teacher();assert [r['condition'] for r in rec['metadata']]==m.CONDITIONS;assert not any('D1' in r['condition'] or 'D3' in r['condition'] for r in rec['metadata'])
 source=[source_audit(c,i,rec['metadata'][i],t) for i,c in enumerate(m.CONDITIONS)]
 for row in rec['calibration']:assert row['held_guide'] not in row['fitting_guides'] and row['held_cells']>=30
 alpha_scores={a:float(np.mean([np.mean([r['mse_ratio'] for r in rec['calibration'] if r['alpha']==a and r['condition']==c]) for c in m.CONDITIONS])) for a in [.125,.25,.5]};assert min(alpha_scores,key=lambda a:(alpha_scores[a],a))==rec['chosen_alpha']
 screen=json.loads((m.A/'mean_screen.json').read_text());assert screen['plan_sha256']==b.sha256(m.A/'plan.json');assert screen['teacher_sha256']==rec['sha256'];summary={}
 for arm in m.ARMS:
  vals=[]
  for study in ['H1','Orion','CD4_D1']:
   rows=[r for r in screen['records'] if r['arm']==arm and r['study']==study];assert rows;vals.append({k:float(np.mean([r[k] for r in rows])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']})
  summary[arm]={k:float(np.mean([r[k] for r in vals])) for k in vals[0]}
  for k in summary[arm]:assert abs(summary[arm][k]-screen['summary'][arm][k])<1e-12
 routes=[]
 for ctx in plan['contexts']:
  assert ctx['family'] not in m.CONDITIONS
  if ctx['family'].startswith('CD4'):assert ctx['excluded']==['CD4']
  if ctx['family'] in ['HCT116','HEK293T']:assert set(ctx['excluded'])=={'HCT116','HEK293T'}
  x,meta,measured,src=m.raw(ctx['family'],m.data()['genes']);query=np.asarray(x[ctx['left']].astype(float).multiply((1e6/x[ctx['left']].sum(1).A1)[:,None]).mean(0)).ravel();query[~measured]=np.nan;effect,base,route=m.predict(query,ctx['excluded'],'guide');own=np.array([np.flatnonzero(m.data()['genes']==g)[0] for g in m.data()['targets']]);np.testing.assert_array_equal(effect[np.arange(300),own],base[np.arange(300),own]);routes.append({'family':ctx['family'],**route})
 out={'passed':True,'plan_sha256':b.sha256(m.A/'plan.json'),'teacher_sha256':rec['sha256'],'raw_guide_reconstruction':source,'held_guide_calibration_roles_disjoint':True,'alpha_recomputed':rec['chosen_alpha'],'development_summary_arithmetic_recomputed':True,'development_routes':routes,'auditor_sha256':b.sha256(Path(__file__)),'limits':'Reconstructs selected guide teachers and source roles. Mean-screen arithmetic verified; no claim of independent-study replication, generated-count/Wilcoxon or fresh-donor success.'};m.write(m.A/'independent_audit.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
