"""Prospective control-only false-DE guard for V81 (only after development selects)."""
import json,gc
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v81_residual_temperature as m
from v57_pie import read,write,sha
from v25_methods import raw
from v6_methods import moments
from v16_methods import emit
from v12_cpu_metrics import de
from audit_v76_counts import spot_ranks

def main():
 m.check();selection=read(m.A/'selection.json');selected=selection['selected'];assert selected,'No qualifying development selection'
 protocol=read(m.A/'null_protocol.json');genes=m.data()['genes'];atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();records=[]
 for ctx in m.contexts():
  name=ctx['family'];out=m.A/'nulls';out.mkdir(exist_ok=True);p=out/(name+'.json')
  if p.exists():records.append(read(p));continue
  x,meta,measured,src=raw(name,genes);spec=ctx['null'];ids=np.asarray(spec['fit_controls']);assert (meta.gene_target.to_numpy(str)[ids]=='Non-Targeting').all();libs,cp,bulk,_=moments(x[ids]);template=x[spec['template']];depth=np.asarray(template.sum(1)).ravel().astype(np.int64);ref=x[spec['reference']];actual=x[spec['actual']]
  parent,_=emit(template,np.arange(len(depth)),np.zeros(len(genes)),None,spec['seed']+81,'residual10',(depth,cp,bulk),atlas);count,diag=m.transform(parent,None,None,selected,spec['seed']+8100000)
  dd=de(count,ref,genes);bd=de(parent,ref,genes);td=de(actual,ref,genes);checks=spot_ranks(count,ref,dd);bc=int((bd['q'][bd['eligible']]<.05).sum());cc=int((dd['q'][dd['eligible']]<.05).sum());tc=int((td['q'][td['eligible']]<.05).sum())
  path=p.with_suffix('.npz');np.savez_compressed(path,pred_p=dd['p'],pred_q=dd['q'],pred_eligible=dd['eligible'],baseline_p=bd['p'],baseline_q=bd['q'],baseline_eligible=bd['eligible'],truth_p=td['p'],truth_q=td['q'],truth_eligible=td['eligible'])
  r=dict(context=name,selected=selected,baseline_false_calls=bc,candidate_false_calls=cc,real_control_false_calls=tc,allowed_false_calls=1.1*bc+10,passes=cc<=1.1*bc+10,diag=diag,manual_rank_coordinates=checks,arrays_sha256=sha(path),scope='Single reused NTC-vs-independent-NTC case; not biological replication or formal FDR calibration.');write(p,r);records.append(r);print(json.dumps(r),flush=True);del x,parent,count;gc.collect()
 write(m.A/'null_assessment.json',dict(selected=selected,passed=all(r['passes'] for r in records),records=records,protocol_sha256=sha(m.A/'null_protocol.json')))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
