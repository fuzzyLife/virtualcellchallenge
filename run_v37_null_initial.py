"""Sampling-only correction to V37 nulls: use the original frozen three-way split."""
import gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v35_methods as m
from v37_ensemble import emit,ARMS
from run_v37_ensemble import A,check,read,write,contexts
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked

def main():
 check();initial=A/'null_primary_results.json';active=A/'null_results.json'
 if not initial.exists():initial.write_bytes(active.read_bytes())
 dest=A/'null_corrected';dest.mkdir(exist_ok=True);plan=dest/'amendment.json'
 if not plan.exists():
  write(plan,{'parent_plan_sha256':m.sha(A/'plan.json'),'initial_null_results_sha256':m.sha(initial),'adapter_sha256':m.sha(m.ROOT/'run_v37_null_corrected.py'),'known_outcomes':'All44 raw mean errors and initial null outcomes known; development treatment run may be in progress. No treatment result was used to set this correction.','reason':'Original V37 null used first400 left controls and whole-left fitting moments. Reuse the already frozen V34 three-way fit/template/reference split, rather than a potentially unrepresentative prefix.','changes':'Only null row identities, corresponding fitting moments, null seed and reference; candidates, coefficients, thresholds, treatment cases and selection rule unchanged.','deduplication':'Zero-response mean arms share exactly the same emitter and seed; zero-response population mixture chooses identical parent counts. Compute this common population once, save all arm receipts with identical numeric arrays. Noise mix is generated separately.','original_results_preserved':True})
 assert read(plan)['adapter_sha256']==m.sha(m.ROOT/'run_v37_null_corrected.py')
 data=m.old.old.read_cache();genes=data['genes'];zero=np.zeros(len(genes),np.float32);from evaluate_atlas_transfer import load_atlas_model
 atlas,_=load_atlas_model();records=[]
 for ctx in contexts():
  name=ctx['family'];null=ctx['null'];left=set(null['fit_controls']);right=set(null['reference']);templateids=set(null['template']);assert not left&right and not left&templateids and not right&templateids
  completed=[dest/f'{name}__-1__0__{arm}.json' for arm in ('v24',)+ARMS]
  if all(p.exists() for p in completed):
   records.extend(read(p) for p in completed);continue
  x,meta,measured,source=m.raw(name,genes);template=x[null['template']];actual=x[null['actual']];reference=x[null['reference']];_,cp,bulk,_=moments(x[null['fit_controls']]);libs=template.sum(1).A1.astype(np.int64);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual);cc,cj=profile_jackknife_chunked(reference)
  write(dest/(name+'_source.json'),{'source':source,'split':null,'three_way_control_sets_disjoint':True,'fitting_controls':len(left)})
  common=None
  for arm in ('v24',)+ARMS:
   stem=f'{name}__-1__0__{arm}';arrays=dest/(stem+'.npz')
   if common is None or arm=='mean50_noise_mix':
    count,diag=emit(template,zero,zero,zero,None,null['seed'],(libs,cp,bulk),atlas,arm);np.testing.assert_array_equal(count.sum(1).A1,libs);dd=de(count,reference,genes);profile,jack=profile_jackknife_chunked(count);metrics=direction_metrics(dd,truth,genes,None)
    values={'pred':profile,'jack':jack,'truth':tp,'truth_jack':tj,'control':cc,'control_jack':cj,'pred_q':dd['q'],'pred_lfc':dd['lfc'],'pred_eligible':dd['eligible'],'truth_q':truth['q'],'truth_lfc':truth['lfc'],'truth_eligible':truth['eligible']}
    details={**diag,**metrics,'null_DE_calls':int(np.sum(dd['eligible']&(dd['q']<.05))),'actual_null_DE_calls':int(np.sum(truth['eligible']&(truth['q']<.05)))}
    if common is None:common=(values,details)
    del count;gc.collect()
   else:values,details=common
   np.savez_compressed(arrays,**values);r={'family':name,'study':ctx['excluded'][0] if name in ('H1',) else ('CD4' if name.startswith('CD4') else 'Orion'),'target':'Non-Targeting','case_index':-1,'seed_offset':0,'arm':arm,'cells':len(libs),'row_libraries_exact':True,'arrays_path':str(arrays.relative_to(A)),**details,'arrays_sha256':m.sha(arrays)};write(dest/(stem+'.json'),r);records.append(r);print('V37 corrected null',name,arm,r['null_DE_calls'],flush=True)
  del x;gc.collect()
 result={'panels':[],'records':records,'arms':['v24',*ARMS],'scope':'Corrected original three-way control-null split; original prefix-based null outcomes preserved separately.','amendment_sha256':m.sha(plan),'initial_null_results_sha256':m.sha(initial)};write(A/'null_corrected_results.json',result);write(active,result);print('Corrected nulls complete; selection reads corrected split',flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
