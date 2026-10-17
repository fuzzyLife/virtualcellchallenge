"""Reuse mathematically identical zero-effect tests and V24 comparator counts.

Only file/receipt namespaces change. Model-dependent treatment populations
are never copied. Exact same source splits, seeds and emitter implementation.
"""
import argparse,shutil
from pathlib import Path
import v38_response_hurdle as m
import v39_kolf_hurdle as extra
def main():
 extra.install();from run_v38_response_hurdle import check
 check();selected=m.read(extra.A/'mean_results.json')['selected'];parent=m.ROOT/'experiments/v38_response_hurdle/analysis'
 out=extra.A/'null';out.mkdir(exist_ok=True);r=m.read(parent/'null_results.json');records=[]
 for rr in r['records']:
  oldarm=rr['arm'];newarm='v24' if oldarm=='v24' else selected+('_transport' if oldarm.endswith('_transport') else '')
  stem=f"{rr['family']}__-1__0__";source=parent/'null'/(stem+oldarm+'.npz');assert m.b.sha256(source)==rr['arrays_sha256'];destination=out/(stem+newarm+'.npz');shutil.copyfile(source,destination);assert m.b.sha256(destination)==rr['arrays_sha256'];new={**rr,'arm':newarm};m.write(out/(stem+newarm+'.json'),new);records.append(new)
 m.write(extra.A/'null_results.json',{'records':records,'panels':[],'arms':['v24',selected,selected+'_transport'],'scope':'Exact zero-effect reuse from V38: null runner does not call any response model; same zero vector, seeds, control splits, template, means and emitter implementation.','parent_null_results_sha256':m.b.sha256(parent/'null_results.json')})
 count=0;dest=extra.A/'development';dest.mkdir(exist_ok=True)
 for source in sorted((parent/'development').glob('*__v24.json')):
  rr=m.read(source);arr=source.with_suffix('.npz');assert m.b.sha256(arr)==rr['arrays_sha256'];shutil.copyfile(arr,dest/arr.name);shutil.copyfile(source,dest/source.name);count+=1
 m.write(extra.A/'identical_comparator_reuse.json',{'adapter_sha256':m.b.sha256(m.ROOT/'reuse_v39_comparators.py'),'parent_V38_plan_sha256':m.b.sha256(parent/'plan.json'),'new_V39_plan_sha256':m.b.sha256(extra.A/'plan.json'),'null_cases_copied':len(records),'V24_treatment_cases_copied':count,'response_model_dependent_treatments_copied':0,'proof':'V38/V39 plans pin identical v38_response_hurdle.py, v16_methods.py, run_v38_response_hurdle.py and original V34 cases. V24 comparator effects come from the same SHA-checked V37 cache; null effects are exactly zero. Reused arrays checked against original receipts.'});print('Reused',len(records),'null and',count,'identical V24 comparator cases',flush=True)
if __name__=='__main__':main()

