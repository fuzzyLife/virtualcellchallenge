"""Record V49 numerical rejection, preserving all completed partial count evidence."""
import v49_balanced_population as m
def main():
 m.check();folder=m.A/'development';records=[];panels=[]
 for f in folder.glob('*.json'):
  rr=m.read(f)
  if f.name.endswith('__panel.json'):panels.append(rr)
  elif 'arrays_sha256' in rr:
   assert m.b.sha256(f.with_suffix('.npz'))==rr['arrays_sha256'];records.append(rr)
 attempts=[{'bins':32,'mean_reference_L1':.0011898812981397961},{'bins':64,'mean_reference_L1':.000886847098715499},{'bins':128,'mean_reference_L1':.0006114200180131397}]
 failure={'error':'Bounded projection misses CPM contract <=0.0005','attempts':attempts,'context_has_no_completed_result':'HCT116','partial_records':len(records),'completed_panels':len(panels),'scope':'Exact integer row/gene totals and bounded rounding tests passed, but full raw-cell validation did not complete. Partial CD4 fidelity still declines; corrected rounding is not sufficient for a better model. No confirmation/reserve/generation.'}
 assert not (m.A/'HCT116_development_results.json').exists()
 m.write(m.A/'development_partial_results.json',{'records':records,'panels':panels,'failure':failure,'complete':False})
 m.write(m.A/'selection.json',{'selected':'v24','eligible':[],'relative_gain':{},'checks':{'full_development_completed':False,'numerical_CPM_contract':False},'failure':failure,'script_sha256':m.b.sha256(m.ROOT/'record_v49_failure.py'),'confirmed_improvement':False})
 print('V49 numerically rejected; partial records retained',len(records),flush=True)
if __name__=='__main__':main()

