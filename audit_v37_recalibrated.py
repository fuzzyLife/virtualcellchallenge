"""Verify effect aliases and baseline parity for the separately calibrated mixture."""
import argparse
import numpy as np
import v35_methods as m
from run_v37_ensemble import A as PARENT,read,write,contexts
A=PARENT/'recalibrated';NAME='recalibrated_population50'

def aliases():
 rows=[]
 for ctx in contexts():
  name=ctx['family'];source=PARENT/'cache'/name/'effects.npz';alias=A/'effects'/name/'effects.npz'
  with np.load(source,allow_pickle=False) as a,np.load(alias,allow_pickle=False) as b:
   for key in ('v24','retrieval','cp','bulk'):np.testing.assert_array_equal(a[key],b[key])
   np.testing.assert_array_equal(a['mean50'],b[NAME])
  rows.append({'family':name,'source_sha256':m.sha(source),'alias_sha256':m.sha(alias)})
 return rows

def main(seal=False):
 rows=aliases();manifest=A/'effect_alias_manifest.json'
 if seal:
  assert not manifest.exists();write(manifest,{'rows':rows,'effect_aliases_exact':True,'scope':'Effect alias equality sealed before this corrected arm treatment scoring. Initial null outcomes and prior primary outcomes known.'});print('Effect aliases verified and sealed',flush=True);return
 assert read(manifest)['rows']==rows
 import audit_v37_ensemble as audit
 from run_v37_recalibrated import check
 audit.A=A;audit.check=check;audit.main()
 parity=0
 for stage in ('null','development'):
  if not (A/(stage+'_results.json')).exists():continue
  result=read(A/(stage+'_results.json'));parent=read(PARENT/(stage+'_results.json'));prior={(r['family'],r['case_index'],r['seed_offset']):r for r in parent['records'] if r['arm']=='v24'}
  for r in result['records']:
   if r['arm']!='v24':continue
   q=prior[r['family'],r['case_index'],r['seed_offset']];here=A/stage/f"{r['family']}__{r['case_index']}__{r['seed_offset']}__v24.npz";before=PARENT/q['arrays_path'] if 'arrays_path' in q else PARENT/stage/f"{r['family']}__{r['case_index']}__{r['seed_offset']}__v24.npz"
   with np.load(here,allow_pickle=False) as z,np.load(before,allow_pickle=False) as zz:
    assert set(z.files)==set(zz.files)
    for k in z.files:np.testing.assert_array_equal(z[k],zz[k])
   parity+=1
 result=read(A/('development_results.json' if (A/'development_results.json').exists() else 'null_results.json'));rr=[r for r in result['records'] if r['arm']==NAME];improve={k:float(np.mean([r[k] for r in rr])) for k in ('precalibration_mean_L1','precalibration_bulk_L1','mean_probability_L1','bulk_probability_L1')}
 auditrec=read(A/'independent_audit.json');auditrec.update({'baseline_comparisons_exactly_match_primary':parity,'effect_alias_manifest_sha256':m.sha(manifest),'post_calibration_diagnostics':improve,'adapter_auditor_sha256':m.sha(m.ROOT/'audit_v37_recalibrated.py')});write(A/'independent_audit.json',auditrec);print('Recalibrated alias/baseline parity passed',parity,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--seal',action='store_true');p.add_argument('--namespace',choices=['recalibrated','noise_recalibrated'],default='recalibrated');args=p.parse_args();A=PARENT/args.namespace;NAME='recalibrated_noise_mix' if args.namespace=='noise_recalibrated' else 'recalibrated_population50'
 if args.namespace=='noise_recalibrated':
  import run_v37_noise_recalibrated
 main(args.seal)
