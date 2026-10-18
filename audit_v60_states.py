"""Reuse independent V59 arithmetic/boundary auditor on separately activated V60."""
import argparse,shutil
import v60_depth_adjusted_states as v

def main(stage):
 v.activate();m=v.m;p=m.plan()
 for name,key in [('pre_score_correction.json','parent_CPM_correction_sha256'),('pre_transfer_mask_correction.json','parent_mask_correction_sha256')]:
  source=v.PARENT_A/name;assert m.sha(source)==p[key];dest=v.A/name
  if dest.exists():assert m.sha(dest)==m.sha(source)
  else:shutil.copyfile(source,dest)
 import audit_v59_state_composition as audit
 audit.activate=lambda:None
 if stage=='boundary':audit.boundary()
 else:audit.final()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['boundary','final']);main(p.parse_args().stage)
