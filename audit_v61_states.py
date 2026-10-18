"""Independent boundary/array audit for V61, preserving inherited correction receipts."""
import argparse,shutil
import v61_pooled_state_response as v

def main(stage):
 v.activate();m=v.m
 for name in ['pre_score_correction.json','pre_transfer_mask_correction.json']:
  source=v.parent.PARENT_A/name;dest=v.A/name
  if dest.exists():assert m.sha(dest)==m.sha(source)
  else:shutil.copyfile(source,dest)
 import audit_v59_state_composition as audit
 audit.activate=lambda:None
 if stage=='boundary':audit.boundary()
 else:audit.final()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['boundary','final']);main(p.parse_args().stage)
