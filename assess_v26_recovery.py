"""Operational Path coercion for the frozen V26 assessor's hash reporting only."""
from pathlib import Path
import argparse,json
import assess_v26_experiment as frozen
from v26_methods import ROOT,A,load_plan,write
import build_submission as b

def main(stage,family):
 plan=load_plan();original=b.sha256;assert original(ROOT/'assess_v26_experiment.py')==plan['dependencies']['assess_v26_experiment.py'];tag='null_'+family if stage=='null' else stage;assert not (A/f'{tag}_result.json').exists();before={p.name:original(p) for p in A.glob(tag+'*.csv')};b.sha256=lambda path:original(Path(path));frozen.main(stage,family)
 for name,h in before.items():assert original(A/name)==h
 write(A/'assessment_recovery.json',{'operation':'Coerce path arguments to Path for the unchanged sha256 helper while running frozen assessor','original_error':"AttributeError: str has no attribute open when hashing __file__ after checks completed; no decision file was written",'frozen_assessor_sha256':plan['dependencies']['assess_v26_experiment.py'],'adapter_sha256':original(Path(__file__)),'stage':stage,'family':family,'result_sha256':original(A/f'{tag}_result.json'),'numeric_model_counts_samples_thresholds_unchanged':True,'input_tables_sha256':before});print('V26 assessment reporting recovered; frozen numerics unchanged',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','development']);p.add_argument('--family');args=p.parse_args();main(args.stage,args.family)
