"""Run unchanged frozen V22 scoring, waiting only for complete fitted model files."""
import sys,time,json
from pathlib import Path
import pandas as pd
pd.options.future.infer_string=False
import run_v22_experiment as runner
from v22_methods import A
original=runner.load_model

def wait_model(family):
 start=time.monotonic();last=0
 while not ((A/'models'/f'{family}.npz').exists() and (A/'models'/f'{family}.json').exists()):
  if time.monotonic()-start>3600:raise TimeoutError('Model not ready: '+family)
  if time.monotonic()-last>30:print('Waiting for frozen V22 fitted model',family,flush=True);last=time.monotonic()
  time.sleep(5)
 return original(family)
runner.load_model=wait_model
if __name__=='__main__':runner.main()
