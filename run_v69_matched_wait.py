"""Overlap fixed matched controls with remaining context jobs; wait for inputs.

Only file readiness is changed. The frozen matched-control evaluator and all
count formulas/seeds remain unchanged.
"""
import sys,time
from pathlib import Path
from threadpoolctl import threadpool_limits
import run_v69_matched_control as helper
from v57_pie import read as original_read,write,sha
READS=[]


def ready_read(path):
    path=Path(path);start=time.monotonic();notified=False
    while not path.exists():
        if not notified:print('Waiting for completed comparison input',path.name,flush=True);notified=True
        if time.monotonic()-start>1800:raise TimeoutError(str(path))
        time.sleep(2)
    value=original_read(path)
    if notified:READS.append(dict(file=path.name,seconds=time.monotonic()-start))
    return value


if __name__=='__main__':
    stage=sys.argv[1];assert stage=='confirmation'
    helper.read=ready_read
    with threadpool_limits(limits=2):helper.run(stage)
    write(helper.m.A/'matched_control_wait_receipt.json',dict(wrapper_sha256=sha(Path(__file__)),fixed_helper_sha256=sha(Path(helper.__file__)),waits=READS,scientific_computation_unchanged=True))
