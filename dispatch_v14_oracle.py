"""Run unchanged oracle evaluation as each frozen source model becomes available.

Only schedules cache reads; numerical fitting/emission/scoring are unchanged.
"""
from pathlib import Path
import time
from threadpoolctl import threadpool_limits
import run_v14_experiment as run
import build_submission as build
original=run.model_cache


def wait_for_model(data,family):
    side=run.A/'models'/f'{family}.json'
    while not side.exists():time.sleep(5)
    return original(data,family)


if __name__=='__main__':
    record={'script_sha256':build.sha256(Path(__file__)),'evaluator_sha256':build.sha256(run.ROOT/'run_v14_experiment.py'),
            'scope':'Scheduling only: wait for existing fitted model sidecars before calling unchanged model_cache. Does not refit in parallel or change any numerical method.'}
    path=run.A/'oracle_dispatch.json'
    if path.exists():
        import json
        assert json.loads(path.read_text())==record
    else:build.write_json(path,record)
    run.model_cache=wait_for_model
    with threadpool_limits(limits=2):run.evaluate('oracle',3)
