"""Run the unchanged V69 evaluator with content-addressed, bounded memoization.

Only repeated DE/jackknife computations are reused. The complete CSR contents,
shape and dtype enter each key; cached arrays are read-only. No fitted model,
sample, statistic, candidate rule or original script is changed.
"""
import hashlib
from collections import OrderedDict
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import run_v69_experiment as runner
from v57_pie import sha,write
DE=runner.de;PROFILE=runner.profile_jackknife_chunked
CACHE={};HITS={'de':0,'profile':0};MISSES={'de':0,'profile':0}


def signature(x):
    x=sparse.csr_matrix(x);h=hashlib.sha256(str((x.shape,x.dtype.str)).encode())
    for a in (x.data,x.indices,x.indptr):
        h.update(a.dtype.str.encode());h.update(memoryview(np.ascontiguousarray(a)))
    return h.digest()


def immutable(v):
    if isinstance(v,np.ndarray):v.setflags(write=False)
    elif isinstance(v,dict):
        for z in v.values():immutable(z)
    elif isinstance(v,tuple):
        for z in v:immutable(z)
    return v


def memo(kind,key,compute):
    cache=CACHE.setdefault(kind,OrderedDict())
    if key in cache:
        HITS[kind]+=1;v=cache.pop(key);cache[key]=v;return v
    MISSES[kind]+=1;v=immutable(compute());cache[key]=v
    if len(cache)>16:cache.popitem(last=False)
    return v


def de(x,control,genes):
    gene_digest=hashlib.sha256('\0'.join(map(str,genes)).encode()).digest()
    return memo('de',(signature(x),signature(control),gene_digest),lambda:DE(x,control,genes))


def profile(x):return memo('profile',signature(x),lambda:PROFILE(x))


def test():
    rng=np.random.default_rng(691);a=sparse.csr_matrix(rng.poisson(3,(70,24)));b=sparse.csr_matrix(rng.poisson(3,(110,24)));genes=np.array([f'gene{i}' for i in range(24)])
    original=DE(a,b,genes);cached=de(a,b,genes);again=de(a.copy(),b.copy(),genes)
    for k in original:np.testing.assert_array_equal(original[k],cached[k]);np.testing.assert_array_equal(original[k],again[k])
    for x,y in zip(PROFILE(a),profile(a)):np.testing.assert_array_equal(x,y)
    old=signature(a);a.data[0]+=1;assert signature(a)!=old
    changed=de(a,b,genes);fresh=DE(a,b,genes)
    for k in fresh:np.testing.assert_array_equal(changed[k],fresh[k])
    print('Content-cache equivalence and mutation-key checks passed',flush=True)


if __name__=='__main__':
    import sys
    with threadpool_limits(limits=2):
        if sys.argv[1:]==['--self-test']:test()
        else:
            runner.de=de;runner.profile_jackknife_chunked=profile
            runner.main()
            stage=sys.argv[1];source=sys.argv[sys.argv.index('--source')+1] if '--source' in sys.argv else 'all'
            write(runner.m.A/f'cache_receipt_{stage}_{source}.json',dict(wrapper_sha256=sha(Path(__file__)),frozen_runner_sha256=sha(Path(runner.__file__)),hits=HITS,misses=MISSES,cache_key='Complete CSR data/indices/pointers, dtype/shape and gene axis; no identity-only cache',statistic_unchanged=True))
