"""Run frozen V58 counts with a sparse-row-count compatibility fix in logging.

The frozen module's only len(ref) call is its null-summary metadata. Its numerical
functions remain byte-identical. Existing completed comparisons are never refit.
"""
import argparse,builtins
from functools import lru_cache
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_v57_pie_evidence as evidence
from run_v58_cached_stage import ReadOnlyNumpy
import evaluate_v58_pie_counts as frozen
from v58_pie import *

def row_count(value):
 return value.shape[0] if sparse.issparse(value) else builtins.len(value)

def main(source):
 assert (ROOT/'evaluate_v58_pie_counts.py').read_text().count('len(ref)')==1
 assert row_count(sparse.csr_matrix((17,5)))==17 and row_count([1,2])==2
 evidence.np=ReadOnlyNumpy();evidence.read=lru_cache(maxsize=None)(evidence.read)
 frozen.len=row_count
 with threadpool_limits(limits=2):frozen.run(source)
 write(A/'null_logging_compatibility.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'frozen_evaluator_sha256':sha(ROOT/'evaluate_v58_pie_counts.py'),'wrapper_sha256':sha(Path(__file__)),'scope':'Module-local len on sparse matrices returns row count; ordinary len unchanged. Only frozen len(ref) call writes null-summary reference_cells. No array, count formula, seed, model or comparison changed.'})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',choices=SOURCES,required=True);main(p.parse_args().source)
