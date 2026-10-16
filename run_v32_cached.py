"""Run frozen V32 using byte-verified baseline caches, preserving arithmetic."""
from threadpoolctl import threadpool_limits
import v32_methods as m
import v32_runtime as runtime
from run_v32_experiment import main
if __name__=='__main__':
 m.write(m.A/'runtime_cache_amendment.json',{'parent_plan_sha256':m.sha(m.A/'plan.json'),'wrapper_sha256':m.sha(m.ROOT/'run_v32_cached.py'),'cache_code_sha256':m.sha(m.ROOT/'v32_runtime.py'),'reason':'Original V24 teacher loader rehashes multi-gigabyte original inputs on every baseline call. Cache exact verified existing teacher arrays and deterministic baseline outputs; no changes to model, split, scores or thresholds. Interrupted before any inner scores/outer results were written. Earlier partial source fit may be deterministically repeated.','scope':'Computational optimization only; frozen scientific code unchanged.'})
 runtime.install()
 with threadpool_limits(limits=2):main()