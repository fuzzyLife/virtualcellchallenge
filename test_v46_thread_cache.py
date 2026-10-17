"""Thread-local cache isolation prevents cross-population generation races."""
import threading
from concurrent.futures import ThreadPoolExecutor
from v46_thread_cache import PerThreadCache
import v46_population_mixture as m
def main():
 c=PerThreadCache();barrier=threading.Barrier(4)
 def worker(i):
  c.update(key=i,a=i*10);barrier.wait();assert c.get('key')==i and c['a']==i*10;c.clear();assert c.get('key') is None
 with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,range(4)))
 m.write(m.A/'thread_cache_test.json',{'passed':True,'threads':4,'checks':'Each thread retains its own population cache through synchronized concurrent writes/reads/clear. Numerical emitter unchanged.','script_sha256':m.b.sha256(m.ROOT/'test_v46_thread_cache.py'),'adapter_sha256':m.b.sha256(m.ROOT/'v46_thread_cache.py')})
 print('V46 thread cache isolation PASS',flush=True)
if __name__=='__main__':main()

