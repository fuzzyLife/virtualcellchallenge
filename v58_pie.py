"""V58: V57 with CD4 excluded from every source-calibration baseline."""
from v57_pie import *
PARENT_A=A
PARENT_CACHE=CACHE
A=ROOT/'experiments/v58_pie_holdstudy_guard/analysis'
CACHE=Path(tempfile.gettempdir())/'vcc_v58_pie'

@lru_cache(1)
def plan():
 p=read(A/'plan.json')
 assert sha(PARENT_A/'plan.json')==p['parent_plan_sha256']
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,'Frozen V58 dependency changed: '+n
 return p
