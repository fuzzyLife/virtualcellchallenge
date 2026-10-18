"""V57 common paths and scientific contracts; no mutation of prior experiments."""
from pathlib import Path
import hashlib,json,tempfile,os,math
from functools import lru_cache
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v57_pie_cross/analysis'
CACHE=Path(tempfile.gettempdir())/'vcc_v57_pie'
R=ROOT/'reference/v18_pie'
SOURCES=('H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2')
MATCH={'H1':('arc_vcc_25',None),'HCT116':('orion','hct116'),'HEK293T':('orion','hek293t')}
SEED=57057

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(16*1024**2),b''):h.update(x)
 return h.hexdigest()

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if hasattr(x,'tolist'):return clean(x.tolist())
 if isinstance(x,float) and not math.isfinite(x):return None
 return x
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(x),indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
def split(t):return int(hashlib.sha256(('v57-calibration:'+str(t)).encode()).hexdigest()[:8],16)%5

@lru_cache(1)
def plan():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,'Frozen V57 code changed: '+n
 return p
