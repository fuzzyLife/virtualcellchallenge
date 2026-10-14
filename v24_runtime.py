"""Guard exact numerical regeneration; official packaging uses its separate runtime."""
from pathlib import Path
import json,sys,importlib.metadata

def require_science():
 record=json.loads((Path(__file__).resolve().parent/'experiments/v24_cd4_response/analysis/runtime.json').read_text());assert sys.version.split()[0]==record['python'].split()[0],'Use the pinned V24 scientific Python, not the official packaging runtime'
 for package,version in record['versions'].items():assert importlib.metadata.version(package)==version,'Numerical runtime differs for '+package+': use requirements-v23-science.txt'
 return record
