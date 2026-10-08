"""Generate with frozen V11 numerical code, then validate before safe packaging."""
import argparse
import build_submission as build
from build_v11_submission import generate,load,OUT,A

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 if not 1<=args.workers<=8:p.error('workers must be 1..8')
 generate(False,args.workers)
 manifest,genes,targets,model,effect=load()
 checked=build.validate_h5ad(OUT/'prediction.h5ad',genes,targets,manifest['contexts'],manifest['cells_per_pert'])
 build.write_json(A/'prediction_validation.json',checked)
 print('Frozen V11 count matrix generated and fully validated',flush=True)
