"""Offline smoke test of CLI 0.3 final-manifest metadata; never upload toy predictions."""
from pathlib import Path
import tempfile,json,subprocess,sys,importlib.metadata
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from vcc.vccfile import read_vcc_meta,validate_vcc
from vcc.panel import perts_fingerprint
import build_submission as b
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def main():
 with tempfile.TemporaryDirectory(prefix='final_v91_cli_smoke_') as tmp:
  p=Path(tmp);pd.DataFrame({'gene_name':['ACLY','VIM','GAPDH']}).to_csv(p/'gene_names.csv',index=False);pd.DataFrame({'target_gene':['ACLY']}).to_csv(p/'pert_counts.csv',index=False);(p/'manifest.json').write_text(json.dumps(dict(contexts=['D','E','F'],partition='test',panel_id='local-final-smoke')),encoding='utf-8');obj=ad.AnnData(sparse.csr_matrix(np.array([[10,100,900],[20,200,800],[30,300,700]],np.int32)),obs=pd.DataFrame({'context':pd.Categorical(['D','E','F']),'target_gene':pd.Categorical(['ACLY']*3)},index=['smoke_D','smoke_E','smoke_F']),var=pd.DataFrame(index=['ACLY','VIM','GAPDH']));obj.write_h5ad(p/'smoke.h5ad');cmd=[sys.executable,'-m','vcc','--json','prep',str(p/'smoke.h5ad'),'-g',str(p/'gene_names.csv'),'--perts',str(p/'pert_counts.csv'),'--cells-per-pert','1','-o',str(p/'smoke.vcc')];run=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace');assert run.returncode==0,run.stderr;validate_vcc(p/'smoke.vcc');meta=read_vcc_meta(p/'smoke.vcc');assert meta['contexts']==['D','E','F'] and meta['panel_id']=='local-final-smoke';assert meta['perts_fingerprint']==perts_fingerprint(str(p/'pert_counts.csv'));b.write_json(A/'cli_final_smoke.json',dict(passed=True,cli_version=importlib.metadata.version('vcc-cli'),automatic_manifest_context_inference=True,metadata=meta,uploaded=False,production_prediction=False,scope='Offline metadata round-trip only; toy file is deleted and must never be submitted.',script_sha256=b.sha256(Path(__file__))));print('Final CLI offline metadata smoke passed',flush=True)
if __name__=='__main__':main()
