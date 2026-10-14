"""Prove the optimized Stack mean readout matches the original author get_prediction."""
from pathlib import Path
import os,sys,json,time,hashlib,warnings
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';A=ROOT/'experiments/v20_stack_transfer/analysis';os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,pandas as pd,torch,anndata as ad
from scipy import sparse
from stack.models.finetune import ICL_FinetunedModel
from stack_v20_mean_adapter import predict_means,mixed_windows
from smoke_v20_stack import sha

def main():
 pd.options.future.infer_string=False;torch.set_num_threads(4);torch.manual_seed(20020);torch.use_deterministic_algorithms(True);assert sha(R/'bc_large_aligned.ckpt')=='f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36';ck=torch.load(R/'bc_large_aligned.ckpt',map_location='cpu',weights_only=False);model=ICL_FinetunedModel(**ck['hyper_parameters']['model_config']);model.load_state_dict({k[6:]:v for k,v in ck['state_dict'].items()},strict=True);del ck;model.to('cuda:0').eval();folder=A/'smoke_inputs';base=sparse.load_npz(folder/'source_control.npz');query=sparse.load_npz(folder/'query_control.npz');genes=json.loads((R/'genes.json').read_text());t=time.time();mu,theta=predict_means(model,base,query);elapsed=time.time()-t;repeat,repeat_theta=predict_means(model,base,query);np.testing.assert_array_equal(mu,repeat);np.testing.assert_array_equal(theta,repeat_theta);windows,keep=mixed_windows(base,query);inputs=ad.AnnData(sparse.csr_matrix(windows.reshape(-1,15012)),var=pd.DataFrame(index=genes));torch.use_deterministic_algorithms(True,warn_only=True)
 original=model.get_prediction(inputs,str(R/'genes.pkl'),cell_ratio=.25,context_ratio=.4,batch_size=4,num_workers=0,show_progress=False,filter_organism=False,random_seed=20020)
 np.testing.assert_array_equal(mu,original[0][keep]);np.testing.assert_array_equal(theta,original[1][keep]);assert np.isfinite(mu).all() and (mu>=0).all() and np.isfinite(theta).all() and (theta>0).all();result={'passed':True,'mean_bitwise_matches_author':True,'dispersion_bitwise_matches_author':True,'mean_and_dispersion_bitwise_repeat':True,'shape':list(mu.shape),'seconds_optimized_forward':elapsed,'peak_gpu_bytes':torch.cuda.max_memory_allocated(),'optimized_path_deterministic_algorithms':True,'sampler_omitted':'Original author mean/dispersion outputs retained; no stochastic NB sample or median is needed for the mean expert. Final submission counts still come from validated integer emitter.','adapter_sha256':sha(ROOT/'stack_v20_mean_adapter.py'),'test_sha256':sha(Path(__file__))};(A/'mean_adapter_parity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
