"""Control-only variance diagnostic; no held treatment outcomes consumed."""
import numpy as np
from scipy import sparse
from v6_methods import moments
from run_v35_v36 import contexts
from run_v37_ensemble import checked_folder
import build_submission as b
def main():
 records=[]
 for c in contexts():
  folder,_=checked_folder(c['family']);x=sparse.load_npz(folder/'0_template.npz');libs,p,bulk,var=moments(x);pois=p*np.mean(1/libs);keep=p>5e-6;ratio=var[keep]/np.maximum(pois[keep],1e-20)
  records.append({'context':c['family'],'cells':x.shape[0],'library_quantiles':np.quantile(libs,[0,.5,1]).tolist(),'native_to_Poisson_variance_quantiles':np.quantile(ratio,[.1,.5,.9]).tolist(),'eligible_outputs':int(keep.sum()),'cache_receipt_sha256':b.sha256(folder/'receipt.json')})
 b.write_json(b.ROOT/'experiments/v44_response_dispersion/analysis/control_dispersion.json',{'records':records,'scope':'Controls only, normalized composition variance relative to Poisson sampling at observed depths. Similar medians across contexts; does not justify assigning arbitrary variance coefficients by guessed cell identity.','script_sha256':b.sha256(b.ROOT/'analyze_v44_control_dispersion.py')})
 print('Control dispersion diagnostic saved',flush=True)
if __name__=='__main__':main()

