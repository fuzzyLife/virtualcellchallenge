"""Tests for reference alignment, donor independence and source fallback."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
import build_cd4_model as model
import build_submission as build

def metadata(targets):
 n=len(targets)
 return pd.DataFrame({'target_contrast_gene_name':targets,'target_contrast':targets,
  'culture_condition':['Rest']*n,'n_cells_target':[200]*n,'ontarget_significant':[True]*n,
  'distal_offtarget_flag':[False]*n,'neighboring_gene_KD':[False]*n,
  'single_guide_estimate':[False]*n,'low_target_gex':[False]*n})
class CD4Tests(unittest.TestCase):
 def test_quality_filters_reject_each_known_failure(self):
  m=metadata(['T']*7)
  for i,col in enumerate(['distal_offtarget_flag','neighboring_gene_KD','single_guide_estimate','low_target_gex']):m.loc[i,col]=True
  m.loc[4,'ontarget_significant']=False;m.loc[5,'n_cells_target']=99
  np.testing.assert_array_equal(model.eligible(m),[False]*6+[True])
 def test_alignment_handles_reordered_axes_and_excludes_own_target(self):
  lm=metadata(['T']*20);lm['target_contrast']=[f'e{i}' for i in range(20)]
  rm=lm.iloc[::-1].reset_index(drop=True)
  ld={'genes':np.array(['T','G']),'log_fc':np.tile([9.,.4],(20,1)),'adj_p_value':np.zeros((20,2))}
  rd={'genes':np.array(['G','T']),'log_fc':np.tile([.4,-9.],(20,1)),'adj_p_value':np.zeros((20,2))}
  rows,left,right,_,_,valid=model.align_disjoint((lm,ld),(rm,rd))
  self.assertEqual(len(rows),20);np.testing.assert_allclose(left[:,0],right[:,0]);self.assertFalse(valid[:,1].any())
 def test_hybrid_averages_conditions_and_preserves_unmeasured_fallback(self):
  with tempfile.TemporaryDirectory() as name:
   root=Path(name);directory=root/'reference/cd4_probe';directory.mkdir(parents=True)
   output=root/'analysis/model_audit';output.mkdir(parents=True);(output/'cd4_extraction.json').write_text('[]')
   m=metadata(['T','T']);m.loc[1,'culture_condition']='Stim8hr';m['source_row']=[0,1]
   m.to_csv(directory/'merged.csv',index=False)
   np.savez_compressed(directory/'merged.npz',genes=np.array(['G','T']),source_rows=np.array([0,1]),log_fc=np.array([[1.,-2.],[1.,-2.]]),adj_p_value=np.array([[.001,.001],[1.,.001]]))
   k=(np.array([[[.5,.2],[.1,.8]]]),np.array([[True,True]]),np.array([[2.,3.]]),np.array([0,2]),{'shrinkage':.7},['K562'],{})
   with patch.object(model,'select_settings',return_value={'fdr_threshold':.05,'transfer_shrinkage':.7}),patch.object(build,'reference_model',return_value=k):
    fitted=model.fit_model(root,np.array(['G','T','F','U']),np.array(['T','U']),output,.7)
   np.testing.assert_allclose(fitted[0][0,0],[.35,-1.4,.14,0.])
   np.testing.assert_allclose(fitted[0][0,1],[.07,0.,.56,0.])
   self.assertTrue(fitted[1].all())
 def test_donor_selection_uses_downstream_agreement_not_on_target(self):
  with tempfile.TemporaryDirectory() as name:
   root=Path(name);directory=root/'reference/cd4_probe';directory.mkdir(parents=True)
   output=root/'analysis/model_audit';output.mkdir(parents=True)
   names=['D1_D2','D3_D4']
   (output/'cd4_extraction.json').write_text(json.dumps([{'subset':value+'.npz'} for value in names]))
   targets=[f'T{i}' for i in range(20)];genes=np.array(targets+['G'])
   for side,label in enumerate(names):
    m=metadata(targets);m['source_row']=np.arange(20);m.to_csv(directory/(label+'.csv'),index=False)
    lfc=np.zeros((20,21));lfc[:,-1]=1.;lfc[np.arange(20),np.arange(20)]=9. if side==0 else -9.
    q=np.ones((20,21));q[:,-1]=.001;q[np.arange(20),np.arange(20)]=.001
    np.savez_compressed(directory/(label+'.npz'),genes=genes,source_rows=np.arange(20),log_fc=lfc,adj_p_value=q)
   selected=model.select_settings(directory,output)
   self.assertEqual(selected['transfer_shrinkage'],1.)
   self.assertAlmostEqual(selected['relative_lfc_mse'],0.)
 def test_missing_basal_avoids_invented_context_correlations(self):
  data=(np.array([[[1.,0.,0.]]]),np.array([[True]]),np.zeros((1,3)),np.arange(3),{'shrinkage':1.,'temperature':0.},['hybrid'],{'basal_profile_available':False})
  factors,diag=build.context_factors(np.array([20.,5.,30.]),data,np.array(['G','T','F']),np.array(['T']),.9)
  np.testing.assert_allclose(factors,[[2.,.1,1.]])
  self.assertIsNone(diag['source_profile_correlations']['hybrid'])
if __name__=='__main__':unittest.main(verbosity=2)