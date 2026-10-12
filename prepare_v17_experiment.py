"""Freeze V17 comparisons before reading their outcomes."""
import json
import numpy as np,pandas as pd
import build_submission as b
from v17_methods import ROOT,A,CANDIDATES,RAW_FAMILIES

def main():
 A.mkdir(parents=True,exist_ok=True);path=A/'plan.json'
 if path.exists():raise FileExistsError('Preserve frozen V17 plan')
 old=json.loads((ROOT/'experiments/v16_control_preserving/analysis/plan.json').read_text());contexts=old['contexts']
 verification={'H1':['MED15','ANKZF1'],'HCT116':['PATL1','RNF141'],'HEK293T':['ADNP','GLI3','CAT','MYRF']}
 for fi,c in enumerate(contexts):
  c['context']=c['family']+'_v17';c['cells']=[r for r in c['cells'] if r['split']=='development'];c['null_cases']=[]
  file=ROOT/('reference/v5_h1' if c['family']=='H1' else 'reference/v7_sources/'+c['family'])/'selected_cells.csv'
  meta=pd.read_csv(file);labels=meta['target_gene' if c['family']=='H1' else 'gene_target'].to_numpy(str)
  for seed in old['seeds']:
   rng=np.random.default_rng(np.random.SeedSequence([17012,fi,seed]));rows=rng.choice(c['left'],600,replace=False);c['null_cases'].append({'seed':seed,'template':rows[:400].tolist(),'actual':rows[400:].tolist()})
   for target in verification[c['family']]:
    actual=np.flatnonzero(labels==target);assert len(actual)>=400
    actual=rng.choice(actual,400,replace=False).tolist();template=rng.choice(c['left'],400,replace=False).tolist()
    c['cells'].append({'target':target,'seed':seed,'actual':actual,'template':template,'reference':c['right'],'split':'confirmation'})
 deps=['v17_methods.py','prepare_v17_experiment.py','run_v17_experiment.py','assess_v17_experiment.py','test_v17_heterogeneity.py','collect_v17_resources.py','v16_methods.py','fit_v7_response.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py']
 record={'version':17,'contexts':contexts,'seeds':old['seeds'],'candidates':CANDIDATES,'eligible_finalists':['learned_half','learned_full'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(ROOT/'experiments/v16_control_preserving/analysis/plan.json'),'public_metadata_sha256':b.sha256(A/'public_data_audit.json'),
 'mean_response':'Exact V16/V13 response, full held family excluded from response and variation fitting. Production effects copied bitwise from V16. Library sizes, Poisson weight .65, residual weight .10 and own KD .10 fixed.',
 'variation':'Source raw H1/HCT116/HEK293T; batch-matched controls, two-fold cross-fitting of response direction vs projection variance. Exclude all 300 panel genes. >=80 cells/target, >=8 controls/batch/half; subtract batch-matched control means, pool treated halves across batches; subtract control-center estimation noise; top256 standardized response genes above .1, >=10 genes. Excess treatment/control projected variance divided by squared mean shift; cap CV2 .5, shrink n/(n+100)*SNR/(SNR+10), abstain SNR<4. Equal batch/fold then source family averages; absent targets CV2=0.',
 'release':'V16 paired comparator. Old available raw guard thresholds, >=3 improved metrics including DE; no zero-imputation of missing NMAE. Learned finalist must beat generic and shuffled controls on >=2 fidelity/reach/Jaccard measures and reduce response-projection Wasserstein distance. Same criterion on once-only verification, no retuning. Null emitter must be bitwise identical to V16.',
 'selection':'Among passing learned_half/full, equal improvement-rank average across six raw metrics; lower Wasserstein tie-break. Generic and target-permuted CV2 are negative controls, never finalists.',
 'verification':'Eight algorithmically reserved target/family pairs, two seeds each. Identities selected from cell counts only. Public source biology previously available and partly inspected: this is algorithm verification, not blind biological validation. Old V16 verification names spent and not reused here.',
 'interpretation':'Shared latent strength models joint response variation, not measured causal cell-state conditioning. Training original raw selected counts; new essential screen metadata has zero target overlap and is not trained as if it supplied responses. No leaderboard/private anchor fitting.'}
 b.write_json(path,record);print('V17 plan frozen',b.sha256(path),flush=True)
if __name__=='__main__':main()
