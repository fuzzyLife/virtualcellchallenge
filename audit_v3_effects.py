"""Compare modeled V2/V3 effects and report mixture gene-resampling stability."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from build_submission import write_json

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v3_context_ensemble/analysis'


def main():
    with np.load(OUT/'reference_model_arrays.npz',allow_pickle=False) as z:
        effects=z['effects'].astype(float); genes=z['genes']; targets=z['targets']
    with np.load(ROOT/'experiments/v2_cd4_hybrid/analysis/reference_model_arrays.npz',allow_pickle=False) as z:
        old=z['effects'][0].astype(float)
    records=[]; summary=[]
    for ci,context in enumerate(['A','B','C']):
        for ti,target in enumerate(targets):
            mask=genes!=target
            a=effects[ci,ti,mask]; b=old[ti,mask]
            norm=np.linalg.norm(a)*np.linalg.norm(b)
            records.append({'context':context,'target_gene':target,'v2_effect_rms':float(np.sqrt(np.mean(b*b))),
                'v3_effect_rms':float(np.sqrt(np.mean(a*a))),'v2_v3_cosine':float(a@b/norm) if norm>0 else None,
                'v2_downstream_above_abs_lfc_0p1':int((np.abs(b)>.1).sum()),'v3_downstream_above_abs_lfc_0p1':int((np.abs(a)>.1).sum())})
        summary.append({'context':context,'targets_with_nonzero_downstream_estimate':int(np.any(effects[ci]!=0,axis=1).sum()),
            'downstream_pairs_above_abs_lfc_0p01':int((np.abs(effects[ci])>.01).sum()),
            'downstream_pairs_above_abs_lfc_0p1':int((np.abs(effects[ci])>.1).sum()),
            'mean_target_effect_rms':float(np.sqrt(np.mean(effects[ci]*effects[ci],axis=1)).mean())})
    pd.DataFrame(records).to_csv(OUT/'v2_v3_effect_diagnostics.csv',index=False)
    write_json(OUT/'effect_diagnostics.json',{'contexts':summary,'scope':'Modeled effects only; more changes are not evidence of accuracy'})
    boot=pd.read_csv(OUT/'mixture_gene_bootstraps.csv').query("basis=='cell_profiles'")
    boot.groupby('context')[['K562','HCT116','HEK293T','H1','CD4']].quantile([.025,.5,.975]).to_csv(OUT/'mixture_bootstrap_quantiles.csv')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
