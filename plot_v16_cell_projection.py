"""Plot observed controls and sampled V16 predictions in the frozen full-gene rank PCA."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import time
import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cluster_full_rank_matrix import rank_unit_rows

ROOT=Path(__file__).resolve().parent
BASE=ROOT/'analysis/full_rank_clustering'
OUT=ROOT/'analysis/v16_cell_projection'
CONTEXTS=['A','B','C']
COLORS={'A':'#4477aa','B':'#ee7733','C':'#228833'}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def h5_strings(node):
    # AnnData supports categorical and nullable-string encodings as well as datasets.
    return np.asarray(read_elem(node),dtype=str)


def project(raw,model):
    return (rank_unit_rows(raw)-model['pca_mean'])@model['pca_components'].T


def read_rows(handle,lo,hi,pointers,n_genes):
    begin,end=int(pointers[lo]),int(pointers[hi])
    return sparse.csr_matrix((handle['X/data'][begin:end],handle['X/indices'][begin:end],
                              pointers[lo:hi+1]-begin),shape=(hi-lo,n_genes))


def plots(cells,variance):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    titles=['Observed controls: 55,200 cells','V16 predictions: 18,000 sampled cells','Overlay: same frozen axes']
    for row,(px,py) in enumerate([('PC1','PC2'),('PC3','PC4')]):
        bounds=[]
        for col,source in enumerate(['Observed control','V16 prediction',None]):
            ax=axes[row,col]
            for context in CONTEXTS:
                for kind in ['Observed control','V16 prediction']:
                    if source is not None and kind!=source:continue
                    part=cells[(cells.context==context)&(cells.source==kind)]
                    predicted=kind=='V16 prediction'
                    ax.scatter(part[px],part[py],s=3 if predicted else 2,
                               alpha=.28 if predicted else (.10 if source is None else .24),
                               c=COLORS[context],marker='x' if predicted else 'o',linewidths=.4 if predicted else 0,
                               rasterized=True,label=context+(' predicted' if predicted else ' control'))
            ax.set_title(titles[col] if row==0 else titles[col].split(':')[0]+' — PCs 3/4')
            ax.set_xlabel(f'{px} ({variance[int(px[2:])-1]:.2%} of control rank variance)')
            ax.set_ylabel(f'{py} ({variance[int(py[2:])-1]:.2%} of control rank variance)')
            ax.grid(alpha=.12)
            if row==0:ax.legend(markerscale=3,fontsize=8,loc='best')
            bounds.append(ax)
        # Include every point; use identical limits for the three source views.
        for axis in ['x','y']:
            limits=[getattr(ax,'get_'+axis+'lim')() for ax in bounds]
            for ax in bounds:getattr(ax,'set_'+axis+'lim')((min(x[0] for x in limits),max(x[1] for x in limits)))
    fig.suptitle('Observed A/B/C controls and V16 predictions\n18,533 genes • average tied raw-count ranks • PCA fitted on controls only',fontsize=15)
    fig.savefig(OUT/'cells_observed_vs_v16.png',dpi=180)
    fig.savefig(OUT/'cells_observed_vs_v16.pdf')
    plt.close(fig)


def interactive(cells,variance):
    import plotly.graph_objects as go
    items=[]
    fig=go.Figure()
    for context in CONTEXTS:
        for source in ['Observed control','V16 prediction']:
            part=cells[(cells.context==context)&(cells.source==source)]
            payload={'context':context,'source':source,'pc':part[['PC1','PC2','PC3','PC4']].round(8).to_numpy().tolist(),
                     'meta':part[['cell_id','target_gene','total_counts','detected_genes']].to_numpy().tolist()}
            items.append(payload)
            fig.add_trace(go.Scattergl(x=part.PC1,y=part.PC2,mode='markers',name=context+' — '+source,
                customdata=payload['meta'],marker=dict(color=COLORS[context],size=4 if source=='V16 prediction' else 3,
                opacity=.5 if source=='V16 prediction' else .2,symbol='x' if source=='V16 prediction' else 'circle'),
                hovertemplate='Cell: %{customdata[0]}<br>Target: %{customdata[1]}<br>Counts: %{customdata[2]}<br>Detected genes: %{customdata[3]}<extra>%{fullData.name}</extra>'))
    fig.update_layout(template='plotly_white',height=760,title='Observed controls and V16 predictions — all targets',
                      xaxis_title=f'PC1 ({variance[0]:.2%})',yaxis_title=f'PC2 ({variance[1]:.2%})',
                      legend=dict(orientation='h',y=-.12),margin=dict(t=70,b=110))
    script='''
const allRows=PAYLOAD;
const targets=TARGETS;
const panel=document.createElement('div');
panel.style.cssText='font-family:Arial,sans-serif;padding:16px;background:#f3f5f7;border-radius:8px;margin:12px';
panel.innerHTML='<b>Target:</b> <select id="target"><option value="all">All 300 targets</option></select> &nbsp; <b>Axes:</b> <select id="plane"><option value="0">PC1 / PC2</option><option value="2">PC3 / PC4</option></select><p>Circles: observed unperturbed controls. Crosses: V16 predictions (20 / 400 cells per target/context). Click legend entries to hide/show groups; drag to zoom. No real treated A/B/C cells are available here.</p>';
const gd=document.getElementById('{plot_id}');gd.parentNode.insertBefore(panel,gd);
const target=document.getElementById('target'),plane=document.getElementById('plane');
for(const name of targets){const o=document.createElement('option');o.value=name;o.textContent=name;target.appendChild(o);}
function refresh(){
 const offset=Number(plane.value),selected=target.value;
 const xs=[],ys=[],meta=[];
 for(const group of allRows){
  const keep=group.meta.map((m,i)=>i).filter(i=>group.source==='Observed control'||selected==='all'||group.meta[i][1]===selected);
  xs.push(keep.map(i=>group.pc[i][offset]));ys.push(keep.map(i=>group.pc[i][offset+1]));meta.push(keep.map(i=>group.meta[i]));
 }
 Plotly.restyle(gd,{x:xs,y:ys,customdata:meta});
 Plotly.relayout(gd,{'title.text':'Observed controls and V16 predictions — '+(selected==='all'?'all targets':selected),
                    'xaxis.title.text':'PC'+(offset+1),'yaxis.title.text':'PC'+(offset+2),'xaxis.autorange':true,'yaxis.autorange':true});
}
target.onchange=refresh;plane.onchange=refresh;
'''.replace('PAYLOAD',json.dumps(items,separators=(',',':'))).replace('TARGETS',json.dumps(sorted(cells.loc[cells.source=='V16 prediction','target_gene'].unique())))
    fig.write_html(OUT/'cells_observed_vs_v16.html',include_plotlyjs=True,full_html=True,
                   post_script=script,config={'responsive':True,'displaylogo':False})


def main(per_group=20,seed=20261012):
    started=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
    if per_group!=20:raise ValueError('This descriptive release fixes 20 cells per target/context.')
    report=json.loads((BASE/'analysis_report.json').read_text())
    for name in ['cluster_models.npz','gene_axis.csv','per_cell_full_rank_clusters.csv','pca_variance.csv']:
        assert sha(BASE/name)==report['outputs'][name]['sha256'],name
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    np.testing.assert_array_equal(pd.read_csv(BASE/'gene_axis.csv').gene_name,genes)
    model=np.load(BASE/'cluster_models.npz')
    controls=pd.read_csv(BASE/'per_cell_full_rank_clusters.csv')
    assert len(controls)==55200 and not controls.cell_id.duplicated().any()
    checks=[]
    # Independently check that the cached observed coordinates use this exact transform.
    with threadpool_limits(limits=6):
        for context in CONTEXTS:
            frame=controls[controls.context==context]
            assert len(frame)==18400
            with h5py.File(ROOT/f'context_{context}.h5ad','r') as handle:
                np.testing.assert_array_equal(h5_strings(handle['var/_index']),genes)
                np.testing.assert_array_equal(h5_strings(handle['obs/_index']),frame.cell_id)
                pointers=handle['X/indptr'][:]
                raw=read_rows(handle,0,8,pointers,len(genes)).toarray()
                scores=project(raw,model)
                np.testing.assert_allclose(scores,frame[[f'PC{i+1}' for i in range(8)]].iloc[:8],rtol=0,atol=2e-6)
                checks.append({'context':context,'cached_control_projection_rows_verified':8,
                               'max_absolute_error':float(abs(scores-frame[[f'PC{i+1}' for i in range(8)]].iloc[:8].to_numpy()).max())})
        controls['source']='Observed control';controls['prediction_row']=-1
        path=ROOT/'experiments/v16_control_preserving/prediction.h5ad'
        fingerprint={'bytes':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns}
        audit=json.loads((path.parent/'analysis/independent_validation.json').read_text())
        assert audit['passed'] and audit['production_responses_recomputed_exactly']
        rng=np.random.default_rng(seed);parts=[]
        with h5py.File(path,'r',rdcc_nbytes=128*1024*1024) as handle:
            np.testing.assert_array_equal(h5_strings(handle['var/_index']),genes)
            context=h5_strings(handle['obs/context']);target=h5_strings(handle['obs/target_gene']);ids=h5_strings(handle['obs/_index'])
            assert len(ids)==360000 and len(set(ids))==len(ids)
            pointers=handle['X/indptr'][:]
            boundaries=np.r_[0,np.flatnonzero((context[1:]!=context[:-1])|(target[1:]!=target[:-1]))+1,len(ids)]
            assert len(boundaries)==901 and np.all(np.diff(boundaries)==400)
            panel=set(pd.read_csv(ROOT/'pert_counts.csv').iloc[:,0].astype(str))
            for n,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
                selected=np.sort(rng.choice(hi-lo,per_group,replace=False))
                matrix=read_rows(handle,int(lo),int(hi),pointers,len(genes))[selected]
                raw=matrix.toarray()
                assert np.isfinite(raw).all() and (raw>=0).all() and np.equal(raw,np.floor(raw)).all()
                scores=project(raw,model)
                if n==0:
                    # Same transform must be independent of batching.
                    np.testing.assert_allclose(scores,np.vstack([project(row[None,:],model) for row in raw]),rtol=0,atol=2e-6)
                rows=lo+selected
                frame=pd.DataFrame({'cell_id':ids[rows],'context':context[rows],'target_gene':target[rows],
                                    'total_counts':raw.sum(1,dtype=np.int64),'detected_genes':(raw>0).sum(1),
                                    'source':'V16 prediction','prediction_row':rows})
                for i in range(8):frame[f'PC{i+1}']=scores[:,i]
                parts.append(frame)
                if (n+1)%100==0:print(f'Projected {n+1}/900 groups ({(n+1)*per_group:,} predicted cells)',flush=True)
            predictions=pd.concat(parts,ignore_index=True)
            counts=predictions.groupby(['context','target_gene']).size()
            assert len(counts)==900 and (counts==per_group).all()
            assert set(predictions.context)==set(CONTEXTS) and set(predictions.target_gene)==panel
        assert fingerprint=={'bytes':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns}
    columns=['cell_id','context','target_gene','source','prediction_row','total_counts','detected_genes']+[f'PC{i+1}' for i in range(8)]
    cells=pd.concat([controls[columns],predictions[columns]],ignore_index=True)
    cells.to_csv(OUT/'cell_coordinates.csv',index=False)
    summary=[]
    for (source,context),frame in cells.groupby(['source','context']):
        pcs=frame[[f'PC{i+1}' for i in range(8)]].to_numpy()
        summary.append({'source':source,'context':context,'cells':len(frame),'mean_total_counts':frame.total_counts.mean(),
                        'mean_detected_genes':frame.detected_genes.mean(),'pc8_total_variance':np.var(pcs,axis=0,ddof=1).sum(),
                        **{f'PC{i+1}_mean':pcs[:,i].mean() for i in range(4)}})
    pd.DataFrame(summary).to_csv(OUT/'projection_summary.csv',index=False)
    variance=pd.read_csv(BASE/'pca_variance.csv').explained_variance_ratio.to_numpy()
    plots(cells,variance);interactive(cells,variance)
    outputs=['cell_coordinates.csv','projection_summary.csv','cells_observed_vs_v16.png','cells_observed_vs_v16.pdf','cells_observed_vs_v16.html']
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),'seed':seed,
            'observed_cells':55200,'sampled_prediction_cells':len(predictions),'sampled_per_target_context':per_group,
            'total_prediction_cells':360000,'genes':18533,'groups':900,'sampling':'Simple random sample without replacement separately within every 400-cell target/context group.',
            'embedding':'Frozen eight-PC PCA of per-cell average tied raw-count ranks, centered and unit-normalized across all 18,533 genes. Fitted only on 55,200 original controls; no prediction refit, no pseudocounts.',
            'control_projection_checks':checks,'prediction_projection_batch_check_passed':True,
            'prediction_file':str(path.relative_to(ROOT)),'prediction_file_fingerprint':fingerprint,
            'expected_prediction_sha256_from_prior_full_audit':audit['prediction_sha256'],
            'prediction_integrity_scope':'Axes, all group sizes and sampled counts checked; file fingerprint unchanged during plotting. Full-file hash inherited from prior audit, not recomputed.',
            'frozen_embedding_inputs':{n:sha(BASE/n) for n in ['analysis_report.json','cluster_models.npz','gene_axis.csv','per_cell_full_rank_clusters.csv','pca_variance.csv']},
            'outputs':{n:sha(OUT/n) for n in outputs},'runtime_seconds':time.perf_counter()-started,
            'limits':'Observed cells are non-targeting controls, not real treated truth. Prediction sampling covers all targets but is not a full-cell census. A two-dimensional projection cannot establish perturbation accuracy; pooled target heterogeneity affects prediction dispersion.'}
    (OUT/'analysis_report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    doc=f'''# Observed controls and V16 predicted cells

Plots use all **55,200 observed non-targeting controls** and **18,000 V16 predicted cells**: a seeded simple random sample of 20/400 cells from every one of the 900 target/context groups. Every target is represented in A, B and C. Source counts stay unchanged.

![Cells colored by context](analysis/v16_cell_projection/cells_observed_vs_v16.png)

The top row shows PC1/PC2; the bottom shows PC3/PC4. Columns show observed controls, predicted cells, and the overlay, with identical axes across each row. Blue=A, orange=B, green=C; circles=controls, crosses=predictions. The frozen axes explain {variance[0]:.2%}, {variance[1]:.2%}, {variance[2]:.2%}, and {variance[3]:.2%} of the observed-control rank variance.

## Method and checks

Reuse the original full-matrix rank PCA fitted only on controls. Each cell is ranked across the same 18,533 raw-count genes with average ties and zeros retained, then centered and unit-normalized. Predictions are transformed with the original control mean/components; there is no refit, pseudocount, gene imputation, or separate alignment. This retains our earlier rank-correlation strategy. All 900 groups, gene axes, unique prediction IDs, sampled nonnegative integer counts, and coverage are checked. Eight cached observed coordinates per context are independently replayed against the original counts; predicted projection is checked for batching invariance.

The standalone interactive HTML supports all 300 targets, source/context legend toggles, hover, zoom, and switching PC1/PC2 to PC3/PC4. Observed circles remain unperturbed controls even when a target is selected. There are no measured A/B/C treated responses here. Separation/overlap is descriptive and is not an accuracy score. The prediction pool includes different interventions, so its spread is not directly equivalent to within-treatment biological variance. `projection_summary.csv` reports descriptive eight-PC spread and count summaries, not full-gene variance.

## Reproduce

Run from the root in Windows Command Prompt after restoring the original controls, V16 prediction H5AD and the existing full-rank clustering cache:

```cmd
.\\.competition_env\\Scripts\\python.exe plot_v16_cell_projection.py
```

Outputs: `analysis/v16_cell_projection/cells_observed_vs_v16.png`, `.pdf`, `.html`, `cell_coordinates.csv`, `projection_summary.csv`, and `analysis_report.json`. Large regenerable HTML/coordinate files remain local. Inputs and outputs are recorded with checksums; the prediction file's prior full-audit hash is inherited and explicitly not rehashed by this plotting script. No model or prediction is changed.
'''
    (ROOT/'V16_CELL_PROJECTION.md').write_text(doc,encoding='utf-8',newline='\n')
    print(pd.DataFrame(summary).to_string(index=False),flush=True)
    print(f'Plots complete in {time.perf_counter()-started:.1f}s: {OUT}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int,default=20261012)
    main(seed=p.parse_args().seed)
