"""Audit every supplied ZIP value; cluster 300-target controls and rank markers."""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import csv
import io
import zipfile
import hashlib
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata, mannwhitneyu
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cluster_full_rank_matrix import rank_unit_rows, cluster
from analyze_target_rank_correlations import ROOT, digest, table

OUT=ROOT/'analysis/target_archive_audit'
ZIP=ROOT/'vcc_2026_300_target_matrix.zip'
CONTEXTS=['A','B','C']
PANELS={
 'T_cell':['CD3D','CD3E','TRAC','CD2','LCK','CD4','CD8A','CD8B','IL7R','CCR7','TCF7'],
 'other_immune':['PTPRC','LYZ','LST1','MS4A1','NKG7'],
 'epithelial':['EPCAM','CDH1','KRT8','KRT18','KRT19','KRT7','KRT14','CLDN4'],
 'mesenchymal_associated':['VIM','COL1A1','COL1A2','DCN','LUM','FN1','ACTA2','TAGLN','S100A4','PDGFRB'],
 'proliferation':['MKI67','TOP2A','PCNA','MCM2','MCM5','CCNB1','CDC20','UBE2C','BIRC5'],
 'chromatin_and_requested':['EZH2','SMARCA5','HDAC1','CHD3','S100A11','FXYD5','E2F1']}


def read_archive(path, targets):
    with zipfile.ZipFile(path) as archive:
        files=[x for x in archive.infolist() if not x.is_dir()]
        if len(files)!=1 or not files[0].filename.endswith('.csv'):
            raise ValueError('Expected one unambiguous CSV member; no archive extraction is used')
        member=files[0]
        if archive.testzip() is not None:
            raise ValueError('ZIP CRC failure')
        with archive.open(member) as handle:
            header=next(csv.reader(io.TextIOWrapper(handle,encoding='utf-8-sig')))
        if len(header)!=301 or header[0]!='context' or len(set(header))!=len(header):
            raise ValueError('Expected context plus 300 unique named gene columns')
        missing=sorted(set(targets)-set(header[1:]));unexpected=sorted(set(header[1:])-set(targets))
        if missing or unexpected:
            raise ValueError(f'Incorrect gene list: missing={missing}, unexpected={unexpected}')
        sha=hashlib.sha256()
        with archive.open(member) as handle:
            for block in iter(lambda:handle.read(8*1024**2),b''):sha.update(block)
        with archive.open(member) as handle:frame=pd.read_csv(handle)
        if len(frame)!=55200 or frame.context.value_counts().to_dict()!={'A':18400,'B':18400,'C':18400}:
            raise ValueError('Unexpected row count or context labels')
        raw=frame[list(targets)].to_numpy(dtype=np.float64)
        if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any():
            raise ValueError('Expected finite nonnegative integer-valued raw counts')
        return frame, raw, {'file':Path(path).name,'sha256':digest(path),'member':member.filename,
                'member_sha256':sha.hexdigest(),'member_bytes':member.file_size,'crc32':member.CRC,
                'zip_crc_valid':True,'rows':len(frame),'gene_columns':len(targets),
                'missing_targets':missing,'unexpected_genes':unexpected,'duplicate_columns':False,
                'gene_order_matches_target_list':header[1:]==list(targets),
                'finite_nonnegative_integer_values':True,'cell_ids_in_archive':False}


def marker_statistics(counts, libraries, groups, genes):
    """Descriptive gene markers using full-library CPM, without any added count."""
    counts=np.asarray(counts,dtype=np.float64);libraries=np.asarray(libraries,dtype=float)
    groups=np.asarray(groups)
    if not np.isfinite(counts).all() or np.any(counts<0) or np.any(libraries<=0) or counts.shape[0]!=len(groups):
        raise ValueError('Invalid marker inputs')
    cpm=counts/libraries[:,None]*1e6
    names=sorted(set(groups));means={g:cpm[groups==g].mean(axis=0) for g in names}
    detection={g:(counts[groups==g]>0).mean(axis=0) for g in names}
    aucs={}
    for i,left in enumerate(names):
        for right in names[i+1:]:
            x=cpm[groups==left];y=cpm[groups==right]
            auc=np.asarray(mannwhitneyu(x,y,axis=0,method='asymptotic').statistic)/(len(x)*len(y))
            aucs[left,right]=auc;aucs[right,left]=1-auc
    ranks=rankdata(cpm,axis=0,method='average')
    rows=[]
    for group in names:
        mask=groups==group;n=int(mask.sum());rest=int((~mask).sum())
        other=[g for g in names if g!=group]
        denominator=np.max(np.array([means[g] for g in other]),axis=0)
        fold=np.divide(means[group],denominator,out=np.full(len(genes),np.nan),where=denominator>0)
        fold[(denominator==0)&(means[group]>0)]=np.inf
        minimum=np.min(np.array([aucs[group,g] for g in other]),axis=0)
        max_detection=np.max(np.array([detection[g] for g in other]),axis=0)
        rest_auc=(ranks[mask].sum(axis=0)-n*(n+1)/2)/(n*rest)
        frame=pd.DataFrame({'group':group,'gene_name':genes,'cells':n,
                'mean_raw_count':counts[mask].mean(axis=0),'median_raw_count':np.median(counts[mask],axis=0),
                'detection_fraction':detection[group],'mean_cpm':means[group],
                'median_cpm':np.median(cpm[mask],axis=0),'max_other_mean_cpm':denominator,
                'fold_over_highest_other_group':fold,'log2_fold_no_pseudocount':np.log2(fold,where=fold>0,out=np.full(len(genes),np.nan)),
                'auc_vs_rest':rest_auc,'minimum_pairwise_auc':minimum,'max_other_detection_fraction':max_detection,
                'panel_enriched':(minimum>=.8)&(fold>=2)&(detection[group]>=.25),
                'detection_restricted':(detection[group]>=.5)&(max_detection<=.05)})
        for name in other:frame['auc_vs_'+name]=aucs[group,name]
        frame=frame.sort_values(['minimum_pairwise_auc','fold_over_highest_other_group','gene_name'],ascending=[False,False,True])
        frame['marker_rank']=np.arange(1,len(frame)+1)
        rows.append(frame)
    return pd.concat(rows,ignore_index=True)


def raw_summary(raw,libraries,groups,genes):
    rows=[]
    for group in CONTEXTS:
        mask=np.asarray(groups)==group;x=raw[mask];cpm=x/np.asarray(libraries)[mask,None]*1e6
        rows.append(pd.DataFrame({'context':group,'gene_name':genes,'mean_raw_count':x.mean(axis=0),
            'median_raw_count':np.median(x,axis=0),'detection_fraction':(x>0).mean(axis=0),
            'mean_cpm':cpm.mean(axis=0),'median_cpm':np.median(cpm,axis=0)}))
    return pd.concat(rows,ignore_index=True)


def make_plot(cells,composition,markers,depth_markers,broad,explained):
    fig,axes=plt.subplots(2,3,figsize=(17,10),layout='constrained')
    shuffled=cells.sample(frac=1,random_state=20261004);colors=['#4477aa','#ee7733','#228833']
    for context,color in zip(CONTEXTS,colors):
        sub=shuffled[shuffled.context==context]
        axes[0,0].scatter(sub.PC1,sub.PC2,s=1,alpha=.2,c=color,rasterized=True,label=context)
        sub=shuffled[shuffled.assigned_context==context]
        axes[0,1].scatter(sub.PC1,sub.PC2,s=1,alpha=.2,c=color,rasterized=True,label=context)
    axes[0,0].legend(markerscale=7,title='Original context');axes[0,1].legend(markerscale=7,title='Blind cluster name')
    axes[0,0].set_title('All 55,200 cells: supplied context');axes[0,1].set_title('Clustering uses only the 300 targets')
    image=axes[0,2].scatter(shuffled.PC1,shuffled.PC2,c=np.log10(shuffled.full_library_counts),s=1,alpha=.2,cmap='viridis',rasterized=True)
    fig.colorbar(image,ax=axes[0,2],label='log10(full-library raw counts)');axes[0,2].set_title('Count-depth variation')
    for ax in axes[0]:ax.set_xlabel(f'PC1 ({explained[0]:.1%})');ax.set_ylabel(f'PC2 ({explained[1]:.1%})')
    matrix=composition.pivot(index='context',columns='assigned_context',values='cells').reindex(index=CONTEXTS,columns=CONTEXTS)
    axes[1,0].imshow(matrix,cmap='Blues',aspect='auto')
    for i in range(3):
        for j in range(3):axes[1,0].text(j,i,f'{matrix.iloc[i,j]:,}',ha='center',va='center',color='white' if matrix.iloc[i,j]>10000 else 'black')
    axes[1,0].set_xticks(range(3),CONTEXTS);axes[1,0].set_yticks(range(3),CONTEXTS)
    axes[1,0].set_xlabel('Assigned cluster name');axes[1,0].set_ylabel('Supplied context');axes[1,0].set_title('Complete cell counts; labels attached after fitting')
    selected=markers[markers.panel_enriched].groupby('group',sort=True).head(5).gene_name.drop_duplicates().tolist()
    matrix=markers.pivot(index='gene_name',columns='group',values='mean_cpm').loc[selected,CONTEXTS]
    scaled=(matrix.to_numpy()-matrix.to_numpy().mean(axis=1,keepdims=True))/np.maximum(matrix.to_numpy().std(axis=1,keepdims=True),1e-12)
    axes[1,1].imshow(scaled,cmap='RdBu_r',vmin=-1.5,vmax=1.5,aspect='auto')
    axes[1,1].set_xticks(range(3),CONTEXTS);axes[1,1].set_yticks(range(len(selected)),selected,fontsize=8)
    axes[1,1].set_title('Top panel markers: gene-wise standardized mean CPM')
    requested=['VIM','S100A11','FXYD5','EZH2','SMARCA5']
    for i,gene in enumerate(requested):
        sub=depth_markers[depth_markers.gene_name==gene]
        for context,color in zip(CONTEXTS,colors):
            part=sub[sub.context==context]
            axes[1,2].plot(part.depth_quintile,part.mean_cpm/1000,color=color,linestyle=['-','--',':','-.','-'][i],alpha=.7,label=f'{gene} {context}')
    axes[1,2].set_xlabel('Within-context full-depth quintile');axes[1,2].set_ylabel('Mean CPM / 1,000')
    axes[1,2].set_title('Requested markers across depth strata');axes[1,2].set_xticks(range(1,6))
    axes[1,2].legend(fontsize=6,ncol=3)
    fig.suptitle('Verified 300-target matrix: raw-rank clustering and descriptive gene markers\nZeros retained; no pseudocounts; marker CPM uses the full 18,533-gene library',fontsize=14)
    fig.savefig(OUT/'target300_clusters_markers.png',dpi=180);fig.savefig(OUT/'target300_clusters_markers.pdf');plt.close(fig)
    selected=['CD3D','CD3E','TRAC','VIM','KRT8','KRT18','EPCAM','CDH1','S100A11','FXYD5','EZH2','SMARCA5','MKI67','TOP2A','CCNB1','UBE2C']
    fig,axes=plt.subplots(1,2,figsize=(13,7),layout='constrained')
    values=broad.pivot(index='gene_name',columns='context',values='mean_cpm').loc[selected,CONTEXTS].to_numpy()
    scaled=np.divide(values,values.max(axis=1,keepdims=True),out=np.zeros_like(values),where=values.max(axis=1,keepdims=True)>0)
    axes[0].imshow(scaled,cmap='viridis',vmin=0,vmax=1,aspect='auto');axes[0].set_title('Expression relative to each gene maximum\nFull-matrix lineage/proliferation check')
    detection=broad.pivot(index='gene_name',columns='context',values='detection_fraction').loc[selected,CONTEXTS].to_numpy()
    axes[1].imshow(detection,cmap='Blues',vmin=0,vmax=1,aspect='auto');axes[1].set_title('Fraction of cells with measured count > 0')
    for ax in axes:ax.set_xticks(range(3),CONTEXTS);ax.set_yticks(range(len(selected)),selected)
    for i in range(len(selected)):
        for j in range(3):
            axes[0].text(j,i,f'{values[i,j]:.1f}',ha='center',va='center',fontsize=8,color='black' if scaled[i,j]>.6 else 'white')
            axes[1].text(j,i,f'{detection[i,j]:.1%}',ha='center',va='center',fontsize=8,color='white' if detection[i,j]>.5 else 'black')
    fig.suptitle('Biological interpretation uses the original full matrix\nThese additional markers do not enter the 300-gene clustering; left labels are mean CPM',fontsize=12)
    fig.savefig(OUT/'broader_lineage_markers.png',dpi=180);fig.savefig(OUT/'broader_lineage_markers.pdf');plt.close(fig)


def run(seed=20261004,threads=12):
    OUT.mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    frame,raw,audit=read_archive(ZIP,targets)
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy();position={g:i for i,g in enumerate(genes)}
    requested=list(dict.fromkeys(g for names in PANELS.values() for g in names))
    broader_genes=np.array([g for g in requested if g in position])
    audit['broader_markers_unmeasured']=[g for g in requested if g not in position]
    cells=pd.DataFrame({'archive_row_index':np.arange(len(frame)),'context':frame.context.to_numpy()})
    cells['cell_id']='';cells['ntc_id']='';cells['full_library_counts']=0
    broad_raw=np.empty((len(raw),len(broader_genes)))
    old=json.loads((ROOT/'analysis/full_rank_clustering/analysis_report.json').read_text())
    input_records=[];total_mismatches=0
    for context in CONTEXTS:
        path=ROOT/f'context_{context}.h5ad';sha=digest(path)
        if sha!=next(x['sha256'] for x in old['inputs'] if x['file']==path.name):raise ValueError('Original input changed')
        obj=ad.read_h5ad(path,backed='r');indices=np.flatnonzero(frame.context.to_numpy()==context)
        mismatches=0;max_difference=0
        try:
            if obj.shape!=(18400,18533) or not np.array_equal(obj.var_names,genes):raise ValueError('Original axes differ')
            cells.loc[indices,'cell_id']=obj.obs_names.to_numpy();cells.loc[indices,'ntc_id']=obj.obs.ntc_id.astype(str).to_numpy()
            for start in range(0,len(obj),512):
                end=min(start+512,len(obj));original=sparse.csr_matrix(obj.X[start:end],dtype=np.float64)
                expected=original[:,[position[g] for g in targets]].toarray()
                actual=raw[indices[start:end]]
                mismatches+=int(np.sum(actual!=expected));max_difference=max(max_difference,float(np.max(np.abs(actual-expected))))
                broad_raw[indices[start:end]]=original[:,[position[g] for g in broader_genes]].toarray()
                cells.loc[indices[start:end],'full_library_counts']=np.asarray(original.sum(axis=1)).ravel().astype(np.int64)
            input_records.append({'file':path.name,'sha256':sha,'context':context,'checked_gene_values':len(obj)*len(targets),
                                  'mismatched_values':mismatches,'max_absolute_difference':max_difference})
            total_mismatches+=mismatches
            print(f'Checked {context}: {len(obj)*len(targets):,} values; mismatches={mismatches}',flush=True)
        finally:obj.file.close()
    audit.update({'original_inputs':input_records,'values_checked':int(raw.size),'mismatched_values':total_mismatches,
                  'status':'passed' if total_mismatches==0 else 'failed',
                  'row_identity_policy':'CSV contains no cell IDs; exact source row order verified independently within each context; original IDs attached from that ordering'})
    (OUT/'archive_validation.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    if total_mismatches:raise ValueError('ZIP values differ from original counts; inspect archive_validation.json')
    libraries=cells.full_library_counts.to_numpy()
    z=rank_unit_rows(raw);rng=np.random.default_rng(seed);sample=np.sort(rng.choice(len(z),1500,replace=False))
    metrics=[];models={};histories={}
    with threadpool_limits(limits=threads):
        distance=np.clip(1-z[sample]@z[sample].T,0,2);np.fill_diagonal(distance,0)
        for k,repeat in [(2,0),(3,0),(4,0),(6,0),(3,1)]:
            labels,centers,history,done=cluster(z,k,seed+repeat)
            if not done:raise ValueError('Target-only spherical clustering failed to converge')
            key=f'k{k}_repeat{repeat}';models[key]=centers;histories[key]=history;cells['cluster_'+key]=labels
            metrics.append({'k':k,'repeat':repeat,'seed':seed+repeat,'converged':done,'iterations':len(history),
                   'context_ARI':float(adjusted_rand_score(cells.context,labels)),
                   'context_NMI':float(normalized_mutual_info_score(cells.context,labels)),
                   'sample_silhouette_spearman':float(silhouette_score(distance,labels[sample],metric='precomputed')),
                   'minimum_cluster_cells':int(np.bincount(labels,minlength=k).min())})
            print('Clustered',key,'ARI',metrics[-1]['context_ARI'],flush=True)
        labels=cells.cluster_k3_repeat0.to_numpy();cross=pd.crosstab(cells.context,labels).reindex(index=CONTEXTS,columns=range(3),fill_value=0)
        rows,cols=linear_sum_assignment(-cross.to_numpy());order=cols[np.argsort(rows)]
        assignment={int(k):CONTEXTS[i] for i,k in enumerate(order)}
        cells['assigned_context']=[assignment[int(k)] for k in labels]
        centers=models['k3_repeat0'][order];similarity=z@centers.T;sorted_similarity=np.sort(similarity,axis=1)
        cells['best_centroid_similarity']=sorted_similarity[:,-1];cells['centroid_margin']=sorted_similarity[:,-1]-sorted_similarity[:,-2]
        cells['matches_context']=cells.context==cells.assigned_context
        cells['nonzero_target_genes']=(raw>0).sum(axis=1);cells['target_count_sum']=raw.sum(axis=1).astype(np.int64)
        pca=PCA(n_components=8,svd_solver='randomized',iterated_power=3,random_state=seed)
        pca.fit(z);pcs=pca.transform(z)
        for i in range(8):cells[f'PC{i+1}']=pcs[:,i]
        models['pca_mean']=pca.mean_;models['pca_components']=pca.components_;models['context_order_centroids']=centers
        print('Ranking descriptive depth-adjusted gene markers',flush=True)
        cluster_markers=marker_statistics(raw,libraries,cells.assigned_context.to_numpy(),targets)
        context_markers=marker_statistics(raw,libraries,cells.context.to_numpy(),targets)
    composition=pd.crosstab(cells.context,cells.assigned_context).reindex(index=CONTEXTS,columns=CONTEXTS,fill_value=0).rename_axis(index='context',columns='assigned_context').stack().reset_index(name='cells')
    broader=raw_summary(broad_raw,libraries,cells.context,broader_genes)
    group_of={g:group for group,names in PANELS.items() for g in names}
    broader['marker_panel']=[group_of[g] for g in broader.gene_name];broader['in_300_targets']=broader.gene_name.isin(targets)
    requested_targets=['VIM','S100A11','FXYD5','EZH2','SMARCA5','E2F1']
    top=cluster_markers[cluster_markers.panel_enriched].groupby('group',sort=True).head(10)
    selected=list(dict.fromkeys(top.gene_name.tolist()+requested_targets));depth_rows=[]
    for context in CONTEXTS:
        mask=cells.context.to_numpy()==context;idx=np.flatnonzero(mask)
        quintiles=pd.qcut(pd.Series(libraries[idx]).rank(method='first'),5,labels=False).to_numpy()+1
        for quintile in range(1,6):
            rows=idx[quintiles==quintile];sub=raw[rows][:,[list(targets).index(g) for g in selected]];cpm=sub/libraries[rows,None]*1e6
            depth_rows.append(pd.DataFrame({'context':context,'depth_quintile':quintile,'gene_name':selected,
                'cells':len(rows),'mean_cpm':cpm.mean(axis=0),'mean_raw_count':sub.mean(axis=0),'detection_fraction':(sub>0).mean(axis=0)}))
    depth=pd.concat(depth_rows,ignore_index=True)
    gene_summary=raw_summary(raw,libraries,cells.context,targets)
    mismatched=cells[~cells.matches_context].copy()
    original_cells=pd.read_csv(ROOT/'analysis/full_rank_clustering/per_cell_full_rank_clusters.csv',usecols=['cell_id','context','cluster_k3'])
    mismatched=mismatched.merge(original_cells.rename(columns={'context':'original_full_gene_context','cluster_k3':'original_full_gene_cluster'}),on='cell_id',how='left',validate='one_to_one')
    outputs={'per_cell_target300_clusters.csv':cells,'cluster_composition.csv':composition,
            'cluster_metrics.csv':pd.DataFrame(metrics),'target300_cluster_markers.csv':cluster_markers,
            'target300_context_markers.csv':context_markers,'top10_enriched_cluster_markers.csv':top,
            'target300_expression_summary.csv':gene_summary,'broader_lineage_expression_summary.csv':broader,
            'marker_depth_quintiles.csv':depth,'cells_disagreeing_with_context.csv':mismatched,
            'gene_axis.csv':pd.DataFrame({'gene_name':targets}),
            'pca_variance.csv':pd.DataFrame({'PC':np.arange(1,9),'explained_variance_ratio':pca.explained_variance_ratio_})}
    for name,data in outputs.items():data.to_csv(OUT/name,index=False)
    np.savez_compressed(OUT/'target300_rank_models.npz',**models)
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'archive_sha256':audit['sha256'],
       'cells':len(cells),'targets':len(targets),'zip_value_mismatches':total_mismatches,'seed':seed,'threads':threads,
       'context_ARI':float(adjusted_rand_score(cells.context,labels)),'repeat_ARI':float(adjusted_rand_score(cells.cluster_k3_repeat0,cells.cluster_k3_repeat1)),
       'correctly_grouped_cells':int(cells.matches_context.sum()),'disagreeing_cells':len(mismatched),
       'cluster_to_context':assignment,'context_inputs':input_records,'refinement_histories':histories,
       'rank_method':'All 300 genes, raw average tied ranks, centered unit vectors; no pseudocount/log/CPM in clustering',
       'marker_method':'Full-library CPM for descriptive AUC and mean ratios, raw-count detection and means; no pseudocounts; no population p-values',
       'panel_enriched_rule':'Minimum AUC against either other cluster >=0.8, CPM fold >=2 over largest other mean, detection >=0.25',
       'detection_restricted_rule':'Detection >=0.5 in group and <=0.05 in either other group; not tissue specificity',
       'panel_enriched_counts':cluster_markers.groupby('group').panel_enriched.sum().astype(int).to_dict(),
       'detection_restricted_counts':cluster_markers.groupby('group').detection_restricted.sum().astype(int).to_dict(),
       'targets_with_no_enriched_group':int(len(set(targets)-set(cluster_markers.loc[cluster_markers.panel_enriched,'gene_name']))),
       'versions':{name:version(name) for name in ['numpy','scipy','pandas','anndata','scikit-learn','matplotlib']},
       'scripts':{name:digest(ROOT/name) for name in ['analyze_target_archive.py','test_target_archive.py','validate_target_archive.py']}}
    make_plot(cells,composition,cluster_markers,depth,broader,pca.explained_variance_ratio_)
    report['outputs']={p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in OUT.iterdir() if p.is_file() and p.name not in ['analysis_report.json','validation.json']}
    (OUT/'analysis_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    text='# Verified 300-target archive: clustering and markers\n\n'
    text+=f"The ZIP contains every one of the 300 requested target genes, all 55,200 controls, and **{audit['values_checked']:,} original raw count values with zero mismatches**. Float CSV notation such as `2.0` preserves the integer count 2. There are no cell IDs in the archive; original IDs were attached using the verified row order within each context.\n\n"
    text+='## Target-only rank clustering\n\n'+table(composition)+'\n\n'
    text+=f"Correctly grouped cells: **{report['correctly_grouped_cells']:,}/{len(cells):,} ({report['correctly_grouped_cells']/len(cells):.4%})**. Context ARI: {report['context_ARI']:.6f}; independent-seed repeat ARI: {report['repeat_ARI']:.6f}. Cluster names are attached after fitting; neither context labels nor other genes enter the clustering. All assignments and disagreements are saved.\n\n"
    text+='## Strong panel-enriched markers\n\n'+table(top[['group','gene_name','mean_raw_count','detection_fraction','mean_cpm','fold_over_highest_other_group','minimum_pairwise_auc']])+'\n\n'
    text+='Markers are enriched relative to the other supplied clusters under explicit descriptive thresholds, not proven tissue-specific genes. Full 900-row tables include every gene/group, pairwise AUCs, CPM fold ratios without pseudocounts and detection-restriction flags. The 300 CRISPRi targets form a selected panel rather than an unbiased collection of lineage markers.\n\n'
    text+='## The proposed diagnostic genes\n\n'+table(gene_summary[gene_summary.gene_name.isin(requested_targets)])+'\n\n'
    text+='## Broader identity and proliferation check\n\n'+table(broader[broader.gene_name.isin(['CD3D','CD3E','TRAC','KRT8','KRT18','EPCAM','CDH1','MKI67','TOP2A','CCNB1','UBE2C'])][['context','gene_name','mean_raw_count','detection_fraction','mean_cpm','in_300_targets']])+'\n\n'
    text+='The full-matrix CD3D/CD3E/TRAC program supports a T-cell-like interpretation for A. B expresses both VIM and epithelial keratins, supporting a mixed expression program rather than proof of mesenchymal origin or an observed EMT transition. C expresses EPCAM/CDH1 and keratins, providing stronger epithelial evidence than S100A11/FXYD5 alone. Chromatin/proliferation expression is a cell-state signal and does not identify tissue of origin. These are marker-based inferences, not established anonymous identities. Literature and limitations are documented in [TARGET_LINEAGE_INTERPRETATION.md](TARGET_LINEAGE_INTERPRETATION.md).\n\n'
    text+='![Target clustering](analysis/target_archive_audit/target300_clusters_markers.png)\n\n![Broader markers](analysis/target_archive_audit/broader_lineage_markers.png)\n\nReproduction and audit details: [TARGET_ARCHIVE_AUDIT.md](TARGET_ARCHIVE_AUDIT.md).\n'
    (ROOT/'TARGET_ARCHIVE_RESULTS.md').write_text(text,encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['context_ARI','repeat_ARI','correctly_grouped_cells','disagreeing_cells','panel_enriched_counts','detection_restricted_counts','targets_with_no_enriched_group']},indent=2),flush=True)
    print(top[['group','gene_name','minimum_pairwise_auc','fold_over_highest_other_group']].to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--seed',type=int,default=20261004);parser.add_argument('--threads',type=int,default=12)
    args=parser.parse_args();run(args.seed,args.threads)
