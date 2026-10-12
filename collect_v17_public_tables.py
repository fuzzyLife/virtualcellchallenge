"""Collect primary iPSC and KOLF manifests and exact target-coverage audits."""
from pathlib import Path
import json,hashlib,requests
import numpy as np,pandas as pd
import build_submission as b
A=b.ROOT/'experiments/v17_learned_heterogeneity/analysis'
def fetch(url,file):
 if not file.exists():
  r=requests.get(url,timeout=180);r.raise_for_status();file.write_bytes(r.content)
 return file

def main():
 A.mkdir(parents=True,exist_ok=True);_,genes,targets=b.read_spec(b.ROOT);target_set=set(targets);records=[];sources=[]
 for article,name in [(27989294,'ipsc'),(27261219,'kolf')]:
  url=f'https://api.figshare.com/v2/articles/{article}';file=fetch(url,A/f'{name}_figshare_metadata.json');d=json.loads(file.read_text());sources.append({'name':name,'manifest_url':url,'manifest_sha256':b.sha256(file),'doi':d.get('doi'),'version':d.get('version'),'files':d['files']})
  if name=='ipsc':
   for item in d['files']:
    if not item['name'].endswith('Cell-Metadata.tsv.gz'):continue
    p=fetch(item['download_url'],A/item['name']);assert p.stat().st_size==item['size'];assert hashlib.md5(p.read_bytes()).hexdigest()==item['computed_md5']
    m=pd.read_csv(p,sep='\t');gene=m.Guide_Call.str.extract(r'^([A-Za-z0-9.-]+)_[ACGT]{19,23}$',expand=False);counts=gene.value_counts();over=sorted(target_set&set(counts.index));group=m.assign(target=gene);line_power=group[group.target.isin(targets)].groupby(['target','Cell_Line']).size()
    for target in targets:records.append({'dataset':item['name'],'target_gene':target,'present':target in over,'pooled_assigned_cells':int(counts.get(target,0)),'maximum_cells_in_one_donor':int(line_power.loc[target].max()) if target in line_power.index.get_level_values(0) else 0})
    sources.append({'name':item['name'],'url':item['download_url'],'sha256':b.sha256(p),'cells':len(m),'unambiguous_single_guide_cells':int(gene.notna().sum()),'target_overlap':len(over),'competition_target_cells':int(counts.reindex(over).sum()),'targets_with_400_pooled_cells':int((counts.reindex(over)>=400).sum()),'target_cells_median':float(counts.reindex(over).median()),'note':'Unassigned cells are NOT non-targeting controls. Count matrix and matched donor/batch controls needed before training.'})
 treefile=fetch('https://api.github.com/repos/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas/git/trees/main?recursive=1',A/'kolf_repo_tree.json');tree=json.loads(treefile.read_text());sha=tree['sha'];files=['docs/data/perts/index.json','docs/data/panel/all/meta.json','docs/data/clusters/summary.json','psp/notebooks/mde_webapp/webapp_build/build_panel_all_perts.py','psp/notebooks/mde_webapp/webapp_build/build_panel_pert_means.py']
 for name in files:
  path=A/('kolf_'+Path(name).parent.name+'_'+Path(name).name);url=f'https://raw.githubusercontent.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas/{sha}/{name}';fetch(url,path);sources.append({'name':name,'revision':sha,'url':url,'sha256':b.sha256(path)})
 license_url=f'https://raw.githubusercontent.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas/{sha}/LICENSE';license_file=fetch(license_url,A/'kolf_LICENSE.txt');sources.append({'name':'KOLF author MIT license','revision':sha,'url':license_url,'sha256':b.sha256(license_file)})
 meta=json.loads((A/'kolf_all_meta.json').read_text());over=sorted(target_set&set(meta['perts']))
 for t in targets:records.append({'dataset':'KOLF_browser_NTC_z_score','target_gene':t,'present':t in over,'pooled_assigned_cells':None,'maximum_cells_in_one_donor':None})
 sources.append({'name':'KOLF browser','target_overlap':len(over),'output_gene_overlap':len(set(genes)&set(meta['genes'])),'matrix_encoding':meta,'note':'Quantized mean NTC z-scores, not raw counts or log fold changes; cannot estimate cell-level variance from these means. Full unselected AnnData counts 189393177972 bytes; not downloaded.'})
 # Keep compact coverage and provenance; the large cell tables are ignored.
 pd.DataFrame(records).to_csv(A/'new_public_target_coverage.csv',index=False);b.write_json(A/'new_public_data_inventory.json',{'date':'2026-10-12','collector_sha256':b.sha256(Path(__file__)),'coverage_sha256':b.sha256(A/'new_public_target_coverage.csv'),'sources':sources,'used_for_V17_fitting':False,'reason':'Additional-source discovery separate from frozen V17 fitting; no substitution of z-scores for counts, no pooling donors as identical biological replicates.'});print(pd.DataFrame(records).groupby('dataset').present.sum().to_string(),flush=True)
if __name__=='__main__':main()
