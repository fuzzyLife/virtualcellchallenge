"""Collect bounded public research inputs; pin bytes and never modify V11."""
import argparse,hashlib,io,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,requests
import build_submission as build
ROOT=build.ROOT; OUT=ROOT/'analysis/post_v11_research'; DATA=OUT/'data'
URLS={
 'reactome_gmt.zip':'https://reactome.org/download/current/ReactomePathways.gmt.zip',
 'reactome_version.txt':'https://reactome.org/ContentService/data/database/version',
 'viperturb_record.json':'https://zenodo.org/api/records/18460279',
 'replogle_record.json':'https://api.figshare.com/v2/articles/20029387'}

def fetch(name,url,records,limit=8*1024**2,expected_md5=None):
 path=DATA/name
 if path.exists():
  assert name in records and build.sha256(path)==records[name]['sha256'],name
  return path
 with requests.get(url,stream=True,timeout=(15,45)) as response:
  response.raise_for_status();parts=[];size=0
  for chunk in response.iter_content(65536):
   size+=len(chunk)
   if size>limit:raise ValueError('Bounded download exceeded: '+name)
   parts.append(chunk)
  content=b''.join(parts)
  if expected_md5:assert hashlib.md5(content).hexdigest()==expected_md5,name
  path.write_bytes(content)
  records[name]={'url':url,'resolved_url':response.url.split('?')[0],'bytes':size,'sha256':build.sha256(path),
   'retrieved_utc':datetime.now(timezone.utc).isoformat(),'etag':response.headers.get('ETag'),'publisher_md5':expected_md5}
  build.write_json(OUT/'data_manifest.json',records)
 print('Collected',name,size,'bytes',flush=True);return path

def main():
 DATA.mkdir(parents=True,exist_ok=True); mp=OUT/'data_manifest.json'; records=json.loads(mp.read_text()) if mp.exists() else {}
 for name,url in URLS.items():
  try:fetch(name,url,records)
  except requests.RequestException as error:print('Unavailable',name,str(error),flush=True)
 genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy(str)
 if (DATA/'reactome_gmt.zip').exists():
  with zipfile.ZipFile(DATA/'reactome_gmt.zip') as archive:
   members=[n for n in archive.namelist() if n.endswith('.gmt')];assert len(members)==1
   info=archive.getinfo(members[0]);assert info.file_size<32*1024**2
   text=archive.read(members[0]).decode('utf-8-sig')
  memberships=[];pathways=[]
  for line in text.splitlines():
   fields=line.split('\t');name,pid=fields[:2];symbols=set(fields[2:])
   if not pid.startswith('R-HSA-'):continue
   measured=sorted(symbols&set(genes));panel=sorted(symbols&set(targets))
   pathways.append({'pathway_id':pid,'pathway_name':name,'all_genes':len(symbols),'challenge_genes':len(measured),'panel_targets':len(panel)})
   memberships.extend({'pathway_id':pid,'pathway_name':name,'gene_name':g,'is_panel_target':g in set(targets)} for g in measured)
  pd.DataFrame(pathways).to_csv(OUT/'reactome_pathways.csv',index=False)
  frame=pd.DataFrame(memberships);frame.to_csv(OUT/'reactome_membership.csv',index=False)
  count=frame.groupby('gene_name').size();pd.DataFrame({'target_gene':targets,'reactome_pathways':[int(count.get(t,0)) for t in targets]}).to_csv(OUT/'target_annotation_coverage.csv',index=False)
  print('Reactome pathways',len(pathways),'annotated targets',sum(count.get(t,0)>0 for t in targets),flush=True)
 if (DATA/'viperturb_record.json').exists():
  record=json.loads((DATA/'viperturb_record.json').read_text());files=record['files']
  inventory=[{'filename':f['key'],'bytes':f['size'],'checksum':f['checksum'],'url':f['links']['self']} for f in files]
  pd.DataFrame(inventory).to_csv(OUT/'viperturb_inventory.csv',index=False)
  file=next(f for f in files if f['key']=='genome_wide_manifest.txt')
  path=fetch('viperturb_target_manifest.txt',file['links']['self'],records,limit=1024**2,expected_md5=file['checksum'].removeprefix('md5:'))
  # Store the exact original manifest. Exact symbol tokens only, never substring matches.
  manifest=pd.read_csv(path,sep=r'\s+');manifest=manifest[manifest.gene.isin(targets)].copy();assert manifest.gene.is_unique
  coverage=pd.DataFrame({'target_gene':targets}).merge(manifest.rename(columns={'gene':'target_gene'}),on='target_gene',how='left',validate='one_to_one')
  coverage['manifest_exact_token']=coverage.cell_count.notna()
  coverage.to_csv(OUT/'viperturb_target_coverage.csv',index=False)
  build.write_json(OUT/'viperturb_assessment.json',{'doi':record['metadata']['doi'],'license':record['metadata'].get('license'),
   'manifest_exact_target_overlap':int(coverage.manifest_exact_token.sum()),'panel_targets':len(targets),'qc_filtered_panel_targets':int(coverage.is_in_filtered_obj.fillna(False).sum()),'median_available_target_cells':float(coverage.cell_count.median()),'targets_with_400_available_cells':int((coverage.cell_count>=400).sum()),
   'scope':'Manifest exact-token eligibility only. No expression, usable-cell, output-gene or knockdown-efficacy validation. Shared controls across bin objects must be deduplicated.',
   'large_expression_files_downloaded':False})
  print('VIPerturb exact-token overlap',int(coverage.manifest_exact_token.sum()),flush=True)
 fetch('fingerprint_commit.json','https://api.github.com/repos/satijalab/rna-fingerprinting/commits/main',records,limit=1024**2)
 commit=json.loads((DATA/'fingerprint_commit.json').read_text())['sha']
 fetch('fingerprint_tree.json',f'https://api.github.com/repos/satijalab/rna-fingerprinting/git/trees/{commit}?recursive=1',records,limit=1024**2)
 for remote,local in [('README.md','fingerprint_readme.md'),('LICENSE.md','fingerprint_license.md'),('R/fit.R','fingerprint_fit.R'),('R/utilities.R','fingerprint_utilities.R')]:
  fetch(local,f'https://raw.githubusercontent.com/satijalab/rna-fingerprinting/{commit}/{remote}',records,limit=2*1024**2)
 build.write_json(OUT/'collection_run.json',{'script_sha256':build.sha256(Path(__file__)),'data_manifest_sha256':build.sha256(mp),
  'scope':'Small public annotations and source metadata only; no challenge treated labels or model changes.'})

if __name__=='__main__':main()
