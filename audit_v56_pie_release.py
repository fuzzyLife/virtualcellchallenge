"""Reproduce the V56 PIE release/coverage audit; never load weights or execute PIE.

First run --fetch to snapshot bounded public metadata. Subsequent runs verify the
saved hashes and recompute coverage offline. No competition submission is made.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, hashlib, json, subprocess
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/v56_pie_release'
OLD = ROOT / 'analysis/post_v17_research'
CODE = ROOT / 'pie'
WD = '055c7a2cabdaff5659121c64fa2b609c231fc5be'
XD = 'f6da9607226139a1040b6992736118c9b5c7c313'
SPLITS = '396ab9563175ee887750c9eed7ccaea6f5fdbf50'
META_REVS = {
 'replogle_nadig_essential': '20c9faef76fc96fdc809871340bd93499413dfe7',
 'tahoe100m': '597d6ded9f9a148331388484ffb193489121c702',
 'jiang': '99242356dfcf58e049d867703c5743df8d1790bf',
 'arc_vcc_25': 'e7dda6065316959958b87583ddf0a15f6e69c37d',
 'x_atlas_orion': '0950c8aa5ce5b33dc9d6f0c5eade3ae445ca0ffb',
}

def sha(path):
 return hashlib.sha256(path.read_bytes()).hexdigest()

def write(path, obj):
 path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')

def git(*args):
 safe = str(CODE).replace('\\', '/')
 return subprocess.check_output(['git', '-c', 'safe.directory='+safe, '-C', str(CODE), *args], text=True).strip()

def sources():
 u = {
  'github_head.json': 'https://api.github.com/repos/ArcInstitute/pie/commits?per_page=1',
  'wd_current.json': 'https://huggingface.co/api/models/arcinstitute/PIE_replogle_wdataset?blobs=true',
  'xd_current.json': 'https://huggingface.co/api/models/arcinstitute/PIE_replogle_xdataset?blobs=true',
  'wd_train.json': f'https://huggingface.co/datasets/arcinstitute/PIE_splits/resolve/{SPLITS}/replogle_wdataset/unseen_ctx/k562/train.json',
  'xd_train.json': f'https://huggingface.co/datasets/arcinstitute/PIE_splits/resolve/{SPLITS}/replogle_xdataset/train.json',
 }
 for k,v in META_REVS.items():
  u[k+'_meta.json'] = f'https://huggingface.co/datasets/arcinstitute/PIE_{k}/resolve/{v}/preprocessed/meta.json'
 return u

def fetch():
 OUT.mkdir(parents=True, exist_ok=True)
 manifest_path = OUT/'source_manifest.json'
 if manifest_path.exists():
  raise RuntimeError('Snapshot exists; run without --fetch to verify, or use a new audit namespace for updates.')
 entries=[]
 for name,url in sources().items():
  path=OUT/name
  with requests.get(url, timeout=(20,90), stream=True) as r:
   r.raise_for_status(); chunks=[]; size=0
   for chunk in r.iter_content(262144):
    size += len(chunk)
    if size > 20*1024**2: raise ValueError('Metadata exceeds bounded download: '+name)
    chunks.append(chunk)
  data=b''.join(chunks); json.loads(data)
  path.write_bytes(data)
  entries.append({'file':name,'url':url,'bytes':size,'sha256':sha(path)})
  print('Fetched metadata:',name,size,flush=True)
 write(manifest_path,{'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sources':entries})

def main(do_fetch=False):
 if do_fetch: fetch()
 manifest=json.loads((OUT/'source_manifest.json').read_text())
 assert {s['file']:s['url'] for s in manifest['sources']} == sources()
 for s in manifest['sources']:
  p=OUT/s['file']; assert p.stat().st_size==s['bytes'] and sha(p)==s['sha256'],s['file']
 def load(n): return json.loads((OUT/n).read_text())
 def axis(n,col):
  with (ROOT/n).open(encoding='utf-8-sig',newline='') as f: return [r[col] for r in csv.DictReader(f)]
 genes=axis('gene_names.csv','gene_name'); targets=axis('pert_counts.csv','target_gene')
 assert len(genes)==len(set(genes))==18533 and len(targets)==len(set(targets))==300
 gs,ts=set(genes),set(targets)
 train={k:load(k+'_train.json') for k in ['wd','xd']}
 oldstats={k:json.loads((OLD/n).read_text())['stats'] for k,n in [('wd','PIE_wdataset_k562_data_stats.json'),('xd','PIE_replogle_xdataset_data_stats.json')]}
 for k in train: assert sha(OUT/(k+'_train.json'))==oldstats[k]['train_json_sha256'], 'Training split differs from checkpoint'
 meta={k:load(k+'_meta.json') for k in META_REVS}
 rows=[]
 for k,d in meta.items():
  rows.append({'dataset':d['dataset'],'vocabulary_genes':len(d['genes']), 'competition_genes':len(gs&set(d['genes'])), 'perturbation_identifiers':len(d['pert_to_id']), 'competition_target_identifiers':len(ts&set(d['pert_to_id']))})
 def split_summary(k):
  targets_in_split=set().union(*(set(v) for v in train[k].values()))
  per_ds={}
  for context,values in train[k].items():
   dataset=context.split('.',1)[0]
   z=per_ds.setdefault(dataset,{'contexts':0,'pairs':0,'target_identifiers':set()})
   z['contexts']+=1; z['pairs']+=len(values); z['target_identifiers'].update(ts&set(values))
  for z in per_ds.values():
   z['target_identifiers']=sorted(z['target_identifiers']); z['target_count']=len(z['target_identifiers'])
  return {'pairs':sum(map(len,train[k].values())), 'contexts':len(train[k]), 'distinct_perturbations':len(targets_in_split),'competition_target_count':len(ts&targets_in_split),'competition_targets':sorted(ts&targets_in_split),'by_dataset':per_ds}
 split_summaries={k:split_summary(k) for k in train}
 union=set().union(*(set(d['genes']) for k,d in meta.items() if k!='replogle_nadig_essential'))
 gene_text=set(json.loads((OLD/'gene_text_meta.json').read_text())['keys'])
 # Presence in a contributing dataset's axis is not a count of significant training labels.
 with (OUT/'target_coverage.csv').open('w',newline='',encoding='utf-8') as f:
  w=csv.writer(f); w.writerow(['target_gene','within_training','cross_training']+[r['dataset']+'_vocabulary' for r in rows])
  for t in targets:
   w.writerow([t,int(t in split_summaries['wd']['competition_targets']),int(t in split_summaries['xd']['competition_targets'])]+[int(t in meta[k]['pert_to_id']) for k in meta])
 head=git('rev-parse','HEAD'); remote=load('github_head.json')[0]['sha']
 files=git('ls-files').splitlines(); compared=[]; different=[]; line_endings=[]
 for name in files:
  original=ROOT/'reference/v18_pie/code'/name
  if original.is_file():
   a=(CODE/name).read_bytes(); b=original.read_bytes(); compared.append(name)
   if a!=b: line_endings.append(name)
   if a.replace(b'\r\n',b'\n')!=b.replace(b'\r\n',b'\n'): different.append(name)
 assert not different, 'Cloned source differs from prior downloaded source'
 revs={}
 for k,oldrev in [('wd',WD),('xd',XD)]:
  info=load(k+'_current.json')
  ck={s['rfilename']:s.get('lfs',{}) for s in info['siblings'] if s['rfilename'].endswith('.ckpt')}
  revs[k]={'snapshot_revision':info['sha'],'prior_revision':oldrev,'same_revision':info['sha']==oldrev,'checkpoints':ck}
 oldtrain=json.loads((OLD/'PIE_training_meta.json').read_text())
 assert meta['replogle_nadig_essential']==oldtrain
 assert split_summaries['wd']['competition_target_count']==0
 summary={
  'audit_date':'2026-10-18','script_sha256':sha(Path(__file__)),
  'source_manifest_sha256':sha(OUT/'source_manifest.json'),
  'local_code_head':head,'remote_code_head_at_snapshot':remote,
  'same_code_as_v18':head=='21a526add9ff3a81b4f962d9b16bbae570052c63',
  'remote_equals_clone':remote==head,'prior_source_files_compared':len(compared),'source_files_differing':different,'checkout_CRLF_only_differences':line_endings,
  'models':revs,'dataset_axes':rows,'training_split_coverage':split_summaries,
  'cross_contributing_axis_union':{'genes':len(union),'competition_genes':len(gs&union),'competition_genes_with_descriptions':len(gs&union&gene_text),'missing_competition_genes':sorted(gs-union)},
  'within_target_ensembl_map_present':bool(oldtrain.get('pert_ensembl')),
  'within_target_matching':'Exact gene symbols from canonical metadata and actual saved training split; no heuristic aliases',
  'v18_unused_output':'delta_p_pred = mean(log1p CP10k) shift; not rank or detection probability',
  'lfc_contract':'Log2 FC; trained on DE, finite, positive-FC labels only; no hidden output standardization in readout',
  'evidence_contract':'Donor responses from exact saved training split; query response excluded; novel context has no context response block',
  'models_loaded':False,'new_predictions':False,'submission_created':False,
  'limitations':['Gene-axis coverage does not prove tested/supervised coverage for every gene or successful generalization.',
   'Cross-dataset training includes Orion and VCC25; their reused public populations are not unseen-study tests.',
   'Exact-symbol overlap does not count biological synonyms as observed training evidence.',
   'This is a code/metadata audit, not a new biological performance experiment.']}
 summary['local_input_sha256']={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [ROOT/'gene_names.csv',ROOT/'pert_counts.csv',OLD/'PIE_training_meta.json',OLD/'gene_text_meta.json',OLD/'PIE_wdataset_k562_data_stats.json',OLD/'PIE_replogle_xdataset_data_stats.json',ROOT/'infer_v18_pie.py',ROOT/'v18_methods.py',ROOT/'V51_PLATEAU_DIAGNOSIS.md',ROOT/'experiments/v24_cd4_response/analysis/leaderboard_result.json']}
 write(OUT/'audit.json',summary)
 v24=json.loads((ROOT/'experiments/v24_cd4_response/analysis/leaderboard_result.json').read_text())
 # Rounded values transcribed from the user's October 18 leaderboard paste.
 leader={'pds':.878,'mse':.786,'nmae':.438,'fid':.172,'reach':.463,'jac':.088}
 gap={k:{'V24':v24['metrics'][k],'leader_user_rounded':v,'scaled_gap':v-v24['metrics'][k],'overall_gap_contribution':(v-v24['metrics'][k])/6} for k,v in leader.items()}
 total=sum(v['scaled_gap'] for v in gap.values())
 write(OUT/'leaderboard_gap.json',{'source':'User-provided leaderboard paste, Illumina AI / iCell, 87 submissions, overall display .4709; component scores rounded to 3 decimals. Not a fresh live API snapshot.','V24_entry':v24['entry_id'],'V24_overall':v24['overall'],'leader_display_overall':.4709,'metrics':gap,'expression_LFC_reach_share_of_component_gap':sum(gap[k]['scaled_gap'] for k in ['mse','nmae','reach'])/total,'interpretation':'Descriptive score decomposition, not an explanation of the leading team method or a performance forecast.'})
 print(json.dumps({k:summary[k] for k in ['local_code_head','remote_code_head_at_snapshot','prior_source_files_compared','dataset_axes','cross_contributing_axis_union'] if k!='cross_contributing_axis_union'},indent=2))
 print('Model revisions:',json.dumps({k:v['snapshot_revision'] for k,v in revs.items()}))
 print('Actual training split coverage:',json.dumps({k:{n:v[n] for n in ['pairs','competition_target_count']} for k,v in split_summaries.items()}))
 print('Cross-dataset contributing axes:',len(union),'competition genes:',len(gs&union))
 print('Audit passed; no weights loaded, predictions generated, or upload made.')

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--fetch',action='store_true'); main(p.parse_args().fetch)
