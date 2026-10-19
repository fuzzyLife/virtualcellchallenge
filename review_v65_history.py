"""Rebuild the complete saved publication trajectory; ranks are historical snapshots."""
from pathlib import Path
import json,csv,hashlib,re,math
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v65_closed_response/analysis'
METHODS={1:'Public CRISPR transfer',2:'Donor-calibrated CD4 hybrid',4:'Lineage-guided response coverage',5:'Atlas-calibrated dual moments',6:'Fidelity-guarded lineage response',7:'V4 fallback resampling',9:'Activation-preserving profile calibration',10:'Pooled intervention population transport',11:'Consensus coverage and dispersion',13:'Assay-aware knockdown transport',16:'Response residual dispersion',20:'Stack-calibrated response ensemble',24:'Context-matched donor response',36:'Retrieval-contrast correction',50:'Empirical response prior',55:'Fixed V24/V36 half-blend',61:'Pooled state-response transfer',63:'Identity-guarded assay response'}
METRICS=('pds','mse','nmae','fid','reach','jac')

def main():
 A.mkdir(parents=True,exist_ok=True);rows=[];inputs={}
 paths=list((ROOT/'experiments').glob('*/analysis/leaderboard_result.json'))+list((ROOT/'experiments/v50_response_prior/analysis').glob('published_result.json'))
 for p in paths:
  r=json.loads(p.read_text(encoding='utf-8-sig'))
  if 'scaled_scores' in r:r=dict(r,overall=r['scaled_scores']['overall'],metrics={k:r['scaled_scores'][k] for k in METRICS},entry_id=r['entry'],rank_at_user_report=r['rank_at_reporting'])
  match=re.match(r'v(\d+)_',p.parents[1].name)
  if not match:continue
  v=int(match[1]);assert v in METHODS,v
  assert math.isclose(sum(r['metrics'][k] for k in METRICS)/6,r['overall'],abs_tol=1e-12)
  rows.append(dict(version=v,method=METHODS[v],overall=r['overall'],rank_snapshot=r.get('rank_at_user_report',r.get('rank_at_api_check')),entry_id=r['entry_id'],panel=r['panel'],anchor_set=r['anchor_set'],**r['metrics']))
  inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 # Earliest receipt, if it was saved outside experiments.
 for p in [ROOT/'analysis/leaderboard_result.json']:
  if p.exists():
   r=json.loads(p.read_text(encoding='utf-8-sig'));rows.append(dict(version=1,method=METHODS[1],overall=r['overall'],rank_snapshot=r.get('rank_at_user_report'),entry_id=r['entry_id'],panel=r['panel'],anchor_set=r['anchor_set'],**r['metrics']));inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 rows.sort(key=lambda r:r['version']);assert len({r['version'] for r in rows})==len(rows)
 best=-float('inf');prev=None
 for r in rows:
  assert r['panel']=='vcc2026-val-1' and r['anchor_set']=='vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4'
  r['change_from_previous_submission']=None if prev is None else r['overall']-prev
  r['change_from_previous_best']=None if best==-float('inf') else r['overall']-best
  r['new_project_best']=r['overall']>best
  best=max(best,r['overall']);prev=r['overall']
 with (A/'published_history.csv').open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (A/'published_history.json').write_text(json.dumps(dict(records=rows,inputs_sha256=inputs,scope='All saved canonical published receipts, one per version. Unsubmitted V numbers are experiments, not leaderboard results. Historical rank changes are confounded by other teams. Successive versions often change multiple components; association is not causal attribution.'),indent=2)+'\n',encoding='utf-8')
 for r in rows:print(f"V{r['version']:2d} {r['overall']:+.9f} {'NEW BEST' if r['new_project_best'] else ''} {r['method']}")
if __name__=='__main__':main()
