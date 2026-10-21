"""Record the user's V82 receipt without changing immutable preparation evidence."""
import argparse,csv,hashlib,json
from decimal import Decimal as D
from pathlib import Path
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'analysis/v82_published'
BASE=ROOT/'experiments/v66_assay_strength/analysis/leaderboard_result.json';READY=ROOT/'experiments/v82_v66_recovery/analysis/submission_ready.json'
VALUES=dict(overall='0.1508821875498449',pds='0.6517009514165119',mse='0.09657486167945828',nmae='0.08527061361861947',fid='-0.02260859548432163',reach='0.07962563498187135',jac='0.014729659086930079')
METRICS=('pds','mse','nmae','fid','reach','jac')
BEGIN='<!-- v82-published:begin -->';END='<!-- v82-published:end -->'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encoded(v):return (json.dumps(v,indent=2,allow_nan=False)+'\n').encode('utf-8')
def main(check=False):
 base=read(BASE);ready=read(READY);assert ready['byte_identical_to_V66'] and not ready['new_predictor'];assert ready['sha256']==base['package_sha256'];assert (ROOT/'experiments/v82_v66_recovery/prediction.vcc').stat().st_size==ready['bytes']
 panel=dict(partition='val',panel='vcc2026-val-1',anchor_set='vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4');assert all(base[k]==v for k,v in panel.items())
 rows=[]
 for k in ('overall',)+METRICS:
  old=D(str(base['overall'] if k=='overall' else base['metrics'][k]));new=D(VALUES[k]);delta=new-old;assert abs(delta)<=D('1e-15');rows.append(dict(metric=k,v66=str(old),v82=str(new),difference=str(delta)))
 assert abs(sum(D(VALUES[k]) for k in METRICS)/6-D(VALUES['overall']))<D('1e-15')
 result=dict(source='Exact published CLI receipt supplied by user; no API query.',entry_id='sMObaQ7bhzac3RfNINSv',model_name='V66 recovery after V69 - identical predictions',final_status='published',rank_at_user_report=550,overall=float(VALUES['overall']),metrics={k:float(VALUES[k]) for k in METRICS},exact_reported_values=VALUES,**panel,scoring_job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-smobaq7bhzac3rfninsv',upload_duration_seconds=965,local_package_bytes=ready['bytes'],local_package_sha256=ready['sha256'],artifact_link_evidence='User command identifies the audited local recovery file; server receipt provides no checksum.',same_predictions_as_v66=True,new_predictor=False,all_reported_scores_match_v66_within_1e_15=True,comparison_entry_id=base['entry_id'],previous_v66_rank_snapshot=base['rank_at_user_report'],rank_change_is_not_a_score_change=True,submitted_utc=None,raw_metrics=None)
 report='''# V82 published: V66 recovery reproduced

The user-supplied CLI receipt reports **published** entry `sMObaQ7bhzac3RfNINSv`, model **V66 recovery after V69 - identical predictions**, overall **0.1508821875498449**, rank **550** at the time of the receipt.

All six reported scores and overall match V66 within 1e-15 on the same validation panel and anchors. The minuscule differences in overall and NMAE are floating-point precision, not evidence of a predictive change. This is the same prediction package, not a newly learned model. V66's old rank 518 and V82's rank 550 are different leaderboard snapshots; the score was reproduced while competitors' standings can change.

|Scaled score (higher is better)|V66|V82|Difference|
|---|---:|---:|---:|
'''+ '\n'.join(f"|{r['metric']}|{r['v66']}|{r['v82']}|{r['difference']}|" for r in rows)+'''

The recovery restores the best published project predictions after V69's 0.09526284481907588. It does not establish an improvement beyond V66, statistical significance, a new biological result or final-test performance. V70-V81 remain rejected under their original public criteria.

The local package SHA-256 is `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`, 4,073,666,560 bytes. Prior preparation verified complete byte equality with V66 and official container validity. This recorder checks the saved evidence and current size; it does not rehash the large package or query the server. The receipt contains no server checksum or exact submission timestamp, so neither is invented.

The scoring job is `projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-smobaq7bhzac3rfninsv`. Upload duration was 965 seconds. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

Reproduce with `python record_submission_result_82.py`, then `python record_submission_result_82.py --check`. Exact decimal strings and comparison CSV are retained in `analysis/v82_published`; publication input/output hashes are recorded there. Original preparation files keep their historical not-uploaded state; this publication record supersedes that status. No upload, API request, training or GitHub push is performed by the recorder.
'''
 block='''<!-- v82-published:begin -->
## V82 published: V66 score reproduced

Recovery entry `sMObaQ7bhzac3RfNINSv` is **published**, overall **0.1508821875498449**, rank **550** at the user-supplied receipt. All six component scores match V66 within floating-point precision on the same panel and anchors. **These are the same V66 predictions, not a new model or an improvement beyond V66.** The old rank 518 and current 550 are different leaderboard snapshots; the score stayed the same.

[Published receipt, exact score comparison and provenance](V82_LEADERBOARD_RESULTS.md). Older preparation blocks retain their historical status and are superseded by this publication record.
<!-- v82-published:end -->'''
 inputs={p.relative_to(ROOT).as_posix():sha(p) for p in [BASE,READY,Path(__file__)]};paths=[OUT/'leaderboard_result.json',OUT/'leaderboard_comparison.csv',ROOT/'V82_LEADERBOARD_RESULTS.md']
 if check:
  audit=read(OUT/'publication_audit.json');assert audit['inputs_sha256']==inputs;assert read(paths[0])==result
  for name,h in audit['outputs_sha256'].items():assert sha(ROOT/name)==h,name
  assert paths[2].read_bytes()==report.encode('utf-8')
  with paths[1].open(encoding='utf-8',newline='') as f:assert list(csv.DictReader(f))==rows
  for p in [ROOT/'README.md',ROOT/'V82_RECOVERY.md']:assert p.read_bytes().decode('utf-8-sig').count(BEGIN)==1 and block in p.read_bytes().decode('utf-8-sig')
  print('V82 publication verified: same panel, exact decimals, score equivalence and metadata hashes.');return
 OUT.mkdir(parents=True,exist_ok=True);paths[0].write_bytes(encoded(result))
 with paths[1].open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 paths[2].write_bytes(report.encode('utf-8'))
 for p in [ROOT/'README.md',ROOT/'V82_RECOVERY.md']:
  text=p.read_bytes().decode('utf-8-sig');assert text.count(BEGIN)==text.count(END)
  if BEGIN in text:
   start=text.index(BEGIN);stop=text.index(END)+len(END);text=text[:start]+block+text[stop:]
  else:text=block+'\n\n'+text
  p.write_bytes(text.encode('utf-8'))
 (OUT/'publication_audit.json').write_bytes(encoded(dict(inputs_sha256=inputs,outputs_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in paths},all_scores_match_v66_within_1e_15=True,same_panel_and_anchors=True,original_preparation_evidence_unchanged=True,large_binary_rehashed=False)))
 print('Recorded V82 published; all scores reproduced within floating-point precision.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
