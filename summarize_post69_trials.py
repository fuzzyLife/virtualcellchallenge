"""Summarize completed post-regression trials without inventing hosted scores."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/post_v69_review/analysis'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 A.mkdir(parents=True,exist_ok=True);rows=[];proof={};count_records=0;rank_coordinates=0;shape_tables=[]
 for v in range(70,82):
  folder=next((ROOT/'experiments').glob(f'v{v}_*'))/'analysis';selection=read(folder/'selection.json');evidence=[folder/'selection.json'];row=dict(version=v,selected=selection.get('selected'),new_predictor_uploaded=False)
  if v<=76:
   p=folder/'development_mean_results.json';d=read(p);row['development_profile_error_ratios']=d['summary'];evidence.append(p)
   if (folder/'broader_assessment.json').exists():
    p=folder/'broader_assessment.json';b=read(p);row['broad']={k:b[k] for k in ['profile_ratio','pds_change','target_context_pairs','checks','passes']};evidence.append(p)
   if v==76:
    p=folder/'count_audit.json';audit=read(p);assert audit['passed'];count_records+=audit['comparison_records'];rank_coordinates+=audit['manual_rank_coordinates'];evidence.append(p)
  elif v==77:
   p=folder/'development_results.json';d=read(p);row['development']=d['summary'];evidence.append(p)
  else:
   p=folder/'development_0_assessment.json';d=read(p);row['summary']=d['summary'];row['checks']=d['checks'];evidence.append(p)
   p=folder/'development_0_independent_audit.json';audit=read(p);assert audit['passed'];count_records+=audit['comparison_records'];rank_coordinates+=audit['manual_rank_coordinates'];evidence.append(p)
   suffix='_opposite' if v==81 else '_shuffled';base=d['summary']['v66']
   for arm,s in d['summary'].items():
    if arm=='v66' or arm.endswith(suffix):continue
    shape_tables.append(f"|V{v} {arm}|{s['fidelity']-base['fidelity']:+.6f}|{s['jaccard']-base['jaccard']:+.6f}|{s['reach']-base['reach']:+.6f}|{all(d['checks'][arm].values())}|")
  row['release_qualified']=False
  if v in [75,76]:
   assert row['broad']['passes'] is False
  else:assert selection['selected'] is None,(v,'A qualifying selection requires separate confirmation')
  row['evidence']={str(p.relative_to(ROOT)):sha(p) for p in evidence};proof.update(row['evidence']);rows.append(row)
 recovery=read(ROOT/'experiments/v82_v66_recovery/analysis/submission_ready.json');assert recovery['byte_identical_to_V66'] and not recovery['new_predictor'] and not recovery['uploaded']
 result=dict(experiments=rows,new_hypothesis_families=len(rows),fully_audited_count_comparisons=count_records,manual_tied_rank_coordinates=rank_coordinates,biological_cases_reused=True,new_qualifying_predictor=False,recovery=recovery,evidence_sha256=proof,scope='Public raw statistics differ from hosted anchor-scaled metrics. Repeated populations/configurations are not independent biological samples. No new hosted score or final-test gain is predicted.')
 (A/'comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 text='''# After V69: tests, failures and a transparent recovery release

V69 fell from V66's published 0.150882188 to 0.095262845. Its intended mean effects were frozen, but the emitted count profiles were not protected; discrimination and expression accuracy both regressed. This run tested twelve hypothesis families (V70-V81) rather than changing one mixture and uploading immediately. No new predictor met its frozen public requirements.

## What the experiments showed

| Versions | Change | Result |
|---|---|---|
|V70-V71|Batch-matched raw-count teacher; remove identity cap; center or confidence-project response|No material expression gain.|
|V72|Recover small batches by aggregating batch-offset count likelihoods|Source coverage improves; predictive error does not.|
|V73|Null-centered count MAP and separate assay/bulk inversion|No material gain; strongest confidence-filtered correction loses to shuffled evidence.|
|V74|Amplify downstream responses in assay space|Expression error worsens by 2.5%-8.3%.|
|V75|Shrink downstream responses while retaining own knockdown|Broader expression error improves 2.36%, PDS falls 0.04781. Rejected.|
|V76|Dilute the whole normalized response, including own knockdown|Broader expression improves 1.73%, PDS falls 0.02682; generated count PDS/reach/Jaccard also worsen.|
|V77|Remove the shared response across targets|No discrimination gain on 55 metadata-selected cases. Remaining600-pair confirmation not run.|
|V78|Restore control variance at exact V66 gene and cell totals|Small changes give weak gains below thresholds; large changes damage DE precision/overlap.|
|V79|Signed-underdispersion-aware molecule thinning plus innovations|Original prospective criteria fail.|
|V80|Empirical zero/detection hurdle and positive count calibration|Original prospective criteria fail.|
|V81|Opposite-direction residual contraction, with expansion controls|Original prospective criteria fail.|

The new batch-teacher initially dropped many tiny matched batches; V72 fixes that coverage defect and records the larger coverage. That defect was in the new teacher, not an established bug in V66. Correcting it did not improve prediction. V75 also exposed that a three-development/one-confirmation split left singleton PDS panels in H1/HCT116. V76 states its reused all-case PDS scope; V78-V81 reserve two cases/context, leaving at least two for confirmation. Failed/partial jobs are retained and are not counted as completed confirmation.

## Exact means are necessary but not sufficient

V78-V81 retain every realized V66 row library and every gene total, hence bitwise-identical raw profiles and PDS. The unchanged integer projector also limits cell-normalized mean drift. Cell distributions can still change Wilcoxon statistics, DE means and the sampling correction used by expression accuracy. The following are equal-study *raw public* metric differences from paired V66, not hosted scaled scores:

|Arm|FID change|JAC change|REACH change|All development checks|
|---|---:|---:|---:|---|
'''+ '\n'.join(shape_tables)+f'''

Independent audits replayed **{count_records} completed baseline/candidate/control count comparisons** and **{rank_coordinates} raw tied-rank coordinates**, plus all stored BH arrays, four DE metrics and PDS. These are repeated configurations of reused public cases, not {count_records} independent biological samples. V76 uses a different primary cohort from V78-V81, so differences across their tables are not comparative estimates. Plans, selections, failures, source exclusions and exact input/code hashes are preserved. The individual model criteria were frozen before their outcomes; none was relaxed to promote a failed candidate. Remaining-target tests and second seeds were conditional on success and were not run for failed development arms. V81's control-only guard was specified but not claimed as passed without execution.

## What is ready to upload

`experiments/v82_v66_recovery/prediction.vcc` is a **byte-identical copy of V66**, our best published model. It is 4,073,666,560 bytes, SHA-256 `{recovery['sha256']}`. The official container validator passes; equality to the previously fully audited payload protects all 360,000 cells and 900 groups. This is a recovery release after V69, **not a new model, new prediction experiment or improvement beyond V66**. The prior V66 public guards were also imperfect and remain recorded as failed. Its published score is evidence for recovery, not a guarantee of current rank or final-test performance. Re-uploading uses an allowance and is only useful if restoring the current leaderboard entry is desired.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v82_v66_recovery\\prediction.vcc -m "V66 recovery after V69 - identical predictions" --wait
```

Nothing was uploaded automatically. All original sealed artifacts were preserved. Do not upload a V70-V81 file: no production file was generated for those rejected models.

## Reproduction and next evidence

`summarize_post69_trials.py` regenerates this report from saved result/audit records. `prepare_v82_recovery.py` refuses to copy the champion unless all twelve pilots have completed and failed release, verifies both copies by SHA-256 and calls the official validator. Run it with the official CLI environment. Run numerical experiments/audits with the separately pinned scientific environment; see `MARGIN_VARIANCE_STRATEGY.md` and the per-version strategies/results for commands. Large arrays remain local and ignored by Git; code and compact provenance are committed locally.

The remaining gap is transferable intervention information and validation transfer, not finding another arbitrary temperature. A useful next test needs genuinely unused response data: for example a new CD4 donor with **all CD4 evidence excluded from fitting**, a fixed target/cohort selection, and a mean/decoder crossed comparison fixed before seeing its treatment values. That donor was not downloaded or tested during this run. Reusing the current cases for further tuning cannot establish that a tiny public improvement will transfer to the challenge. This is a recommendation for new evidence, not a claimed successful model or an implemented future experiment.
'''
 (ROOT/'POST_V69_RESULTS.md').write_text(text,encoding='utf-8');print('Report complete',len(rows),'families',count_records,'comparisons',rank_coordinates,'raw-rank coordinates')
if __name__=='__main__':main()
