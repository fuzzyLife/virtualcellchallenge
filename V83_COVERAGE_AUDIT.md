# V83 coverage audit: new context names did not mean new target observations

The attempt to directly add HepG2/Jurkat/RPE1 effects was rejected before any recipient treatment evaluation. The pinned Arc essential-screen DE tables have **zero target-name overlap with the300 competition targets**. Source vocabularies contain1253 HepG2,1439 Jurkat and1453 RPE1 perturbations, with6642 output genes. We checked complete identifier columns, not only Parquet row-group previews or a potentially broken filter.

The local `reference/means_hepg2.parquet` and `means_jurkat.parquet` each have zero rows. They are empty for this target restriction; this is **not evidence of an ingestion bug**. An older V17 metadata audit had already recorded zero overlap at the same publisher revision. That earlier result should have been carried into subsequent planning. The initial direct-pooling idea repeated an already known incompatibility; it was corrected before evaluating donor responses. The publisher's biological or filtering reason for the absence is not established here, and we do not infer that every competition target is nonessential.

Source release: `arcinstitute/PIE_replogle_nadig_essential`, revision`91f86154aa1761a2dd0d885b979b9169f3b9983f`. Full Parquet downloads were verified against publisher LFS SHA-256. The published README specifies CP10k-normalized raw-count means; inspecting sums confirms the released axis sums to10000. There are target-specific reference cohorts, so a single identical-reference assertion in the first parser was invalid. The original parser and correction are retained in`analysis/engineering_history`. No recipient outcome was used in that correction.

V83 consequently has no new direct-response predictor to evaluate or upload. Its source bank correctly abstains on every competition target. This preserves the failed hypothesis rather than manufacturing biological zero votes from absent targets.

The source data remain useful for a different task: **learn how a measured response changes across contexts using other shared interventions**. V84 implements that hypothesis on638 shared essential targets, with513 fitting and125 held targets, before applying a frozen operator to existing K562 evidence for our targets. This is an independent response-transfer prototype, not an implementation of a published foundation model.

## Reserved donor preparation

The previously deferred Marson2025 D3 Rest download was resumed and fully verified against MD5`e8c066b72f12caea34d26985f49474ea`. Only annotations were inspected before freezing the14 metadata-qualified targets. The source has423747 cells; extracting10000 NTCs and those targets gives18463 aligned cells,18533 competition columns,69347099 nonzero integer counts. Missing genes are recorded as unmeasured, not asserted biological zeros.

Fit/reference controls were split within lane into5006/4994 disjoint cells. Every target has400 distinct actual cells and400 distinct lane-matched NTC templates per technical seed. The14 targets are ADNP,AKAP1,MAP2K2,MARCKSL1,MYRF,PCGF1,PRDX3,SLC39A1,SUPT7L,USP28,ZC3HAV1,ZNF25,ZNF280C,ZNF514. Two seeds may overlap cells and are not independent biological replication.

All CD4 response data, including pooled donor summaries, must be excluded from donor prediction. The donor belongs to a previously used public study and therefore is a project-unseen donor, not a new study or final-test surrogate. Treatment values were extracted but no effect or model score was opened while selecting the cohort.

## Reproduce

- `collect_v83_sources.py`: retrieve the three pinned new-context DE tables; publisher SHA-256 checks.
- `collect_v84_sources.py`: also retrieve their paired K562 essential-screen table.
- `extract_v83_donor.py`: freeze cohort from annotations first, then verify and align raw integer counts. It deliberately refuses to overwrite an existing extraction.
- `v83_context_expansion.py sources`: direct-response bank; coverage rejection is recorded explicitly.

Run numerical/data scripts using the project's pinned scientific Python environment. Immutable plans, download receipts, source axes, metadata cohorts and code histories are retained under`experiments/v83_context_expansion/analysis`. Large data/arrays are ignored by Git; code and compact evidence are committed locally.
