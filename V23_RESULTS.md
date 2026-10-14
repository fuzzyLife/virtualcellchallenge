# V23 measured results and release decision

V23 passed its development checks but **failed the once-only new raw-cell reserve**. No V23 prediction package was generated. V16 remains the best published submission (overall0.1416138568460195).

Development compared28 paired400-cell cases at six arms (168 generated populations). Mean-preserving conditional variation changed raw fidelity by+0.00434615, reach by+0.00439498 and significant-set Jaccard by+0.000484023. NMAE changed by-0.000020207 and paired expression error by-0.1687%. Centered response-projection Wasserstein distance fell by0.00215199. These are exploratory public raw metrics, not anchored leaderboard scores; HCT116 NMAE is unavailable. Full count hashes match both the exact V22 conditional populations and all28 old V16 comparators.

Confirmation used eight metadata-selected donor2 resting CD4 targets,400 unique real treatment cells per target,400 unique templates, and disjoint lane-stratified NTC reference halves. All CD4 family responses were excluded from fitting. An independent audit checked all25,281 retained cell IDs, labels, lanes, libraries and quality flags, plus every aligned count in128 original source rows. The author's full-file MD5 passed.

| Raw metric | V16 comparator | V23 mean preserved | Change |
|---|---:|---:|---:|
| PDS |0.696429|0.696429|0|
| Fidelity |0.369179|0.348700|-0.0204792|
| Reach |0.000409|0.000521|+0.000112158|
| Jaccard |0.188562|0.182984|-0.00557761|
| NMAE (lower better) |0.979397|0.979397|-0.000000766|
| Centered projection distance (lower better) |0.424034|0.399596|-0.0244381|

Fidelity and Jaccard violate their frozen drop limits, and fewer than three primary metrics improve. All primary metrics are available in this reserve. Improved projection shape does not compensate for poorer gene-level direction and significance. The original uncorrected conditional arm was diagnostic only and cannot be selected after opening confirmation.

The experiment identifies an important limitation: conditional population shape cannot reliably fix a poorly transferred mean response. V24 therefore tests a new mean-response expert with its own prospective plan and fresh condition reserve, rather than relaxing V23's criteria.

Six no-effect cases exactly match the inherited V16 emitter; approximately1,487 significant calls per null population still show that identity alone does not establish calibrated null DE. Numerical code, arrays, selection and confirmation receipts remain frozen. See V23_STRATEGY.md for exact runtime, limits, source caveats and reproduction.
