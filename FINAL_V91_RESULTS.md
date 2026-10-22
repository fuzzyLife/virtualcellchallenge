# Final V91 submission ready

The first final-round candidate is officially packaged with **VCC CLI 0.3.0** and independently validated. It applies the method of our best published validation submission (V90, Overall 0.15214377976638518) to the actual final contexts and targets. **Final biological accuracy and ranking remain unknown.** No final upload occurred automatically.

## Upload file

`experiments/final_v91/prediction.vcc`

- Bytes: 3,929,180,160.
- SHA256: `7595e024a9a4ceac08937e0b673dd056be58bb005afc8fbfe2b7f35f7565bb88`.
- Contexts: **D/E/F**; panel: **vcc2026-test-1**.
- Cells: 360,000; genes: 18,533; target/context groups: 900; cells/group: 400.
- Stored nonzeros: 2,283,563,644; maximum per-cell count total: 53,865.
- Exact final perturbation-list fingerprint is embedded in the official package. It differs from the validation fingerprint.

## What changed for the final round

There are 171 new targets and 129 targets shared with validation. We reused 638 genuine public target/context response rows bitwise and added 849 newly extracted public rows. 298 of 300 final targets have at least one usable source response; among new targets, 169 of 171 do. Full family coverage and missing masks are in `target_coverage.csv` and `source_audit.json`. A target with one measured source is not equivalent to a target measured in several lineages. Unsupported targets: **ANG, DPY30**. Their downstream response stays at the control-based baseline; the fixed on-target factor remains an assumption. All targets still have the required prediction cells.

The source estimator choices, correction strengths, 80% identity guards, on-target factor and residual10 count emitter remain fixed from V90. Source response pooling/centering, the donor-supported CD4 teacher, rank axes and finite-update covariance priors were refit for the actual final target panel. New D/E/F controls supply basal profiles, original cell templates and library totals. We did not relabel A/B/C predictions or predict only the shared targets.

The prior uses its original finite 100-update budget and normal-moment interval approximation. Its convergence and interval coverage are not established. Source-only fitting does not validate transfer to anonymous final contexts.

## Verification

- All nine archived V24/V66/V90 response arrays replayed **bitwise** under the working Python3.14 runtime; three archived full 400-cell V90 populations also regenerated exactly.
- Input axes, labels, every control raw count, guide sizes and source hashes checked.
- Newly acquired Orion expression offsets and full-source library totals verified; NTC counts match the original source selections exactly. H1 controls were included once across public splits.
- Final source masks/own coordinates and exact reused response rows audited.
- Fitted source likelihoods checked against the original numerical implementation, mixture weights checked, all final own-transfer coordinates zero, both identity guards passed.
- Twelve final full-population preflight groups passed, followed by all 360,000 cells.
- Official prep checked final genes, targets, D/E/F and the 400-cell requirement. Lossless staging compared every count/index/pointer.
- Independent validation compared the entire packaged payload, all original selected cell identities/library totals and twelve exactly regenerated full final groups.

These establish correct construction and reproducibility, not unseen treatment accuracy. The final data contain no perturbed-cell measurements to score biological predictions locally.

## Basal compatibility

Pairwise finite-gene Spearman rank correlations, without added pseudocounts:

| Final context | Nearest available external control | Spearman |
|---|---|---:|
| D | CD4_Stim48hr | 0.850588 |
| E | H1 | 0.776744 |
| F | H1 | 0.802042 |

These descriptive similarities do not identify the cell lines or prove response transfer.

## Submit from CMD

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\final_v91\prediction.vcc -m "V90 source-refit D E F final v91" --wait
```

The final phase is blind: a successful receipt should identify the **test** partition, with no public score or rank until announcement. Only the team's **latest final submission** is considered. Upload before **November 5, 2026, 23:59 UTC** (November 6, 00:59 Europe/Oslo). See the [official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/).

If interrupted, use `vcc submit --resume <entry-id> --wait` with this project's Python, preserving the same file. Do not start a duplicate entry while one is in progress.

See [FINAL_ROUND.md](FINAL_ROUND.md) and the root `run_final_pipeline.py` plus its stage scripts for complete reproduction. Large controls/source arrays/predictions remain local and ignored by Git; compact provenance and code are committed.
