# FAQ audit — 2026-10-09

Reviewed all 50 answers on the [official FAQ](https://virtualcellchallenge.org/faq), including the supplied screenshot, and compared them with our current code. The important addition is final-round readiness. No model or submission was changed.

## Official guidance relevant to this project

| Topic | FAQ guidance |
|---|---|
| Final round | D/E/F replace A/B/C; the perturbation panel also changes. Controls arrive October 22, 2026; final entry due November 5. |
| Biology | Anonymous cell lines; CRISPRi with 10x Flex; no separate competition training set. |
| External resources | External/proprietary data are allowed with appropriate rights. Arc model use remains subject to applicable licenses. |
| Entries | Two results per day; one entry considered for prizes. |
| Scores | Validation and final scores cannot be directly compared. |
| Density | Stored entries determine the limit; explicit zeros count too. |

Source: [FAQ](https://virtualcellchallenge.org/faq). The release/deadline answer does not specify a cutoff time or timezone. A different final panel does not imply that every target is different: overlap remains unknown until release.

## Consequences from inspecting our code

**V11 is a frozen validation artifact.** `manifest.json` specifies the validation partition and A/B/C. `fit_v11_submission.py` requires the source-cache target axis to equal the current specification. `build_v11_submission.py` requires the stored model targets and genes to match, verifies original control-file hashes, and preserves the frozen output. These protections correctly prevent substituting a different panel or new controls into an existing model.

The lightweight audit verified current model/source target-axis equality, context order, 18 frozen numerical dependencies, 17 original generation snapshots, and every file in the prior evaluation audit manifest. None of those files changed. The current source cache is indexed by the present 300 targets; final coverage cannot be inferred from that cache alone. The model can still inform a new run, but renaming contexts on its predictions will not create a final-round prediction.

Before the final run, load the new round into a separate input directory, preserve the current files, read its manifest/gene/target axes, recompute baseline moments from its controls, rebuild target-keyed source responses and record missing-target coverage. Generate a new artifact with original labels and a separate provenance record. Full source releases or a target-queryable response store would reduce dependence on the present panel. Public-source tests should examine joint context and target generalization before selecting such a method; the current held-family benchmark principally measures transfer for targets with source responses.

**Lineage estimates remain hypotheses.** The current response bank pools measured source effects and gates them using baseline expression. That supports using control similarity as uncertain evidence. It does not verify the identity of A/B/C or establish how D/E/F will respond. Prefer conditioning on measured basal profiles to treating a proposed tissue assignment as known truth.

**No context swap was found.** `build_submission.load_control` checks the file's observed context label against the requested label. The frozen builder verifies the original input hashes and matching gene/target order. Archived independent validation reports complete payload/axis/provenance agreement. This audit reads that record; it does not re-read the multi-gigabyte package. Keep these checks in the next pipeline because a plausible score alone cannot establish label correctness.

## Two FAQ wording details to handle carefully

**Strict per-cell pruning below 1 CPM is ineffective under the stated depth cap.** For positive integer count `x >= 1` and cell depth `L <= 1,000,000`, `CPM = 1,000,000*x/L >= 1`. That threshold removes no positive entries. A cutoff on average gene expression is a different operation. V11's archived validation counted 2,375,850,660 stored entries, about 50.02% of the 4.75-billion limit, so density pruning is unnecessary for this artifact. Removing explicit sparse zeros remains useful; the existing validator already requires every stored count to be positive.

**Absence of an upper clamp at 1 does not imply an infinite attainable score.** With fixed reference anchors, bounded raw metrics still have finite scaled ceilings. The detailed [official specification](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md) gives those limits. Likewise, a zero clamped expression score can include performance worse than the baseline. Use the exact metric definitions when interpreting the FAQ's simplified wording. Our previous [evaluation audit](EVALUATION_CRITERIA_AUDIT.md) documents the scaling, NMAE eligibility and backend gaps.

## Reproduction

From the project root in Windows Command Prompt:

```cmd
.\.competition_env\Scripts\python.exe audit_faq.py
```

Use `--show` to display the extracted answers for manual review. The script reads public JavaScript as text without executing it, verifies its pinned SHA-256, checks 16 FAQ statements, checks current model axes/provenance, and verifies the CPM inequality. It stores source URLs/hashes, findings and code-evidence line numbers in `analysis/faq_audit/audit.json`; full FAQ answers are not republished. A changed upstream asset requires a fresh review before updating the pin.

This is an audit and a prospective final-round preparation recommendation. It neither revises the frozen research plan nor evaluates a new candidate. The remaining priority is still the official-compatible, adequately powered validation layer described in the evaluation audit, followed by isolated emitter and response experiments.
