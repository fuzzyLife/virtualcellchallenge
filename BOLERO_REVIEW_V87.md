# BOLERO review and V87 conditional response test

The user supplied `geneRegPred.pdf` and the `bolero` clone. Before implementation I read the clone's AGENTS/README, conditional factor generation and collapse/forward code, BorzoiLoRA setup/forward/gene head, installation instructions, released model configurations, and the relevant paper results/methods. Rendered PDF pages 4, 10, 11 and 44 were inspected. Extracted text/renderings are local intermediates under `tmp/pdfs/v87/`.

## What the new paper establishes

[Liu et al., *Predicting cell-state-specific gene regulation from DNA across species*, bioRxiv, posted October 20, 2026](https://doi.org/10.64898/2026.10.19.760026) describes a DNA sequence model with cell-state-dependent low-rank weight updates. The backbone is frozen; a conditioning network generates adapter weights. The authors evaluate genomic signals across cell states and species. Their AP-1 example conditions on a measured motif-activity proxy to predict activity-associated chromatin differences. These are author-reported preprint results, not evidence of performance on VCC's 300 knockdowns.

Relevant local PDF sections: architecture/training pp.20-22; held-state comparisons pp.5-7; AP-1 example pp.10-12 and methods p.28; availability pp.34-35. Figure 1 and Extended Data Figure 1 distinguish across-feature from across-sample performance. Strong across-gene reconstruction alone cannot establish accurate intervention responses.

The PDF specifies 10.7 million cells, 33 datasets and six mammals; the cloned README/agent guide mentions 10.8 million/36. The supplied manuscript is the scientific reference for this review.

## Release checked, not assumed unavailable

- Local clone: `0b328b65b04f1927a07001e32750e42c8b422013`.
- [Public models](https://huggingface.co/arcinstitute/bolero-models), inspected revision `47023a28a4fd647cbc4e99004bf76032dc15c9f6`: twelve real checkpoint files, approximately 4.82-4.85 GB each. This is measured file metadata; an older registry's 12.2 GB estimate is not used as the download size.
- [Public atlas](https://huggingface.co/datasets/arcinstitute/bolero-data), revision `bbd1ef1c282007ababf5bace499b1151aea8be87`: 2,163 listed files, including cell embeddings and chromVAR artifacts. This is not a new CRISPRi response bank.
- Downloaded and retained only two small gene-model configurations and metadata. Neither upstream checkpoint nor atlas-wide cell matrices were downloaded, and no upstream model inference was executed.
- Actual gene configurations specify `gene_count`, conditional adapters, a separate Flashzoi base checkpoint, dataset records and embedding/condition artifacts. AP-1 gene configuration also references the cross-dataset chromVAR activity file. Arbitrary RNA PCA coordinates are not interchangeable with those saved model embeddings.
- Documented installation targets Linux with NVIDIA CUDA/FlashAttention. This Windows project has an 8 GB GPU. We have not measured an OOM and do not claim CPU/offloaded inference is impossible; an unverified platform port would be separate work.
- Recent clone commits mostly prepare release paths/dependencies/docs and fix region readers. They do not demonstrate a better VCC submission. `inspect_bolero_v87.py` records commit history, hashes, public revisions, checkpoint sizes and configurations.

Code is [MIT licensed](https://github.com/liuhlab/bolero). The preprint is CC-BY-NC. Do not infer checkpoint/data licensing solely from the code license; check artifact terms before redistribution.

## Gap between the paper and our task

1. VCC needs a full RNA count population after silencing any of 300 genes. A DNA variant substitution is not the same intervention as gene knockdown.
2. Basal RNA expression does not supply a calibrated post-knockdown TF activity score. Changing a gene's RNA input is not automatically an identified intervention on the sequence model.
3. Direct application needs correctly aligned DNA windows, gene masks, compatible embedding/condition encoders and a trained mapping from each knockdown to the changed conditioning variables. We have not established that mapping.
4. Our useful supervision consists of a few assay-confounded public intervention contexts. More parameters cannot manufacture independent context-response evidence.

A direct, properly aligned BOLERO regulatory-feature experiment remains possible future work. V87 instead isolates its conditional-weight idea using the actual response supervision already validated here. It is an independent compact adaptation, not a reproduction of the full model, DNA regulatory grammar, AP-1 experiment, or reported benchmarks.

## V87 implementation and frozen experiment

See `experiments/v87_conditional_adapter/analysis/plan.json` for exact parameters, hashes and gates, frozen before fitting/outcomes. Root code:

- `v87_conditional_adapter.py`: source exclusion, complete observed axis, low-rank basis, fixed backbone, conditional factors, inference and V66 identity guard.
- `run_v87_experiment.py`: freeze, fit, recipient prediction, public mean scoring and release decision.
- `test_v87_adapter.py`: zero-adapter equality, explicit per-example collapsed-matrix forward equivalence, gradients/frozen base and monotone-rank tests.
- `audit_v87_experiment.py`: excluded-study poisoning, basis orthogonality/target partition, checkpoint-selection replay and adapter matrix equations.
- `inspect_bolero_v87.py`: public checkpoint/configuration and local clone/paper provenance.

Only complete measured nonpanel outputs are fitted. Missing response rows are omitted. A rank-32 basis is learned from 240 fitting targets, with 60 target identities held for source-only epoch selection. Directed training pairs come from distinct study groups. An ordinary residual response MLP is trained first. Its weights are frozen, and low-rank factors for both layers are generated from source/query control-rank sketches. The layer starts exactly at the frozen backbone because B is zero-initialized. This is distinct from V14 linear context slopes, V30 STRING regressions, V35 attention/gene decoding and V84 linear essential-target transfer.

Three fixed outputs are retained: `conditional25`, `constant25` with conditioner inputs zero, and `conditional25_shuffled` with source response identities shifted before independently fitting the whole pipeline. All use the same bounded quarter correction and existing V66 identity guard. No target outcome from the recipient study selects epochs or hyperparameters. The post-proposal guard preserves profile margins, but closure can change non-own effect factors outside the initially proposed axis; it is not a promise of unchanged final panel factors.

The sole eligible primary must lower equal-study public profile error by at least 0.5%, lose no more than 0.005 PDS, keep every study within 3% of baseline error, and beat both controls. These 655 public pairs have been reused extensively and are development data. A passing primary would need the unchanged D3 count gates and second seed before packaging. Failure stops before generating cells. No upload is automatic.

## Reproduce

Use the scientific runtime recorded in `requirements_v85_v86.txt` (PyTorch, NumPy/SciPy/scikit-learn); V87 uses the same installed versions. Run from the root with `C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe -X utf8`:

```text
inspect_bolero_v87.py
test_v87_adapter.py
run_v87_experiment.py freeze
run_v87_experiment.py fit --outer H1
run_v87_experiment.py fit --outer Orion
run_v87_experiment.py fit --outer CD4
run_v87_experiment.py fit --outer H1 --shuffled
run_v87_experiment.py fit --outer Orion --shuffled
run_v87_experiment.py fit --outer CD4 --shuffled
audit_v87_experiment.py models
run_v87_experiment.py prepare
run_v87_experiment.py means
```

Each line is a Python script invocation using that runtime. Freeze, inspection, fits and predictions intentionally refuse to overwrite existing records. To reproduce from scratch, use a separate checkout/output directory with the documented data and frozen source bytes. Final quantitative results and any remaining limitations belong in `V87_RESULTS.md`.
