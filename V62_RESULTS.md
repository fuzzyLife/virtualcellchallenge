# V62 conserved-response experiment results

Fixed selected candidate: **assay_center50**. Public support gates: **False**. V24 remains the best published submission (0.142997855); no V62 leaderboard score is known.

## Mean-response screening

| Candidate | LFC squared error / V24 | LFC absolute error / V24 | Pooled-profile squared error / V24 |
| --- | ---: | ---: | ---: |
| centered50 | 0.993039 | 0.997613 | 0.994216 |
| assay50 | 0.976121 | 0.990582 | 0.968155 |
| assay_center50 | 0.973947 | 0.989413 | 0.965691 |

Ratios below one improve on the fold-specific V24 mechanism. These are equal-study aggregates over 44 reused public raw populations; they are not official scaled metrics or an independent biological reserve. All held-study source-ingredient poisoning tests passed. A separate pre-score fixture correction changed a floating round-trip assertion tolerance only; original code and plan are retained.

## Generated integer cells

### Development

| Arm | PDS | Fidelity | Reach | Jaccard | NMAE error | Expression error ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| v24 | 0.791667 | 0.581703 | 0.194389 | 0.070027 | 0.912601 | 37.471936 |
| assay_center50 | 0.750000 | 0.581614 | 0.230796 | 0.069467 | 0.904368 | 36.617947 |
| assay_center50_shuffled | 0.680556 | 0.579669 | 0.212597 | 0.067371 | 0.921841 | 36.649599 |

Checks: resolved=True, expression=True, pds=False, fid=True, reach=True, jac=False, nmae=True, three_improve=True, positive_gain=True, beats_shuffle_expression=True, beats_shuffle_pds=True, study_expression_guard=True.

Confirmation was not opened because development failed.

PDS, fidelity, reach and Jaccard above are raw higher-is-better diagnostics. NMAE and expression ratios are lower-is-better errors. Missing NMAE in an individual context remains missing; aggregate available NMAE does not imply complete metric coverage. Small 4/8-target discrimination panels are not the 300-target leaderboard. Expression is evaluated as a ratio of sums within each context with the jackknife cap, then equal-context/equal-study aggregation; it is not an average of per-cell MSEs. The generator, control inputs, row libraries and seeds are paired across arms. A shuffled source-target diagnostic tests target evidence; it is an inference corruption, not separately retrained biology.

## Interpretation

The fixed candidate did not establish a validated replacement for V24. Preserve the negative results and do not describe an experimental package as improved.

The raw development populations, public study labels and many source preprocessing choices have already been reused in earlier experiments. Confirmation holds out generated count cases and seeds within these same studies; it does not remove that reuse. No D3 outcomes are opened. Returning the same zero-response input to the unchanged emitter reproduces its old null behavior by construction, including its known excess null DE; this method does not solve that limitation.

## Reproduce

See [V62_CONSERVED_RESPONSE_STRATEGY.md](V62_CONSERVED_RESPONSE_STRATEGY.md), root v62_conserved_response.py and run_v62_counts.py. Scientific source, input and production-candidate hashes are frozen before count selection. No automatic upload or Git push.
