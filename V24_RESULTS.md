# V24 measured results

**Passed the locked raw-condition reserve.** Locked candidate: `cdf50_de`. The new package passed official preparation and independent full-payload validation.

V16 remains the best published project result, overall0.1416138568460195. V24 has not been uploaded and has no leaderboard score.

## Development

36 cases and eight arms =288 full-count populations. The28 previously studied family-holdout cases have400 cells; eight Rest D1 cases have200 cells and overlap public summaries, so they are training-like calibration. All Stim48hr responses are excluded from teacher and baseline fitting. HCT116 NMAE is unavailable. Baseline is the V16 mechanism refit without Stim48hr, not the exact published package.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | shape_wasserstein |
|---|---|---|---|---|---|---|---|
| amplitude | 0.79434524 | 0.5920563 | 0.37887019 | 0.034428262 | 0.89018428 | 51.204575 | 0.40155777 |
| cdf25 | 0.79434524 | 0.59454231 | 0.3883288 | 0.035429403 | 0.88505305 | 51.194077 | 0.41375171 |
| cdf50 | 0.79434524 | 0.5969615 | 0.38120867 | 0.036214402 | 0.88138783 | 51.180299 | 0.41085798 |
| cdf50_de | 0.79434524 | 0.59771145 | 0.38142231 | 0.036828339 | 0.87979551 | 51.180297 | 0.40177811 |
| shuffled25 | 0.78988095 | 0.58921397 | 0.36423669 | 0.034066642 | 0.89358625 | 51.207956 | 0.40166429 |
| shuffled50 | 0.78988095 | 0.59119806 | 0.38390263 | 0.034162913 | 0.89260598 | 51.20795 | 0.40720587 |
| shuffled50_de | 0.79434524 | 0.59111625 | 0.38444143 | 0.03381933 | 0.89153524 | 51.208061 | 0.40491511 |
| v16 | 0.79434524 | 0.58976165 | 0.38066987 | 0.034345974 | 0.89318488 | 51.207973 | 0.40684736 |

The development-only winner is selected by predeclared metric-improvement ranks. Three finalists pass available guards; the direction-only expert wins. Its raw changes relative to the same-input comparator are:

```json
{
  "pds": 0.0,
  "fidelity": 0.007949800626501663,
  "reach": 0.0007524322781264255,
  "jaccard": 0.0024823646387138257,
  "nmae": -0.013389376914744807,
  "relative_expression_worsening": -0.02869526163297365
}
```

## Once-only fresh condition reserve

Eight metadata-selected targets,400 distinct real treated cells each,400 unique templates and disjoint lane-stratified reference controls. Same public study; condition response rows are withheld, but upstream nuisance preprocessing may overlap. No independent-study or hidden-anchor accuracy claim. All experimental arms are retained for diagnosis; only the locked development winner is assessed for release.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | shape_wasserstein |
|---|---|---|---|---|---|---|---|
| amplitude | 0.71428571 | 0.52782609 | 0.0033687978 | 0.17868782 | 0.97680443 | 1.0002934 | 0.37131158 |
| cdf25 | 0.71428571 | 0.5299906 | 0.0029302013 | 0.17913898 | 0.97498329 | 0.99625345 | 0.37087986 |
| cdf50 | 0.73214286 | 0.53026242 | 0.0033687978 | 0.18095934 | 0.97109044 | 0.98973028 | 0.36798423 |
| cdf50_de | 0.73214286 | 0.53262269 | 0.0033687978 | 0.17882389 | 0.97344275 | 0.98970877 | 0.36578285 |
| shuffled25 | 0.69642857 | 0.52405272 | 0.0035084782 | 0.17844064 | 0.97842462 | 1.0028867 | 0.36122306 |
| shuffled50 | 0.67857143 | 0.52757646 | 0.0035084782 | 0.17775238 | 0.97688058 | 1.0029393 | 0.37213264 |
| shuffled50_de | 0.67857143 | 0.52326072 | 0.0021114672 | 0.18036985 | 0.97703378 | 1.0032587 | 0.35881666 |
| v16 | 0.69642857 | 0.52483183 | 0.0035068745 | 0.17962174 | 0.97772279 | 1.0028448 | 0.35908081 |

Locked candidate checks and changes:

```json
{
  "candidate": "cdf50_de",
  "pass": 1,
  "checks": {
    "pds": 1,
    "fidelity": 1,
    "reach": 1,
    "jaccard": 1,
    "expression": 1,
    "nmae": 1,
    "three_metrics_improve": 1,
    "DE_metric_improves": 1,
    "family_fidelity": 1,
    "family_reach": 1,
    "beats_alpha_matched_target_shuffle": 1
  },
  "changes": {
    "pds": 0.03571428571428581,
    "fidelity": 0.007790860895049878,
    "reach": -0.0001380766731644001,
    "jaccard": -0.0007978518693042003,
    "nmae": -0.004280049042397427,
    "relative_expression_worsening": -0.013098748154846414
  },
  "improvements": {
    "pds": 1,
    "fidelity": 1,
    "reach": 0,
    "jaccard": 0,
    "nmae": 1,
    "expression": 1
  },
  "complete_primary_metrics": 1,
  "missing_nmae_panels": []
}
```

## Independent checks

Two focused model tests pass (withheld-response removal, bounded correction, zero-compatibility identity and absent-evidence abstention). Source/code hashes, every frozen original cell identity, all count-profile cache hashes, condition exclusions and independent case/formula rechecks pass. The auditor replays both V16 and selected predictions in five full-axis cases (10 generated populations) and checks no-effect identity to the inherited emitter. This does not establish calibrated null DE.

```json
{
  "passed": 1,
  "frozen_case_count": 44,
  "count_populations": 352,
  "cache_hashes_checked": 44,
  "full_axis_replayed_cases": 5,
  "full_axis_two_arm_replays": 10,
  "no_effect_emitter_V16_identity_cases": 5,
  "condition_response_exclusion_verified": 1,
  "teacher_source_rows_and_conditions_verified": 400,
  "independent_raw_input_audit_sha256": "4528e3eba3856c79b77136e696ef893913dc8b5d44af578976d08f2a4051f219",
  "selected_candidate": "cdf50_de",
  "reserve_release_guards_pass": 1,
  "plan_sha256": "1549076a625d6659a05a728ac011e8fc8919500d90c755a1fac746b371e64766",
  "selection_sha256": "155b99c9b86b2491fc394e5f2a380f535b303bae5b8d7d6b237ccb5847631ab2",
  "confirmation_sha256": "f1ab5bcf0985555296b9d0a55ec1bbec31c6f4fa1c5e56accd983be08095d917",
  "auditor_sha256": "a5a54d8b2f98dd6757f3ce554579522116a23f7afe09722ad0615566375888b0",
  "limitations": "Rest development overlaps public summaries; condition reserve is existing study and upstream nuisance preprocessing may share information. Development HCT116 NMAE missing. Baseline refit excludes48hr; published production V16 arrays are separately retained. Null identity does not prove null-DE calibration."
}
```

## Production artifact

`experiments/v24_cd4_response/prediction.vcc`: 4,073,656,320 bytes; SHA256 `42cd3906e17aa2e8e546efea44f1e726a38e577605f211ad96f1564362a20cc6`. All900 groups,360,000 cells and18,533 genes pass official requirements. 207 groups receive new effects; unchanged groups reuse verified V16 count blocks. All three public conditions are used only in final refitting. Every production coefficient is independently reconstructed; every packaged CSR value/index/pointer is compared and12 full groups are regenerated exactly. Original control IDs, guides, unique sampling and libraries are checked for every prediction.

No automatic upload or Git push occurred. See V24_SUBMISSION.md for the exact CMD upload command.
