# V70 work log

1. Reviewed V66 response transfer, count generation, public evaluation and the V69 regression. Kept old sealed releases unchanged.
2. Implemented a raw-count response teacher with batch matching, library-size offsets, uncertainty weighting and capped study votes. Frozen protocol: `analysis/plan.json`.
3. Passed synthetic depth-only, known-effect, uncertainty and zero-count checks using `test_v70_count_teacher.py`.
4. Fitted six public sources: H1 (25 eligible targets), HCT116 (275), HEK293T (287), CD4 Rest D1 (23), Rest D2 (34), Stim48 D2 (34). These are repeated source datasets, not new biological validation.
5. Prepared exact whole-study exclusion audit and fixed paired count-evaluation code before development results. Selection has no failed-guard fallback.
6. Development/confirmation and release decisions pending. D3 remains unopened; no upload or push.
