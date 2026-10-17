# V53: corrected joint-model feature contract

V53 preserves the complete V52 experiment and reuses its exact paired-view mean/rank labels. It changes one preprocessing rule: an imputed numerical feature with fitting-set standard deviation <=1e-8 and its associated missingness indicator are fixed to zero in both fitting and inference. Other channels retain the original transformation. This conservatively drops an associated missingness indicator even if it varied; available-source counts and other explicit support features remain.

The issue was found during implementation review after V52 H1 held-study results and source-validation histories had been inspected. A held-out query can have more donor studies than a fitting recipient, filling donor slots that were completely missing during fitting. The original transform allowed those new values to activate untrained network columns. A small independent fixture reproduces this and verifies that the corrected transformation is invariant to changes in disabled columns. V53 is a disclosed corrective follow-up, not a blind new biological hypothesis.

All V52 source code, plan, models and results remain unchanged. V53 has its own frozen plan and directories. The architecture, training settings, labels, folds, losses, seeds, shuffled controls and pilot checks are identical; all outer neural models and linear heads are retrained. No threshold is changed in response to the held-study scores. The default fitting call also trains an all-source joint checkpoint.

Linked label metadata record both plans and the original NPZ locations. Independent V52 raw-cell label audits remain applicable only because every linked label byte and source metadata field is verified unchanged. The corrected feature normalization and inference receive separate checks. Missing values remain missing evidence, not biological zero votes.

This correction does not solve statistical domain shift, missing response biology or count reconstruction. In particular a query can have different source counts and control distributions from the fitting population even after unused columns are disabled. All mean/rank performance measures remain public diagnostics, not VCC scaled scores or a submission release.

Use the same pinned scientific Python 3.12 runtime and prerequisites as [V52](V52_JOINT_STRATEGY.md):

~~~bat
SCIENCE_PYTHON test_v53_joint_guard.py
SCIENCE_PYTHON run_v53_joint_guard.py freeze
SCIENCE_PYTHON run_v53_joint_guard.py train
SCIENCE_PYTHON run_v53_joint_guard.py evaluate
SCIENCE_PYTHON run_v53_joint_guard.py audit --kind boundaries
SCIENCE_PYTHON run_v53_joint_guard.py audit --kind artifacts
~~~

Freeze only once. `run_v53_finish.py` verifies linked-label audits and completes evaluation/auditing after the training receipts are available. The API remains `predict_v52_joint_response.predict`; first call `run_v53_joint_guard.install()` to select the corrected checkpoint, preprocessing and verified label bank. The original V52 API remains available without that installation.
