# V76 strategy and validation

V75 lowered expression error but lost discrimination, and its singleton confirmation panels could not define PDS. V76 therefore uses complete predicted profiles: p=(1-gain)*control+gain*V66, including effective direct-target suppression, with gain0.75/0.50. Direct-own-factor0.30/0.55 arms isolate the dosage change. Source V66 profiles use factor0.1 once; candidate effects encode their complete own response and candidate emission disables the legacy own override. This is a mean-profile regularizer, not an estimate of observed responder cells.

Input-proportion norm-matched cyclic target controls preserve the recipient's genuine own probability before reclosure. Nonown changes are bounded relative to the closed V66 reference; encoded effects are clipped +/-3. Bounds and rounding can change exact mixture equality after decoding. The existing sampler and control libraries stay fixed.

Development initially ranked candidates by mean error. A stricter joint PDS assessor was specified and applied before new-model confirmation. The original mean-only preselection is preserved in `mse_only_preselection.json`; the active `selection.json` is the joint choice. Mixture050 failed the additional PDS constraint; mixture075 passed with3.05% lower profile error and unchanged three-target development PDS. This is development tightening, not promotion after a confirmation failure.

Count PDS uses all4/4/6/8/8 canonical unique cases (30 total), avoiding singleton panels. Only remaining1/1/3/5/5 cases (15 total) contribute primary profile/DE summaries. All source V66 baselines are regenerated, old arrays checked, and raw p values retained for independent rank/BH/all-four-DE/PDS replay. A repeat seed and production require primary and broad gates; no new-model result is called blind because all public biology has been reused. D3 stays unopened.

Root scripts, frozen plans, explicit own-effect production adapters and compact evidence reproduce every stage. No automatic upload or push. See `V76_RESULTS.md` for the final verdict; V66 remains the best published until a new leaderboard receipt is available.
