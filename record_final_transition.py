"""Record the user-supplied V90 publication receipt separately from frozen preparation evidence."""
from pathlib import Path
import json
R=Path(__file__).resolve().parent;A=R/'experiments/v90_interval_diagnostic/analysis'
def main():
 scores=dict(pds=.6542179228155035,mse=.09572270480857992,nmae=.08600057484213723,fid=-.02200589821017791,reach=.08399308779101305,jac=.014934286551255357)
 baseline=dict(pds=.6517009514165119,mse=.09657486167945828,nmae=.08527061361861961,fid=-.02260859548432163,reach=.07962563498187135,jac=.014729659086930079)
 overall=.15214377976638518;assert abs(sum(scores.values())/6-overall)<1e-14
 record=dict(entry='Fholtkq6MAMAifZSXZlR',model='Conservative interval projection v90 diagnostic',status='published',source='User-pasted CLI receipt; resumed after FileNotFoundError',rank_snapshot=551,overall=overall,scaled_metrics=scores,partition='val',panel='vcc2026-val-1',anchor_set='vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4',baseline_overall=.15088218754984492,overall_change=overall-.15088218754984492,metric_changes={k:scores[k]-baseline[k] for k in scores},local_file_sha256='e042294a250ba151d533edb360c036d36522568240b4986438860ec000f78f81',server_checksum_verified=False)
 p=A/'published_result.json'
 if p.exists():assert json.loads(p.read_text())==record
 else:p.write_text(json.dumps(record,indent=2),encoding='utf-8')
 rows='\n'.join(f'| {k} | {baseline[k]:.9f} | {scores[k]:.9f} | {scores[k]-baseline[k]:+.9f} |' for k in scores)
 (R/'V90_LEADERBOARD_RESULTS.md').write_text(f'''# V90 published result

Entry `{record['entry']}` is published on the validation panel, Overall **{overall:.12f}**, rank **551** at the supplied receipt. The upload resumed successfully from 482.5 MiB. V90 is the best published validation submission in this project, up **{record['overall_change']:.9f}** versus V66/V82 (about {100*record['overall_change']/record['baseline_overall']:.3f}% relative). This is a small single-panel gain, not proof of transfer to the final round.

| Scaled metric | V66/V82 | V90 | Change |
|---|---:|---:|---:|
{rows}

Panel: `val / vcc2026-val-1`. Anchors: `{record['anchor_set']}`. Five scaled metrics improved; expression accuracy declined slightly. Ranks across dates are not controlled comparisons. Original V88/V89 failed gates and V90 pre-upload experimental records remain unchanged. The publication receipt supersedes their historical best-model status, not their test results.

Final inputs use D/E/F and a different panel. Final scores are blind and cannot be compared directly to this validation score.
''',encoding='utf-8')
 head=b'''<!-- final-phase:begin -->\n## Final phase and V90 publication\n\nV90 is published at **0.15214377976638518**, up **0.00126159221654026** versus V66/V82 on the same validation panel and anchors. Five scaled metrics improved; expression accuracy slipped slightly. [Receipt and comparison](V90_LEADERBOARD_RESULTS.md). Older preparation blocks below retain their historical status.\n\nFinal data are under `final/`: new contexts D/E/F, 300 targets (129 overlap the validation panel), and the unchanged 18,533-gene axis. The project-local CLI was upgraded to 0.3.0. Final preparation must rebuild responses and use the final manifest and perturbation list. Final scores remain hidden; the latest final entry counts. See [final preparation notes](FINAL_ROUND.md).\n<!-- final-phase:end -->\n\n'''
 p=R/'README.md';prior=p.read_bytes()
 if not prior.startswith(b'<!-- final-phase:begin -->'):p.write_bytes(head+prior)
if __name__=='__main__':main()
