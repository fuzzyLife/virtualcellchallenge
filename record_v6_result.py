"""Record the published V6 entry, compare scores, and update its publication layer."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import json
import subprocess
import sys
import numpy as np
import pandas as pd
import build_submission as build
from record_v5_result import RAW

ROOT=build.ROOT
A=ROOT/'experiments/v6_fidelity_calibrated/analysis'
ENTRY='6MCjBWpPoCt5Sbw5CMu3'
REPORTED={'overall':-0.05679013942013697,'pds':0.25937425243103096,'mse':0.,
    'nmae':0.0761310724314523,'fid':-0.7302547485948779,
    'reach':0.08353078889187517,'jac':-0.02952220168030233}


def document_publication():
    """Idempotent publication overlay; safe to call after regenerating V6 methods."""
    result=json.loads((A/'leaderboard_result.json').read_text(encoding='utf-8'))
    scaled=pd.read_csv(A/'leaderboard_comparison.csv');raw=pd.read_csv(A/'raw_metric_comparison.csv')
    lines=['# V6 published leaderboard result','',
        f'Entry `{ENTRY}` is **published**, overall **{result["overall"]:.6f}**, rank {result["rank_at_user_report"]} at the user-reported check. '
        'V4 remains the best project submission at **-0.050674**. All five published project entries use the same validation panel and anchor set.','',
        '| Scaled metric | First | V2 | V4 | V5 | V6 | V6 minus V4 | Contribution to overall change |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in scaled.itertuples():
        contribution='' if row.metric=='overall' else f'{row.overall_contribution_change_from_v4:+.6f}'
        lines.append(f'| {row.metric} | {row.first:.6f} | {row.v2:.6f} | {row.v4:.6f} | {row.v5:.6f} | {row.v6:.6f} | {row.change_from_v4:+.6f} | {contribution} |')
    lines+=['',f'Overall changes by **{result["overall_change_from_v4"]:+.6f}** versus V4 and **{result["overall_change_from_v5"]:+.6f}** versus V5. '
        'Three scaled components improve (PDS, LFC accuracy and overlap), two worsen (fidelity and reach), and expression remains at its clipped zero score. '
        f'Fidelity contributes **{result["fidelity_overall_contribution_change_from_v4"]:+.6f}**; the other five together contribute **{result["other_components_overall_contribution_change_from_v4"]:+.6f}**. '
        'This is arithmetic attribution within the six-component mean, not causal attribution.','',
        '| Raw metric | V4 | V5 | V6 | V6 minus V4 | Improved versus V4 |','|---|---:|---:|---:|---:|---|']
    for row in raw.itertuples():
        lines.append(f'| {row.metric} | {row.v4:.6f} | {row.v5:.6f} | {row.v6:.6f} | {row.change_from_v4:+.6f} | {"yes" if row.improved_from_v4 else "no"} |')
    lines+=['',f'Raw expression error rises **{result["raw_expression_error_relative_change_from_v4"]:.2%}**; direction-fidelity yield falls **{-result["raw_fidelity_yield_relative_change_from_v4"]:.2%}**. '
        'A zero scaled expression score does not mean zero prediction error. Fidelity yield uses a response-call/reference budget; it is not simply the fraction of correct signs.','',
        '## What this changes','',
        'V6 retained the V4 stochastic generator and all A mean effects, blending 25% toward V5 only at existing same-sign B/C coordinates. Neural templates and full dual-moment generation were excluded. '
        'The positive public H1 direction/overlap results did not translate to a better aggregate anonymous-context leaderboard score. '
        'The source benchmark changes an H1 contrast gain, whereas production combines lineage estimators and applies an additional sign/zero mask; these are not identical response policies. '
        'Nine primary H1 targets and previously inspected source batches cannot establish anonymous-context transfer validity.','',
        'The control draw/emission seed also changed across all contexts. There is no published unchanged-V4 comparator using that seed, no context-level outcome breakdown, and no repeated-seed leaderboard uncertainty estimate. '
        'The aggregate difference cannot identify which context/target failed, establish statistical significance, or isolate the response blend from sampling variation. '
        'This result is not a test of neural perturbation prediction because the production neural weight was zero.','',
        'Keep V4 as the production reference. Future work should strengthen source-to-context validation and evaluate target/context-specific response reliability and seed stability before adopting another broad blend. '
        'V6 recovers much of V5\'s lost fidelity, but has not surpassed V4.','',
        '## Publication and provenance','',
        f'Submitted at `{result["submitted_utc"]}`. Upload took 13 minutes 29 seconds. Rank {result["rank_at_api_check"]} at the recorded API check; ranks change as other teams submit. '
        f'Partition `{result["partition"]}`, panel `{result["panel"]}`, anchors `{result["anchor_set"]}`.','',
        'The generation-time `submission_ready.json`, frozen model, predictions and original numerical validation records are preserved. '
        'Publication is recorded separately in `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv`. '
        'The package SHA-256 comes from the verified local readiness record; the status API does not independently expose or confirm the uploaded hash.','',
        '```cmd','cd /d L:\\promec\\Animesh\\Download\\Virtual',
        '.\\.competition_env\\Scripts\\python.exe -m vcc status '+ENTRY,
        '.\\.competition_env\\Scripts\\python.exe record_v6_result.py','```','',
        'These inspect the existing entry. [V6_SUBMISSION.md](V6_SUBMISSION.md) preserves generation, validation and the historical upload command.','']
    (ROOT/'V6_LEADERBOARD_RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
    path=ROOT/'V6_SUBMISSION.md';text=path.read_text(encoding='utf-8')
    text=text.replace('# V6: fidelity-guarded submission ready for upload','# V6: verified and published submission')
    text=text.replace('It has **not been uploaded or scored**. V4 remains the best published project submission; no local benchmark establishes a better hidden score.',
        'It is **published**. V4 remains the best project submission; the measured comparison is recorded above.')
    begin,end='<!-- V6 publication begin -->','<!-- V6 publication end -->'
    if begin in text:
        start=text.index(begin);stop=text.index(end,start)+len(end);text=text[:start]+text[stop:].lstrip('\n')
    block=begin+'\n\n'+f'Published as `{ENTRY}`, overall **{result["overall"]:.6f}**, rank {result["rank_at_user_report"]} at the user-reported check. '
    block+='V4 remains best at **-0.050674**. V6 improves PDS, LFC accuracy and overlap but loses fidelity and reach. '
    block+='See [V6_LEADERBOARD_RESULTS.md](V6_LEADERBOARD_RESULTS.md) for raw/scaled comparisons and interpretation.\n\n'+end+'\n\n'
    title,rest=text.split('\n',1);text=title+'\n\n'+block+rest.lstrip('\n')
    text=text.replace('## Submit\n','## Original upload command (historical)\n')
    text=text.replace('This uploads a new entry. Preserve its ID and six metrics to compare against V4 and V5 on the same panel/anchor set.',
        'This historical command created the published entry. Repeating it uploads a duplicate; use `vcc status '+ENTRY+'` to inspect the existing result.')
    path.write_text(text,encoding='utf-8')
    path=ROOT/'README.md';text=path.read_text(encoding='utf-8')
    text=text.replace('Four submissions (first, V2, V4 and V5) have been published.','Five submissions (first, V2, V4, V5 and V6) have been published.')
    text=text.replace('The newest upload-ready candidate is V6, documented in [V6_SUBMISSION.md](V6_SUBMISSION.md). It has not been uploaded or scored.',
        f'The latest published submission is V6 (`{ENTRY}`), with overall **{result["overall"]:.6f}**. V4 remains best. See [V6_LEADERBOARD_RESULTS.md](V6_LEADERBOARD_RESULTS.md) and [V6_SUBMISSION.md](V6_SUBMISSION.md).')
    start=text.index('## V6 fidelity-guarded submission\n');stop=text.index('\n## ',start+3)
    section=text[start:stop].replace('The package has not been uploaded or scored; V4 remains the best published project result.',
        f'Published as `{ENTRY}`, overall **{result["overall"]:.6f}** versus V4 **-0.050674**. Fidelity and reach decline despite gains in PDS, LFC accuracy and overlap. [V6_LEADERBOARD_RESULTS.md](V6_LEADERBOARD_RESULTS.md) records the API-confirmed raw/scaled metrics. V4 remains best.')
    text=text[:start]+section+text[stop:]
    text=text.replace('best of four published project submissions','best of five published project submissions').replace('best of the four published project submissions','best of the five published project submissions')
    text=text.replace('This is a proposal and diagnostic; no new submission has been generated.',
        'That proposal and diagnostic preceded the completed V6 experiment and publication documented above.')
    path.write_text(text,encoding='utf-8')
    path=A/'documentation_manifest.json';manifest=json.loads(path.read_text(encoding='utf-8'))
    manifest['root_code_sha256']['document_v6_submission.py']=build.sha256(ROOT/'document_v6_submission.py')
    manifest['root_code_sha256']['record_v6_result.py']=build.sha256(Path(__file__))
    manifest['document_sha256']=build.sha256(ROOT/'V6_SUBMISSION.md')
    manifest['publication_document_sha256']=build.sha256(ROOT/'V6_LEADERBOARD_RESULTS.md')
    for name in ['leaderboard_result.json','published_api_status.json','leaderboard_comparison.csv','raw_metric_comparison.csv']:
        manifest['artifact_sha256'][name]=build.sha256(A/name)
    manifest['scope']='Generation and verification plus separately API-confirmed published V6 results.'
    build.write_json(path,manifest)


def main(api_file=None):
    if api_file:
        data=Path(api_file).read_bytes();api=json.loads(data.decode('utf-16' if data.startswith((b'\xff\xfe',b'\xfe\xff')) else 'utf-8-sig'))
    else:
        process=subprocess.run([sys.executable,'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45)
        api=json.loads(process.stdout)
    assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
    for key,value in REPORTED.items():np.testing.assert_allclose(api['score_avg' if key=='overall' else 'score_'+key],value,rtol=0,atol=1e-12)
    previous={name:json.loads((ROOT/path/'leaderboard_result.json').read_text(encoding='utf-8')) for name,path in {
        'first':Path('analysis'),'v2':Path('experiments/v2_cd4_hybrid/analysis'),
        'v4':Path('experiments/v4_lineage_response/analysis'),'v5':Path('experiments/v5_atlas_calibrated/analysis')}.items()}
    for old in previous.values():assert (api['partition'],api['panel_id'],api['anchor_version'])==(old['partition'],old['panel'],old['anchor_set'])
    scaled=[]
    for key,value in REPORTED.items():
        row={'metric':key,**{name:old['overall'] if key=='overall' else old['metrics'][key] for name,old in previous.items()},'v6':value}
        row.update(change_from_v4=value-row['v4'],change_from_v5=value-row['v5'],overall_contribution_change_from_v4=None if key=='overall' else (value-row['v4'])/6);scaled.append(row)
    scaled=pd.DataFrame(scaled)
    delta=REPORTED['overall']-previous['v4']['overall'];fid=(REPORTED['fid']-previous['v4']['metrics']['fid'])/6
    np.testing.assert_allclose(scaled.iloc[1:].overall_contribution_change_from_v4.sum(),delta,atol=1e-12,rtol=0)
    np.testing.assert_allclose(np.mean([REPORTED[k] for k in REPORTED if k!='overall']),REPORTED['overall'],atol=1e-12,rtol=0)
    raw=pd.DataFrame([{'metric':key,'v4':previous['v4']['raw_metrics'][key],'v5':previous['v5']['raw_metrics'][key],'v6':api[key],
        'change_from_v4':api[key]-previous['v4']['raw_metrics'][key],'higher_is_better':higher,
        'improved_from_v4':api[key]>previous['v4']['raw_metrics'][key] if higher else api[key]<previous['v4']['raw_metrics'][key]} for key,higher in RAW.items()])
    ready=json.loads((A/'submission_ready.json').read_text(encoding='utf-8'))
    result={'source':'User-reported CLI publication, confirmed by read-only vcc --json status.','reported_date':'2026-10-07','recorded_utc':datetime.now(timezone.utc).isoformat(),
        'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published','submitted_utc':api['submission_date'],
        'rank_at_user_report':1052,'rank_at_api_check':api['rank'],'overall':REPORTED['overall'],'metrics':{k:v for k,v in REPORTED.items() if k!='overall'},
        'raw_metrics':{k:api[k] for k in RAW},'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],
        'package_sha256':ready['sha256'],'package_hash_scope':'Verified local readiness record; status API does not expose the uploaded package hash.',
        'upload_duration_seconds':809,'overall_change_from_v4':delta,'overall_change_from_v5':REPORTED['overall']-previous['v5']['overall'],
        'fidelity_overall_contribution_change_from_v4':fid,'other_components_overall_contribution_change_from_v4':delta-fid,
        'raw_expression_error_relative_change_from_v4':api['expr_mse_unbiased_capped_norm']/previous['v4']['raw_metrics']['expr_mse_unbiased_capped_norm']-1,
        'raw_fidelity_yield_relative_change_from_v4':api['de_wilcoxon_direction_fidelity_yield_raw']/previous['v4']['raw_metrics']['de_wilcoxon_direction_fidelity_yield_raw']-1,
        'interpretation':'PDS, LFC accuracy and overlap improve; fidelity and reach decline; expression is still clipped to zero. V4 remains best. Source gains did not improve aggregate leaderboard performance; seed changes and aggregate-only outcomes prevent causal/context attribution.',
        'readiness_record_scope':'submission_ready.json remains the immutable generation-time snapshot.','recorder_sha256':build.sha256(Path(__file__))}
    build.write_json(A/'published_api_status.json',api);build.write_json(A/'leaderboard_result.json',result)
    scaled.to_csv(A/'leaderboard_comparison.csv',index=False);raw.to_csv(A/'raw_metric_comparison.csv',index=False)
    document_publication();print(scaled.to_string(index=False));print(raw.to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--api-file',type=Path,help='Previously captured vcc --json status output');args=parser.parse_args();main(args.api_file)
