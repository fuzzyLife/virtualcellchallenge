"""Correct V14 freshness claims using older V7 split-file formats.

Does not change frozen identities, numerical code, model fits or outcomes.
Its findings supersede freshness statements in the initial plan/split audit.
"""
from pathlib import Path
import json
from v14_context_methods import ROOT,A
import build_submission as build


def main():
    plan=json.loads((A/'plan.json').read_text());used=set(plan['previous_count_targets']);inputs={}
    future=ROOT/'analysis/post_v11_research/future_experiment_plan.json'
    used.update(json.loads(future.read_text())['excluded_previous_count_targets']);inputs[str(future.relative_to(ROOT))]=build.sha256(future)
    for file in (ROOT/'experiments/v7_multisource_adaptive/analysis').glob('cell_split_*.json'):
        d=json.loads(file.read_text());used.update(d['targets']);inputs[str(file.relative_to(ROOT))]=build.sha256(file)
    old_claim=set(plan['confirmation_targets']);reused=sorted(old_claim&used);fresh=sorted(old_claim-used)
    record={'supersedes':'Initial plan.json/split_audit.json and initial V14_STRATEGY statements claiming all reserved targets were never count-tested.',
            'reason':'The metadata inventory initially read validation_plan.json files but omitted older V7 cell_split_*.json formats. Detected by an independent audit before V14 count outcomes or summary selection were read.',
            'numerical_effect':'None. Original identities, controls, seeds, candidates and methods remain unchanged. Partial K562/CD4 model fits preserved; screen restarted with hash-verified caches.',
            'plan_sha256':build.sha256(A/'plan.json'),'inputs':inputs,'complete_previous_target_exclusions':sorted(used),
            'reserved_targets_previously_count_tested':reused,'reserved_targets_globally_unused':fresh,
            'development_reuse_corrected':{c['family']:[t for t in c['development_targets'] if t in used] for c in plan['contexts']},
            'interpretation':'V14 is a prospective algorithm comparison on retrospective public biology. Reserved count verification uses disjoint target names from development but partly reuses historical biological contrasts; it is not independent blind biological confirmation.',
            'script_sha256':build.sha256(Path(__file__))}
    path=A/'metadata_correction.json'
    if path.exists():assert json.loads(path.read_text())==record
    else:build.write_json(path,record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
