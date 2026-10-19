"""Apply the frozen V69 confirmation guards with explicit unavailable handling.

The original assessor attempts numeric comparisons on None when a whole public
study lacks eligible DE labels. Unavailable results remain unavailable and fail
their guard; they never become zero or get dropped from the equal-study mean.
"""
import numpy as np
import run_v69_experiment as runner
from v57_pie import read,write,sha
import v69_count_bridge as m


def resolved(a,b,k):return a.get(k) is not None and b.get(k) is not None


def main():
    m.check();assert not (m.A/'confirmation_assessment.json').exists()
    s=runner.summary('confirmation');av=s['summary'];ss=s['studies'];a=read(m.A/'selection.json')['selected'];v=av[a];b=av['v66']
    ratios={arm:{g:ss[arm][g]['profile_square_error']/ss['v66'][g]['profile_square_error'] for g in ss[arm]} for arm in av if arm!='v66'}
    checks=dict(profile=float(np.mean(list(ratios[a].values())))<=.99,worst_study=max(ratios[a].values())<=1.03)
    for k,loss in [('pds',.005),('fidelity',.002),('reach',.005),('jaccard',.0005)]:checks[k]=resolved(v,b,k) and v[k]>=b[k]-loss
    checks['nmae']=resolved(v,b,'nmae') and v['nmae']<=b['nmae']*1.01
    gains=sum(resolved(v,b,k) and v[k]>b[k] for k in ('pds','fidelity','reach','jaccard'))
    gains+=int(resolved(v,b,'nmae') and v['nmae']<b['nmae'])+int(v['profile_square_error']<b['profile_square_error'])
    checks['three_gains']=gains>=3
    if a!='mean_bridge':checks['beats_shuffle']=v['profile_square_error']<av['shuffled_shape']['profile_square_error']
    matched=read(m.A/'matched_control_assessment.json');checks['beats_matched_control']=bool(matched['beats_matched_control_profile'])
    unresolved=[k for k in ('pds','fidelity','reach','jaccard','nmae') if not resolved(v,b,k)]
    write(m.A/'confirmation_assessment.json',dict(s,ratios=ratios,selected=a,checks=checks,release_supported=all(checks.values()),unresolved_metrics=unresolved,
        assessment_note='Before confirmation, diagnostic exposed an unavailable NMAE study. This separate assessor implements the original thresholds with undefined comparisons false, preserving the original evaluator. Additional matched-mixture shuffle guard can only make support stricter.',
        assessor_sha256=sha(m.ROOT/'assess_v69_confirmation.py'),frozen_evaluator_sha256=sha(m.ROOT/'run_v69_experiment.py')))
    print('V69 confirmation',a,'support',all(checks.values()),'checks',checks,'unresolved',unresolved,flush=True)

if __name__=='__main__':main()
