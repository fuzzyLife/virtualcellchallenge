"""Reuse byte-identical V65 baseline populations in V66, never candidate outcomes."""
from pathlib import Path
import shutil,time
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;src=ROOT/'experiments/v65_closed_response/analysis';dst=ROOT/'experiments/v66_assay_strength/analysis';assert not (src/'count_assessment.json').exists(),'Stage baselines before V66 starts'
folder=dst/'counts';folder.mkdir(exist_ok=True);records=[]
for p in sorted((src/'counts').glob('*.json')):
 if p.name.endswith('__panel.json') or not any(p.name.endswith('__'+a+'.json') for a in ('v24','v63')) or time.time()-p.stat().st_mtime<5:continue
 r=read(p);arr=p.with_suffix('.npz');assert sha(arr)==r['arrays_sha256'];out=folder/p.name;assert not out.exists();shutil.copyfile(arr,out.with_suffix('.npz'));assert sha(out.with_suffix('.npz'))==r['arrays_sha256']
 if 'reused_source' in r:r['original_upstream_reused_source']=r['reused_source']
 r.update(reused_source=str(p.relative_to(ROOT)),reused_receipt_sha256=sha(p));write(out,r);records.append(dict(file=p.name,source_sha256=sha(p),arrays_sha256=r['arrays_sha256']))
write(dst/'baseline_reuse.json',dict(files=records,populations=len(records),scope='Same immutable V24/V63 effects, controls/templates, library sizes, emitter and seeds. No candidate/shuffled arrays copied. All copied baseline count summaries and profile/DE arrays are byte-identical to the V65 experiment; these are shared computations, not independent biological evidence.'))
print('V66 reused',len(records),'byte-identical V65 baseline populations before count jobs began.')
