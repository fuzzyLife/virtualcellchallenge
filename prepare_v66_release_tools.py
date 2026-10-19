"""Create reviewed V66 production adapters; do not modify predecessor releases."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for name in ('fit','build','package','validate'):
 src=ROOT/f'{name}_v65_submission.py';dst=ROOT/f'{name}_v66_submission.py';assert not dst.exists()
 s=src.read_text(encoding='utf-8-sig').replace('v65_closed_response','v66_assay_strength').replace('V65','V66').replace('v65','v66').replace('available-axis compositional assay response','strength-matched original-assay response')
 if name=='fit':s=s.replace('m.check();assert not','m.check();m.old.check();assert not').replace("m.check();p=read(A/","m.check();m.old.check();p=read(A/")
 dst.write_text(s,encoding='utf-8')
p=ROOT/'verify_v66_effects.py';assert not p.exists();s=(ROOT/'verify_v65_effects.py').read_text(encoding='utf-8-sig').replace('v65_closed_response','v66_assay_strength').replace('V65','V66').replace('v65','v66');s=s.replace('m.check();model=','m.check();m.old.check();model=');p.write_text(s,encoding='utf-8')
s=(ROOT/'audit_v62_boundaries.py').read_text(encoding='utf-8-sig').replace('import v62_conserved_response as m','import v66_assay_strength as m\nm.data=m.old.data; m.bank=m.old.bank; m.GROUP=m.old.GROUP; m.FAMILIES=m.old.FAMILIES').replace('V62','V66');(ROOT/'audit_v66_boundaries.py').write_text(s,encoding='utf-8')
s=(ROOT/'audit_v65_results.py').read_text(encoding='utf-8-sig').replace('import v65_closed_response as m','import v66_assay_strength as m').replace('V65','V66');(ROOT/'audit_v66_results.py').write_text(s,encoding='utf-8')
print('V66 production and independent audit adapters created.')
