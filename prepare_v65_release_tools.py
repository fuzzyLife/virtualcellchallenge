"""Create V65 adapters from the reviewed immutable V63 writer and integrity checker."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for name in ('build','package','validate'):
 source=ROOT/f'{name}_v63_submission.py';dest=ROOT/f'{name}_v65_submission.py'
 assert not dest.exists(),dest
 text=source.read_text(encoding='utf-8-sig').replace('v63_identity_guard','v65_closed_response').replace('V63','V65').replace('v63','v65')
 text=text.replace("'recommended_over_V24'","'recommended_over_V63'")
 dest.write_text(text,encoding='utf-8')
print('Created V65 writer, official package adapter and independent payload validator; prior files unchanged.')
