"""Create separate V70 adapters from reviewed V66 streaming/packaging validators."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for name in ('build','package','validate'):
 destination=ROOT/f'{name}_v70_submission.py';assert not destination.exists()
 source=(ROOT/f'{name}_v66_submission.py').read_text(encoding='utf-8-sig')
 source=source.replace('v66_assay_strength','v70_count_teacher').replace('V66','V70').replace('v66','v70')
 source=source.replace('recommended_over_V63','recommended_over_V66')
 destination.write_bytes(source.encode('utf-8'))
print('Separate V70 adapters created; not executed.')
