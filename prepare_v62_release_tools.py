"""Prepare the V62 production builder without changing earlier implementation files."""
from pathlib import Path
root=Path(__file__).resolve().parent
for name in ['build_v61_submission.py','package_v61_submission.py','validate_v61_submission.py']:
    target=root/name.replace('v61','v62')
    if target.exists():raise FileExistsError(target)
    text=(root/name).read_text(encoding='utf-8-sig').replace('v61','v62').replace('V61','V62').replace('v62_pooled_state_response','v62_conserved_response')
    text=text.replace("'methods_sha256':build.sha256(ROOT/'v62_pooled_state_response.py')", "'methods_sha256':build.sha256(ROOT/'v62_conserved_response.py')")
    text=text.replace("'inherited_transfer_methods_sha256':build.sha256(ROOT/'v59_state_composition.py')", "'inherited_transfer_methods_sha256':build.sha256(ROOT/'v15_ensemble_methods.py')")
    text=text.replace("'candidate':'V62 fixed joint50'", "'candidate':model['candidate']")
    text=text.replace("'experimental':True", "'experimental':not bool(model['available_public_guards_passed'])")
    target.write_text(text,encoding='utf-8',newline='\n')
print('Created separate V62 streaming builder, packager and integrity checker from reviewed V61 infrastructure.')
