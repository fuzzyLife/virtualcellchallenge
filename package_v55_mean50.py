"""V55 adapter around the frozen full-payload-verified V50 official packager."""
from pathlib import Path
import package_v50_submission_v2 as p
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v55_mean50_diagnostic';A=OUT/'analysis'
def main():
 p.OUT=OUT;p.A=A;p.main()
 # Preserve inherited provenance fields, explicitly identify the actual adapter.
 import json
 ready=json.loads((A/'submission_ready.json').read_text())
 ready.update(experiment='V55 fixed mean50 diagnostic',actual_generator_adapter_sha256=p.build.sha256(ROOT/'run_v55_mean50.py'),actual_packaging_adapter_sha256=p.build.sha256(Path(__file__)),actual_emitter='v36_methods.emit with dispersion=None -> inherited V24 residual10',inherited_packager_provenance_fields_describe_V50_template=True,diagnostic_submission=True,original_V37_strict_guards_passed=False)
 p.build.write_json(A/'submission_ready.json',ready)
if __name__=='__main__':main()
