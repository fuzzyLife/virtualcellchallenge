"""Run frozen production inference with the tested pandas metadata compatibility flag."""
from pathlib import Path
import json,hashlib,faulthandler
import pandas as pd
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis'
def main():
 faulthandler.enable();pd.options.future.infer_string=False
 from infer_v20_submission import main as infer
 record={'pandas_future_infer_string':False,'scope':'Metadata-only CSV string compatibility, same flag as original-author smoke/parity test. Production reader explicitly converts axes to Unicode. Numeric inputs, model weights, correction, seeds and count emitter unchanged.','wrapper_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frozen_inference_sha256':hashlib.sha256((ROOT/'infer_v20_submission.py').read_bytes()).hexdigest(),'observed_native_failure_before_flag':'Windows access violation 0xc0000005 in pandas string_arrow._from_sequence while build_submission.read_spec reads CSV, before first target inferred.'}
 (A/'production_runtime_compatibility.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');infer()
if __name__=='__main__':main()
