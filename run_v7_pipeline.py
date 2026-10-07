"""Finish V7 evaluation, fitting, count generation, packaging and independent audit."""
import json,subprocess,sys,time
from pathlib import Path
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v7_multisource_adaptive/analysis'
STAGES=[('original_source_count_replay','validate_v7_original_counts.py',[],ROOT/'analysis/v7_source_inventory/HEK293T_original_count_replay.json'),
 ('evaluation','evaluate_v7_response.py',[],A/'cell_selection.json'),
 ('production_fit','fit_v7_submission.py',[],A/'reference_model.json'),
 ('preflight','build_v7_submission.py',['--preflight','--workers','2'],A/'preflight.json'),
 ('model_audit','validate_v7_submission.py',['--model-only'],A/'model_validation.json'),
 ('generation_and_packaging','build_v7_submission.py',['--workers','3'],A/'submission_ready.json'),
 ('independent_package_audit','validate_v7_submission.py',[],A/'independent_validation.json'),
 ('documentation','document_v7_submission.py',[],ROOT/'V7_SUBMISSION.md')]
if __name__=='__main__':
    while not all((A/f).exists() for f in ['summary_selection.json','neural_manifest.json','ridge_manifest.json']):time.sleep(20)
    records=[]
    for name,file,args,sentinel in STAGES:
        print('V7 pipeline stage:',name,flush=True)
        if not sentinel.exists():
            result=subprocess.run([sys.executable,str(ROOT/file),*args],cwd=ROOT)
            if result.returncode:
                build.write_json(A/'pipeline_status.json',{'status':'failed','stage':name,'returncode':result.returncode,'completed':records});sys.exit(result.returncode)
        records.append({'stage':name,'script':file,'script_sha256':build.sha256(ROOT/file),'result_sha256':build.sha256(sentinel)})
        build.write_json(A/'pipeline_status.json',{'status':'running','completed':records})
    build.write_json(A/'pipeline_status.json',{'status':'verified_ready','completed':records,'uploaded':False})
    print('V7 pipeline verified ready; no upload performed.',flush=True)
