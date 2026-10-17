"""Refresh final documentation and add the exact manual upload command."""
from pathlib import Path
import json
import build_submission as b
import document_v45_v50 as report
import update_v50_readme as readme
def main():
 A=b.ROOT/'experiments/v50_response_prior/analysis';ready=json.loads((A/'final_release_manifest.json').read_text());assert ready['passed'] and not ready['uploaded']
 report.main();readme.main();path=b.ROOT/'V45_V50_RESULTS.md';text=path.read_text(encoding='utf-8')
 text+='\n## Submit manually from Windows CMD\n\n'
 text+='The package is ready; the private score remains unknown. The aggregate public checks pass, with seed-sensitive gains.\n\n~~~bat\n'
 text+=r'cd /d L:\promec\Animesh\Download\Virtual'+'\n'
 text+=r'.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v50_response_prior\prediction.vcc -m "Empirical response-prior calibration v50" --wait'+'\n~~~\n'
 text+='\nThe final release manifest binds every proof receipt and the package SHA256. No upload or Git push was performed.\n';path.write_text(text,encoding='utf-8')
 print('Final V50 documentation updated',flush=True)
if __name__=='__main__':main()

