"""Summarize requested/realized generation moments without claiming biological accuracy."""
from pathlib import Path
import json,hashlib
import pandas as pd
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def summarize(stem,destination,expected):
 path=A/(stem+'.csv')
 if not path.exists():return
 stats=pd.read_csv(path);assert len(stats)==expected;summary={}
 for c,f in stats.groupby('context'):
  summary[str(c)]={k:dict(mean=float(f[k].mean()),max=float(f[k].max())) for k in ['mean_probability_L1','bulk_probability_L1','changed_transcript_fraction']};summary[str(c)]['groups']=len(f);summary[str(c)]['stored_nonzeros']=int(f.nnz.sum())
 result=dict(groups=len(stats),contexts=summary,stored_nonzeros=int(stats.nnz.sum()),source_statistics_sha256=sha(path),summarizer_sha256=sha(Path(__file__)),scope='Requested versus realized count-construction moments only; no held-out final biological accuracy is available.')
 (A/destination).write_text(json.dumps(result,indent=2),encoding='utf-8');print(destination,len(stats),'groups summarized',flush=True)
def main():
 summarize('preflight_group_statistics','preflight_moment_summary.json',12);summarize('prediction_group_statistics','all_group_moment_summary.json',900)
if __name__=='__main__':main()
