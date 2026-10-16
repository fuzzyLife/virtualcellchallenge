"""Lock one candidate across primary and corrective ensemble tests before count confirmation."""
import json
import v35_methods as m
from run_v37_ensemble import A,read,write

def main():
 paths={'primary':A/'selection.json','recalibrated':A/'recalibrated/selection.json'};eligible=[]
 noise=A/'noise_recalibrated/selection.json';paths['noise_recalibrated']=noise if noise.exists() else A/'noise_recalibrated/null_gate.json'
 for namespace,path in paths.items():
  s=read(path)
  for arm in s['eligible']:eligible.append((s['gain_scores'][arm],namespace,arm))
 winner=max(eligible) if eligible else (0.,'primary','v24');file=A/'ensemble_release.json';assert not file.exists()
 write(file,{'selected_arm':winner[2],'namespace':winner[1],'selection_gain':winner[0],'eligible':eligible,'selection_inputs_sha256':{n:m.sha(p) for n,p in paths.items()},'rule':'Same frozen acceptance rules for all candidates; among qualifying arms select greatest oriented relative mean gain. Lock before additional count confirmation. Corrective arm explicitly added after initial diagnostic outcomes, not blind discovery.','confirmation_required':winner[2]!='v24','D3_opened':False,'confirmed_leaderboard_gain':False,'selector_sha256':m.sha(m.ROOT/'select_v37_final.py')});print(json.dumps(read(file),indent=2),flush=True)

if __name__=='__main__':main()
