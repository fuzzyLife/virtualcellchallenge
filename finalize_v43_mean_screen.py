"""Record the V43 mean-control failure without interpreting the shuffled arm as a release."""
import v43_lincs_fingerprints as m
def main():
 m.check();r=m.read(m.A/'mean_results.json');control=r['summary']['shuffled250'];checks={}
 for a in m.ARMS[:2]:checks[a]={k:r['summary'][a][k]<control[k] for k in ['mse_ratio','mae_ratio']}
 assert not any(all(x.values()) for x in checks.values())
 m.write(m.A/'mean_stop_decision.json',{'selected':'v24','release':False,'matched_control_checks':checks,'mean_results_sha256':m.b.sha256(m.A/'mean_results.json'),'count_tests_opened':False,'additional24_count_cases_opened':False,'scope':'Generic mean-stage selected shuffled250 as smallest error, but that is a negative control, never a release candidate. Genuine transfer fails the frozen matched mean-control condition; count/release stages are not justified. Keep all models for future research.','script_sha256':m.b.sha256(m.ROOT/'finalize_v43_mean_screen.py')})
 print('V43 genuine candidates stopped by matched mean controls',flush=True)
if __name__=='__main__':main()

