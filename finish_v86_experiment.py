"""Run the fixed V86 experiment with bounded workers and no automatic upload."""
import subprocess,sys
import v86_joint_posterior as m
from v57_pie import read,write,sha

def batch(jobs):
    for start in range(0,len(jobs),2):
        workers=[]
        for label,args in jobs[start:start+2]:
            log=(m.OUT/(label+'.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',*args],stdout=log,stderr=subprocess.STDOUT);workers.append((label,p,log));print('V86 started',label,flush=True)
        for label,p,log in workers:
            code=p.wait();log.close();assert code==0,(label,code);print('V86 completed',label,flush=True)

def main():
    jobs=[]
    for outer in ['H1','Orion','CD4']:
        for shuffle in [False,True]:jobs.append(('fit_'+outer+('_shuffled' if shuffle else ''),['run_v86_experiment.py','fit','--outer',outer]+(['--shuffled'] if shuffle else [])))
    batch(jobs)
    batch([('prepare',['run_v86_experiment.py','prepare'])])
    subprocess.run([sys.executable,'-X','utf8','run_v86_experiment.py','means'],check=True)
    subprocess.run([sys.executable,'-X','utf8','audit_v86_experiment.py','models'],check=True)
    for offset in [0,620000]:
        if offset and not (read(m.A/'mean_results.json')['passes'] and read(m.A/'count_assessment_0.json')['passes']):break
        batch([(f'counts_{offset}_{chunk}',['run_v86_counts.py','counts','--offset',str(offset),'--chunk',str(chunk),'--chunks','2']) for chunk in range(2)])
        subprocess.run([sys.executable,'-X','utf8','run_v86_counts.py','assess','--offset',str(offset)],check=True)
        subprocess.run([sys.executable,'-X','utf8','audit_v86_experiment.py','counts','--offset',str(offset)],check=True)
    confirmed=(m.A/'count_assessment_620000.json').exists()
    qualified=read(m.A/'mean_results.json')['passes'] and read(m.A/'count_assessment_0.json')['passes'] and confirmed and read(m.A/'count_assessment_620000.json')['passes']
    write(m.A/'experiment_complete.json',dict(completed=True,second_seed_run=confirmed,strict_support=bool(qualified),new_upload_generated=False,runner_sha256=sha(__file__)))
    print('V86 completed, strict support',qualified,flush=True)
if __name__=='__main__':main()
