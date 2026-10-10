"""Run the independently written complete source/cell identity audit on V13 plans."""
from v13_methods import A
import audit_v12_splits as independent

def main():
 independent.A=A
 for folder in ['', 'robust_stage','quality_stage']:independent.main(folder)
if __name__=='__main__':main()
