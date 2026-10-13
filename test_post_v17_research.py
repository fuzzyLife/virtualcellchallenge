"""Meaningful checks for the non-executing vocabulary reader and feature coverage."""
import pickle
import numpy as np,json
from audit_post_v17_feasibility import decode_unicode_list,A

def main():
 for values in [['A1BG','THRA1/BTR'],[np.str_('A1BG'),np.str_('THRA1/BTR')]]:
  assert decode_unicode_list(pickle.dumps(values,protocol=4))==['A1BG','THRA1/BTR']
 try:decode_unicode_list(b'cos\nsystem\n(S"echo should-not-run"\ntR.')
 except ValueError:pass
 else:raise AssertionError('Executable pickle globals must not be interpreted')
 genes=decode_unicode_list((A/'stack_gene_vocabulary.pkl').read_bytes());assert len(genes)==len(set(genes))==15012 and 'THRA1/BTR' in genes
 report=json.loads((A/'feasibility_audit.json').read_text());summary={r['asset']:r for r in report['coverage']};assert summary['Stack published vocabulary']['output_genes']==13078;assert summary['PIE gene_text']['targets']==300 and summary['PIE esm2']['targets']==299
 print('Vocabulary/schema/coverage audit checks passed; no model prediction or biological validation claimed.',flush=True)
if __name__=='__main__':main()
