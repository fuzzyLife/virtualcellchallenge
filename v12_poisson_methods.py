"""Frozen V11 response with a bounded increase of Poisson variation."""
from v11_methods import emit as original_emit,SETTINGS
NOISE={'v11':.65,'noise85':.85,'noise100':1.}
def emit(x,rows,effect,own,seed,variant,prepared=None,atlas=None):
 setting={**SETTINGS['consensus20_noise65'],'noise':NOISE[variant]}
 return original_emit(x,rows,effect,own,seed,setting,prepared=prepared,atlas=atlas)
