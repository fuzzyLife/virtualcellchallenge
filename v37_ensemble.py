"""Fixed V24-anchored ensembles; no treatment labels enter prediction or emission."""
import numpy as np
from scipy import sparse

ARMS = ('mean50', 'consensus', 'center_contrast50', 'population50', 'mean50_noise_mix')

def combine(base, retrieval, population, centered):
    b=np.asarray(base,np.float32);d=np.asarray(retrieval,np.float32)-b;p=np.asarray(population,np.float32)-b
    agree=(d*p>0)
    consensus=b+np.where(agree,np.sign(d)*np.minimum(abs(d),abs(p)),0)
    return {'v24':b.copy(),'mean50':(b+.5*d).astype(np.float32),
            'consensus':np.clip(consensus,-3,3).astype(np.float32),
            'center_contrast50':np.clip(b+.5*d+.5*(centered-b),-3,3).astype(np.float32),
            'population50':(b+.5*d).astype(np.float32),
            'mean50_noise_mix':(b+.5*d).astype(np.float32)}

def emit(template,effect,base,retrieval,own,seed,prepared,atlas,arm):
    from v36_methods import emit as parent
    if arm=='population50':
        a,diag=parent(template,base,own,seed,prepared,atlas)
        b,_=parent(template,retrieval,own,seed,prepared,atlas)
        # Balanced deterministic row assignment independent of expression/labels.
        n=template.shape[0];pick=np.random.default_rng(seed+370037).permutation(n)[:n//2]
        out=a.tolil();out[pick]=b[pick];out=out.tocsr();out.eliminate_zeros()
        libs,p,bulk=prepared;factor=np.exp2(np.clip(effect,-3,3));factor[own]=.1 if own is not None else factor[own]
        wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum()
        rp=(out.multiply(1/np.asarray(libs)[:,None])).mean(0).A1;rb=out.sum(0).A1/np.sum(libs)
        diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'mixture_rows':len(pick)}
        return out,diag
    if arm!='mean50_noise_mix':return parent(template,effect,own,seed,prepared,atlas)
    # Halfway between V11 compact (0 empirical weight) and V16/V24 residual10.
    libs,p,bulk=prepared;libs=np.asarray(libs,np.int64);factor=np.exp2(np.clip(effect,-3,3));factor[own]=.1 if own is not None else factor[own]
    wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum()
    expected=libs[:,None]*wp;empirical=template.toarray().astype(float)*factor
    empirical*=libs[:,None]/empirical.sum(1,keepdims=True)
    initial=.30*expected+.65*np.random.default_rng(seed).poisson(expected)+.05*empirical
    initial/=initial.sum(1,keepdims=True)
    out=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
    rp=(out/libs[:,None]).mean(0);rb=out.sum(0)/libs.sum()
    return sparse.csr_matrix(out),{'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'empirical_weight':.05}
