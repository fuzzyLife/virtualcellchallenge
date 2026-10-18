"""Pool sufficient statistics before shrinkage; preserve matched lane controls."""
import numpy as np
from v61_pooled_state_response import aggregate

def main():
 # Many tiny lanes with the same true conditional change must not be shrunk
 # as if the entire experiment contained only one treated cell per state.
 nc=np.array([20.,20.,20.,20.]);nt=np.array([1.,1.,1.,1.]);cm=np.array([[100.,200.]]*4);tm=cm*2
 lane=(nt,nc,tm,cm,tm.mean(0),cm.mean(0));one=aggregate([lane]);pooled=aggregate([lane]*20)
 assert (pooled['within']>one['within']).all();assert (pooled['within']>.5).all()
 np.testing.assert_allclose(pooled['tilt'],np.zeros(4));np.testing.assert_allclose(pooled['cm'],cm)
 # Zero response across differently sized lanes remains zero.
 lane2=(nt*3,nc,cm,cm,cm.mean(0),cm.mean(0));zero=aggregate([(nt,nc,cm,cm,cm.mean(0),cm.mean(0)),lane2])
 np.testing.assert_allclose(zero['within'],0);np.testing.assert_allclose(zero['tilt'],0)
 print('V61 pooling/zero-response contracts passed',flush=True)
if __name__=='__main__':main()
