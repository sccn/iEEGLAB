import numpy as np, os, sys, time
from picard import picard
H=os.path.expanduser('~')
mode=sys.argv[1]; seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
M=np.load(H+'/ica_matrix.npy').astype(np.float64)[:, ::2]      # 60k samples
if mode=='car':
    M = M - M.mean(0, keepdims=True)
k = M.shape[0] - (1 if mode=='car' else 0)
mu=M.mean(1,keepdims=True); Ac=M-mu
U,s,_=np.linalg.svd(Ac@Ac.T/Ac.shape[1])
W=(U[:,:k]/np.maximum(s[:k],1e-30)**.5).T
t0=time.time()
K_,Wp,_=picard(W@Ac, ortho=False, extended=True, max_iter=100, random_state=seed, tol=1e-5)
unmix=Wp@K_@W
np.savez(H+f'/ica_{mode}{seed}.npz', S=(unmix@Ac).astype(np.float32), A=np.linalg.pinv(unmix), eig=s)
print(mode, 'done %.0fs' % (time.time()-t0), 'eig min/max %.2e' % (s[-1]/s[0]), flush=True)
