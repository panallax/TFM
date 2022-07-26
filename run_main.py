from numba import jit, prange
import numpy as np
import time
import numpy as np
import scipy.io
import scipy.fft

# import mars.tensor as np
import matplotlib.pyplot as plt
plt.jet()

from mat_T import mat_T

@jit(parallel=True)
def run_main(pos,f,freqind,nodes_steps,interfase_points,kt,ka,SL2_r,c_t,zmin,t2,E,frecs,ct):

    sum_FR = np.zeros((len(frecs),len(pos)))           
    for j in range(len(pos)):
        print(j)
        v = np.zeros((1,f))
        for i in prange(f):
            if np.nonzero(freqind==i):
                t = time.time()
    #             E_ = E[i]*np.ones(len(Nodes),1)
    #             T = mat_T(Nodes, points, k[i])
                T1 = mat_T(nodes_steps[:,:,j], interfase_points, kt[i])
    #             T2 = mat_T(interfase_points, points, ka[i])
                T2 = np.exp(-1j*ka[i]*t2)/t2
    #             intern_dist = np.exp(-1j*k[i]*sub)/sub
    #             intern_dist(1:1+size(intern_dist,1):end) = 0
    #             R = SL2[i]**2*intern_dist + SL2[i]*np.eye(n)
    #             R = (SL2[i]*np.ones(1,nt) + c_r*np.ones(1,n-nt))*np.ones(n)
    #             R = SL2_r[i]*np.ones(n)
    #             R = sparse(R)
    #             F = T.T*(R*sum(T,2))
                F = np.matmul(T1.T,T2.T)*SL2_r[i]*np.matmul(T2,np.sum(T1,1))
                v.flat[i] = np.sum(F)*E[i]
                print(time.time()-t)
                
        sum_FR.flat[:,j]=v.T*np.exp(1j*2*np.pi*(frecs.T)/ct*zmin)

    sum_FR = sum_FR*c_t**2
    sum_F_t = scipy.fft.ifft(sum_FR)
    return sum_F_t