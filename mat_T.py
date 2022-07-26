# import mars.tensor as mt
from scipy.spatial.distance import cdist
import time
import numpy as np
from numba import njit
# import mars
# mars.new_session()

@njit(parallel=True)
def e(k,dist):
  return np.exp(-1j*k*dist)/dist

def mat_T(nodes, points, k):
    dst = cdist(nodes, points,"euclidean").T

    T = e(k,dst)
    return T

if __name__ == "__main__":
    for i in range(10)  :
        t = time.time()
        nodes = np.ones((2079,3))*1.3
        points = np.ones((11500,3))*1.1
        k = 2
        print(mat_T(nodes, points, k))
        print(time.time()-t)