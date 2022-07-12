from scipy.spatial.distance import cdist
import numpy as np

def mat_T(nodes, points, k):
    dst = cdist(nodes, points,"euclidean").T
    T = np.exp(-1j*k*dst)/dst
    return T

if __name__ == "__main__":
    nodes = np.ones((3,3))*1.3
    points = np.ones((3,3))*1.1
    k = 2

    print(mat_T(nodes, points, k))