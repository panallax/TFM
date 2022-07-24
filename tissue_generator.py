import numpy as np


def tissue_generator(x_t,y_t,z1,z2, ro):
    tissue_coords = np.array([[x_t[0],x_t[-1]], [y_t[0],y_t[-1]], [z1,z2]]).T

    n = round(ro*(x_t[-1] - x_t[0])*(y_t[-1] - y_t[0])*(z2 - z1))
    tissue = np.zeros((n,3))

    for i in range(tissue_coords.shape[0] + 1):
        tissue[:,i] = tissue_coords[0,i] + (tissue_coords[1,i] - tissue_coords[0,i])*np.random.uniform(0, 1, n)

    return tissue, n

if __name__ == "__main__":
    def linspace(start, stop, step=1.):
        return np.linspace(start, stop, int(np.ceil((stop - start) / step + 1)))

    print(tissue_generator(linspace(-0.8,3.8, 0.04 ), linspace(-2, 2, 0.04), 9.8, 11.2, 5))