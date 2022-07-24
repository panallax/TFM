import numpy as np

def linspace(start, stop, step=1.):
  return np.linspace(start, stop, int(np.ceil((stop - start) / step + 1)))

def interfase_generator(x,y,z):
    xv, yv = np.meshgrid(x, y)
    xy_pairs = np.vstack([xv.reshape(-1), yv.reshape(-1)])
    return np.concatenate([xy_pairs.T, np.ones((xy_pairs.T.shape[0], 1))*z], axis = 1)


if __name__ == "__main__":
    x = linspace(-0.8, 3.8, 0.04)
    y = linspace(-2, 2, 0.04)

    a = interfase_generator(x, y, 11.2)
    print(a)