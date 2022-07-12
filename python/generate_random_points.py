import numpy as np

def generate_random_points(x,y,z,zmin,n):

    dims = [x, y, z]

    coords = np.zeros((2, 3))
    points = np.zeros((n, 3))

    coords[:,0] = np.array([-0.8, x + 0.8])
    coords[:,1] = np.array([-y/2, y/2])
    coords[:,2] = np.array([zmin, zmin+z])

    for i in range(3):

        points[:,i] = coords[0,i] + (coords[1,i] - coords[0,i])*np.random.uniform(0, 1, n)
    
    return points


if __name__ == "__main__":
    a = generate_random_points(5, 2, 6, 10, 1)
    print(a)