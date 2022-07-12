import numpy as np


def generate_mesh(r,h,d):
    dth = d/r
    theta = 0
    nodes = np.zeros((1,3))

    z = 0

    while z <= h:
        z = r*(1 - np.cos(theta))
        Ra = r*np.sin(theta)
        phi = 0

        n = int(np.round(2*np.pi/(10*(d/Ra))))
        m = np.zeros((n,3))

        for i in range(n):

            m[i,0] = Ra*np.cos(phi)
            m[i,1] = Ra*np.sin(phi)
            m[i,2] = z
            phi = phi + 10*d/Ra

        nodes = np.concatenate([nodes, m])
        theta = theta + dth
    
    return nodes

if __name__ == "__main__":
    rc = 14
    r = 3
    h = rc*(1-np.sqrt(1-(r/rc)**2))
    a = generate_mesh(rc,h, 1.5/40)