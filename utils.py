import numpy as np
import scipy.special as sp
import mpmath
mpmath.mp.pretty = True
from scipy.spatial.distance import cdist
from numba import njit

def tissue_generator(x_t,y_t,z1,z2, ro):
    tissue_coords = np.array([[x_t[0],x_t[-1]], [y_t[0],y_t[-1]], [z1,z2]]).T

    n = round(ro*(x_t[-1] - x_t[0])*(y_t[-1] - y_t[0])*(z2 - z1))
    tissue = np.zeros((n,3))

    for i in range(tissue_coords.shape[0] + 1):
        tissue[:,i] = tissue_coords[0,i] + (tissue_coords[1,i] - tissue_coords[0,i])*np.random.uniform(0, 1, n)

    return tissue, n


def mod_reflectores(rof,c,rop,cpl,cps,a,frecs):

    w=2*np.pi*frecs
    sigma=0.5*(1-2*(cps/cpl)**2)/(1-(cps/cpl)**2)

    if sigma==0.5:
        sigma=0.499999999
        cps=0.000000001
    
    k=w/c
    k1=w/cpl
    k2=w/cps

    x=k*a
    x1=k1*a

    x2=k2*a

    nfo = 5        

    sl2 = np.zeros_like(frecs)
    ln2 = np.zeros((nfo+1, frecs.shape[0]), np.complex)

    for n in range(nfo + 1):

        jn_1=np.sqrt(np.pi/2/x1)*sp.jn(n+0.5,x1)
        jnmas1_1=np.sqrt(np.pi/2/x1)*sp.jn(n+1.5,x1)
        jnp_1=n/x1*jn_1-jnmas1_1
        jns_1=(n**2-n-x1**2)/x1**2*jn_1+2*jnmas1_1/x1

        jn_2=np.sqrt(np.pi/2/x2)*sp.jn(n+0.5,x2)
        jnmas1_2=np.sqrt(np.pi/2/x2)*sp.jn(n+1.5,x2)
        jnp_2=n/x2*jn_2-jnmas1_2
        jns_2=(n**2-n-x2**2)/x2**2*jn_2+2*jnmas1_2/x2

        jn=np.sqrt(np.pi/2/x)*sp.jn(n+0.5,x)
        jnmas1=np.sqrt(np.pi/2/x)*sp.jn(n+1.5,x)

        jnp=n/x*jn-jnmas1

        nn=np.sqrt(np.pi/2/x)*np.vectorize(lambda x: float(mpmath.bessely(n+0.5,x)))(x)
        nnmas1=np.sqrt(np.pi/2/x)*np.vectorize(lambda x: float(mpmath.bessely(n+1.5,x)))(x)
        
        nnp=n/x*nn-nnmas1

        hn=jn+1j*nn
        hnp=jnp+1j*nnp
        
        Mn=2*(jn_1-x1*jnp_1)/((n**2+n-2)*jn_2+jns_2*x2**2)

        A11=-rof/x*(x1*jnp_1+Mn*n*(n+1)*jn_2)

        A12=-hnp
        A21=-rop/(1-sigma)*(sigma*jn_1+(1-2*sigma)*(-jns_1+n*(n+1)*Mn*(jn_2-x2*jnp_2)/x1**2))
        A22=-hn
        
        ln2[n,:]=(-1)**n*(2*n+1)*(A11*jn-A21*jnp)/(A11*A22-A21*A12)/k

        sl2=sl2+ln2[n,:]
    
    return sl2


# @njit(parallel=True)
def e(k,dist):
  return np.exp(-1j*k*dist)/dist

def mat_T(nodes, points, k):
    dst = cdist(nodes, points,"euclidean").T

    T = e(k,dst)
    return T


def linspace(start, stop, step=1.):
  return np.linspace(start, stop, int(np.ceil((stop - start) / step + 1)))

def interfase_generator(x,y,z):
    xv, yv = np.meshgrid(x, y)
    xy_pairs = np.vstack([xv.reshape(-1), yv.reshape(-1)])
    return np.concatenate([xy_pairs.T, np.ones((xy_pairs.T.shape[0], 1))*z], axis = 1)


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

def generate_mesh(r,h,d):
    dth = d/r
    theta = 1e-9
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

def open(file):
    
    vrs = np.load(file)
    msg = "Vars loaded: {}".format(" ,".join([str(x) for x in list(vrs.keys())]))
    print(msg)
    return vrs

if __name__ == "__main__":
    def linspace(start, stop, step=1.):
        return np.linspace(start, stop, int(np.ceil((stop - start) / step + 1)))

    print(tissue_generator(linspace(-0.8,3.8, 0.04 ), linspace(-2, 2, 0.04), 9.8, 11.2, 5))

    a = mod_reflectores(1073e-9, 1.638, 1090e-9, 1.672, 0, 5e-3, np.linspace(0.01, 150, 2048))
    print(a)

    nodes = np.ones((2079,3))*1.3
    points = np.ones((11500,3))*1.1
    k = 2
    print(mat_T(nodes, points, k))

    x = linspace(-0.8, 3.8, 0.04)
    y = linspace(-2, 2, 0.04)

    a = interfase_generator(x, y, 11.2)
    print(a)

    a = generate_random_points(5, 2, 6, 10, 1)
    print(a)

    rc = 14
    r = 3
    h = rc*(1-np.sqrt(1-(r/rc)**2))
    a = generate_mesh(rc,h, 1.5/40)