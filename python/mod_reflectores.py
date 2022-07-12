import numpy as np
import scipy.special as sp
import mpmath
mpmath.mp.pretty = True

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

if __name__ == "__main__":
    a = mod_reflectores(1073e-9, 1.638, 1090e-9, 1.672, 0, 5e-3, np.linspace(0.01, 150, 2048))
    print(a)