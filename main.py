import numpy as np
from scipy.spatial.distance import squareform, pdist, cdist
import scipy.io
from scipy.signal.windows import tukey
import matplotlib.pyplot as plt
plt.jet()
import time
from multiprocessing import Pool

from utils import *
from config import *



def runmain(j,f,t2,freqind,nodes_steps,interfase_points,kt,ka,SL2_r,E):
    print(j)
    t2_ = time.time()
    v = np.zeros((f),dtype = 'complex_')
    for i in range(f):
        if np.nonzero(freqind==i)[0].size > 0:
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
            v[i] = np.sum(F)*E[i]
            # print(time.time()-t2_)

    return (j,v*np.exp(1j*4*np.pi*(frecs.T)/ct*z0))
    


if __name__ == "__main__":
    t = time.time()

    Nodes = generate_mesh(Rc,h,d)
    # reflector_points = generate_random_points(x, y, z, zmin, n_r) 
    # points = np.array([[1.5, 0, 13]])

    (tissue_points,n_t) = tissue_generator(x_t,y_t,z_1,z_2, ro_t)
    points = tissue_points
    interfase_points = interfase_generator(x_t,y_t,z_2)
    # points = np.concatenate(tissue_points, interfase_points)
    # sub = squareform(pdist(points))

    # aj = np.concatenate((tissue_points, interfase_points), axis=0)
    # points = load('/Users/alex/Desktop/TFM/TFM/validate clutter/points.mat').points

    # n = len(points)
    n = n_r + n_t

    ############################################### RESTRINGIR FRECUENCIAS 1
    #### Espacio de analisis, 7mm, Tiempo de analisis=7*2/1.5=9.3us
    #### Podemos retrasar 2*9/1.5= 12us
    #### [~,indt]=min(abs(sen.t-9.3))=1396, reducimos a 2048 puntos
    ###################### (nos quedamos con el intervaslo de tiempo
    ###################### estrictamente necesario: 2048 puntos)
    sen = scipy.io.loadmat('sen_Alex.mat')
    frecs = np.linspace(0.01,150,2048)
    f = len(frecs)
    senc = sen["sen"][0][:2048]
    pos = np.arange(0,3,0.016)

    Ev = np.zeros((1,2048), dtype = 'complex_')
    sum_F = np.zeros((len(frecs),len(pos)), dtype = 'complex_')
    sum_FR = np.zeros((len(frecs),len(pos)), dtype = 'complex_')

    ####### ATENUACION EN K ############
    cct = ct*(1+1j*ct/2/np.pi*at*frecs)
    cca = ca*(1+1j*ca/2/np.pi*at*frecs)
    kt = 2*np.pi*frecs/cct
    ka = 2*np.pi*frecs/cca

    SL2 = mod_reflectores(rot,ct,rhon,cpln,cps,a_n,frecs) # tejido
    SL2_r = mod_reflectores(roa,ca,rop,cpl,cps,a_r,frecs) #reflectores

    E = scipy.fft.fft(senc)
    E[round(len(E)/2)-1:] = 0
    ################### RESTRINGIR FRECUENCIAS 2 ##############################

    imax = np.argmax(abs(E[:int(len(E)/2)-1]))
    v = max(abs(E[:int(len(E)/2)-1]))
    freqind0 = np.nonzero(abs(E[:int(len(E)/2)]) > v/5)[0]
    Ev.flat[freqind0[0]:freqind0[-1]+1] = E[freqind0[0]:freqind0[-1]+1]*tukey(max(freqind0)-min(freqind0)+1,0.4).T
    freqind = np.nonzero(abs(Ev) > 0)[1]

    t2 = cdist(interfase_points, points,"euclidean").T

    nodes_steps = np.tile(Nodes[:,:,None], (1, 1, len(pos)))
    nodes_steps = np.transpose(nodes_steps,(0,2,1))
    vec = np.arange(1,len(pos)+1,1)
    M = np.tile(vec, (len(Nodes), 1))
    nodes_steps[:, :, 0] = Nodes[:, 0,None]  + M*(pos[2]-pos[1])
    nodes_steps = np.transpose(nodes_steps, (0,2,1))

    pool = Pool()
    result = pool.starmap_async(runmain, ((j,f,t2,freqind,nodes_steps,interfase_points,kt,ka,SL2_r,E) for j in range(len(pos))))

    for (m,n) in result.get():
        sum_FR[:,m] = n

    sum_FR = sum_FR*c_t**2
    sum_F_t = np.fft.ifft(sum_FR, axis=0)
    a = np.abs(sum_F_t)

    print(time.time()-t)
    plt.imshow(a, extent=[0,pos[-1], z0+(2048-1)/fs/2*ca, z0])
    plt.savefig("Imagenes/python/"+str(time.time())+".png" )

    np.savez("resultados/"+str(time.time()), sum_FR=sum_FR, sum_F_t=sum_F_t, pos=pos, z0=z0, fs=fs, ca=ca, a=a, points=points, interfase_points = interfase_points )

# scatter3(points(:,1), points(:,2), points(:,3), "filled")
# hold on
# scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")
