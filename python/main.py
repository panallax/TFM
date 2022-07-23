import numpy as np
from scipy.spatial.distance import squareform, pdist, cdist
import scipy.oi
import scipy.fft
from scipy.signal.windows import tukey
import matplotlib.pyplot as plt
plt.jet()

from generate_mesh import generate_mesh
from generate_random_points import generate_random_points
from interfase_generator import interfase_generator 
from mat_T import mat_T
from mod_reflectores import mod_reflectores
from tissue_generator import tissue_generator
from config import *



Nodes = generate_mesh(Rc,h,d)
# reflector_points = generate_random_points(x, y, z, zmin, n_r) 
points = np.array([1.5, 0, 13])
(tissue_points,n_t) = tissue_generator(x_t,y_t,z_1,z_2, ro_t)
# interfase_points_1 = interfase_generator(x_t,y_t,z_1)
interfase_points = interfase_generator(x_t,y_t,z_2)
# points = reflector_points
# points = np.concatenate(tissue_points, interfase_points)
sub = squareform(pdist(points))

aj = np.concatenate(tissue_points, interfase_points)
# points = load('/Users/alex/Desktop/TFM/TFM/validate clutter/points.mat').points
# sub = squareform(pdist(tissue_points))

n = len(points)
# n = n_r + n_t

############################################### RESTRINGIR FRECUENCIAS 1
#### Espacio de análisis, 7mm, Tiempo de análisis=7*2/1.5=9.3us
#### Podemos retrasar 2*9/1.5= 12us
#### [~,indt]=min(abs(sen.t-9.3))=1396, reducimos a 2048 puntos
###################### (nos quedamos con el intervaslo de tiempo
###################### estrictamente necesario: 2048 puntos)
sen = scipy.io.loadmat('sen_Alex.mat')
frecs = np.linespace(0.01,150,2048)
senc = sen.sen[1:2048]
pos = np.arange(0,3,0.016)

####### ATENUACIÓN EN K ############
cct = ct*(1+1j*ct/2/np.pi*at*frecs)
cca = ca*(1+1j*ca/2/np.pi*at*frecs)
kt = 2*np.pi*frecs/cct
ka = 2*np.pi*frecs/cca


SL2 = mod_reflectores(rot,ct,rhon,cpln,cps,a_n,frecs) # tejido
SL2_r = mod_reflectores(roa,ca,rop,cpl,cps,a_r,frecs) #reflectores

E = scipy.fft(senc)
E[round(len(E)/2)+1:] = 0

################### RESTRINGIR FRECUENCIAS 2 ##############################

imax = max(abs(E[1:len(E)/2]+1))
v = np.where(abs(E[1:len(E)/2]) == imax)
freqind0 = np.where(abs(E[1:len(E)/2]) > v/5)
Ev = np.zeros((1,2048))
Ev[(freqind0):max(freqind0)] = E[min(freqind0):max(freqind0)]*tukey(max(freqind0)-min(freqind0)+1,0.4).T
freqind = np.where(abs(Ev)>0)

sum_F = np.zeros(len(frecs),len(pos))
sum_FR = np.zeros(len(frecs),len(pos))


f = len(frecs)
t2 = cdist(interfase_points, points,"euclidean").T

nodes_steps = np.tile(Nodes[:,:,None], (1, 1, len(pos)))
nodes_steps = np.transpose(nodes_steps,(0,2,1))
vec = np.arange(1,len(pos)+1,1)
M = np.tile(vec, (len(Nodes), 1))
nodes_steps[:, :, 0] = Nodes[:, 0,None]  + M*(pos[2]-pos[1])
node_steps = np.transpose(nodes_steps, (0,2,1))

for j in range(len(pos)):
    print(j)
    v = np.zeros((1,f))
    for i in range(f):
        if np.nonzero(freqind==i):
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
            F = T1.T*T2.T*SL2_r[i]*(T2*sum(T1,2))
            v[i] = sum(F)*E[i]

    sum_FR[:,j]=v.T*np.exp(1j*2*np.pi*(frecs.T)/ct*zmin)

sum_FR = sum_FR*c_t**2
sum_F_t = scipy.fft.ifft(sum_FR)


plt.imshow(abs(sum_F_t), extent=[0,pos,z+zmin,zmin])
plt.show()


# scatter3(points(:,1), points(:,2), points(:,3), "filled")
# hold on
# scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")
