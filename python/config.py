import numpy as np

####### DATOS NUBE PUNTOS #####
x = 5
y = 2
z = 6
zmin = 10
ro = 10
n_r = round(x*y*z*ro)
n_r = 1
###### ATENUACIÓN #########
at = 2.15e-5
a_r = 5e-3 #
a_t = 15e-3
######## DATOS CASQUETE #######
d = 1.5/20/2                     #Distancia mínima entre nodos (1.5 es ca)
r = 3       #                   #Radio casquete
Rc = 14     #                   #Radio de curvatura
h = Rc*(1-np.sqrt(1-(r/Rc)^2))    #Altura casquete

######### TEJIDOS E INTERFASES ######
x_t = np.arange(-0.8,3.8,0.04)
y_t = np.arange(-2,2,0.04)
z_1 = 9.8
z_2 = 11.2
ro_t = 5


#### FLUIDO 1 y 2
rot = 1073e-9 #*e-9 # densidad tejido
ct = 1.638   # velocidad tejido
zt = rot*ct # impedancia tejido
ca = 1.5      # velocidad agua
roa = 998e-9  #densidada agua
za = roa*ca #impedancia agua
c_r = (za-zt)/(za+zt) #coeficiente reflexión
c_t = 2*za/(zt+za) #coeficiente transmisión


####Propiedades núcleos ####
rhon = 1090e-9
cpln = 1.672
a_n = 5e-3
cps = 0

rop = 1064e-9
cpl = 1.585