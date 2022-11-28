import numpy as np
import matplotlib.pyplot as plt
plt.jet()
from utils import open
import sys
import time
from os import path
import scipy.io

def save_img(npz):
    t = time.localtime()
    vrs = open(npz)
    sum_F_t = vrs["sum_F_t"]
    name = str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    plt.imshow(np.abs(sum_F_t), extent=[0,vrs["pos"][-1], vrs["z0"]+(2048-2)/vrs["fs"]/2*vrs["ca"]-0.15, vrs["z0"]-0.15])
    plt.savefig("Imagenes/python/pruebas/" + name + ".png")

    # ax = plt.axes(projection='3d')
    # ax.scatter3D(vrs["points"][:,0],vrs["points"][:,1],vrs["points"][:,2])
    # plt.savefig("Imagenes/python/pruebas/" + "name1" + ".png")

def join(tejido, reflectores):
    t = time.localtime()
    ref = open("resultados/" + reflectores)
    tej = open("resultados/" + tejido)
    print(np.mean(ref["sum_F_t"]), np.mean(tej["sum_F_t"]))
    sum_F_t = ref["sum_F_t"] + tej["sum_F_t"]
    name = str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    plt.imshow(np.abs(sum_F_t), extent=[0,ref["pos"][-1], ref["z0"]+(2048-2)/ref["fs"]/2*ref["ca"]-0.15, ref["z0"]-0.15])
    # plt.savefig("Imagenes/python/pruebas/completo_" + name + ".png")

def save_txt(file):
    vrs = open("resultados/" + file)
    np.savetxt(file.split(".")[0] + ".txt",vrs["a"])

if __name__ == "__main__":
    join("tejido_11_24_15_27.npz", "reflectores_11_25_20_31.npz")
    # save_img("resultados/tejido_11_24_15_27.npz")
    # save_txt