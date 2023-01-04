import numpy as np
import matplotlib.pyplot as plt
plt.jet()
from utils import open, trat_sen
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
    FACTOR = 4.3498e+03
    t = time.localtime()

    # ref_ = open("resultados/" + reflectores_)
    # vrst_ = open("resultados/" + tejido_)
    # reff_ = np.abs(open("resultados/" + reflectores_)["sum_F_t"])
    # tej_ = np.abs(open("resultados/" + tejido_)["sum_F_t"])
    # sum_F_t_ = reff_/FACTOR + tej_

    # ref = open("resultados/" + reflectores)
    # vrst = open("resultados/" + tejido)
    reff = open("resultados/" + reflectores)["sum_F_t"]
    tej = open("resultados/" + tejido)["sum_F_t"]
    sum_F_t = reff/FACTOR + tej
    
    # name = str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    # im = plt.imshow(sum_F_t, extent=[0,ref["pos"][-1], ref["z0"]+(2048-2)/ref["fs"]/2*ref["ca"]-0.15, ref["z0"]-0.15], vmin= 0, vmax = 1)
    # plt.colorbar(im)
    # plt.savefig("Imagenes/python/pruebas/completo_" + name + ".png")

    return sum_F_t
    # print(np.max(sum_F_t), np.max(sum_F_t_))
    # err = np.sum((sum_F_t - sum_F_t_) ** 2)
    # err /= float(tej.shape[0] * tej.shape[1])
    # print(err)




    # name = "reflectores_" + str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    # plt.imshow(ref, extent=[0,vrsr["pos"][-1], vrsr["z0"]+(2048-2)/vrsr["fs"]/2*vrsr["ca"]-0.15, vrsr["z0"]-0.15])
    # plt.savefig("Imagenes/python/pruebas/" + name + ".png")

    # name = "tejido_" + str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    # plt.imshow(tej, extent=[0,vrst["pos"][-1], vrst["z0"]+(2048-2)/vrst["fs"]/2*vrst["ca"]-0.15, vrst["z0"]-0.15])
    # plt.savefig("Imagenes/python/pruebas/" + name + ".png")

    

def save_mat(file):
    name = "resultados/" + file
    vrs = open(name)
    scipy.io.savemat(name + ".mat", {"sum_F_t": vrs["sum_F_t"]})

def _print(lista):
    for i in lista:
        a  = open("resultados/" + str(i) + ".npz")
        print(str(i) + " : " + str(np.abs(np.max(a["sum_F_t"]))))

def filtro(file):
    t = time.localtime()
    img = open("resultados/" + file + ".npz")

    filtered = trat_sen(img["sum_F_t"], 0.1)
    # filtered = trat_sen(join( file + ".npz", "reflectores_11_25_20_31.npz" ), 0.1)

    name = str(t.tm_mday)  + "_" +  str(t.tm_mon) + "_" + str(t.tm_hour) + "_" + str(t.tm_min) + "_" + str(t.tm_sec)
    im = plt.imshow(filtered, extent=[0,img["pos"][-1], img["z0"]+(2048-2)/img["fs"]/2*img["ca"]-0.15, img["z0"]-0.15])
    plt.colorbar(im)
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.savefig("Imagenes/python/pruebas/" + name + ".png")

if __name__ == "__main__":
    # join( "tejido_11_24_15_27.npz","reflectores_11_25_20_31.npz","tejido_12_27_18_38.npz","reflectores_11_25_20_31.npz")
    # _print( ["tejido_12_14_20_49","reflectores_12_14_20_49"])
    # save_img("resultados/btejido_11_30_17_32.npz")
    # save_img("resultados/atejido_11_30_17_31.npz")
    # save_mat("tejido_12_30_20_5.npz")
    filtro("tejido_1_3_5_10")