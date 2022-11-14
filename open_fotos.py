import numpy as np
import matplotlib.pyplot as plt
from utils import open
import sys
import time

vrs = open(str(sys.argv[-1]))
sum_F = abs(vrs["sum_F_t"])
print(sum_F)
plt.imshow(vrs["a"], extent=[0,vrs["pos"][-1], vrs["z0"]+(2048-1)/vrs["fs"]/2*vrs["ca"], vrs["z0"]])
plt.savefig("Imagenes/python/pruebas/"+str(time.time())+".png")