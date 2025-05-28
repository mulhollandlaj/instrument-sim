import struct
import numpy as np
import matplotlib.pyplot as plt


with open("output/sdata.bin", "rb") as file:
    nx = int(struct.unpack("i", file.read(4))[0])
    length = float(struct.unpack("f", file.read(4))[0]) * 1e3

    h = length / float(nx)
    
    sdata = np.zeros((nx), np.float32)
    xdata = np.zeros((nx), np.float32)

    for i in range(nx):
        sdata[i] = (np.float32) (struct.unpack("f", file.read(4)))[0]
        xdata[i] = h * (i + 0.5)


plt.yscale('log')
plt.plot(xdata, sdata,'k')
plt.grid(which='both')
plt.xlabel("Displacement from top of instrument / mm")
plt.ylabel("Bore cross-section / sq. mm")
plt.xlim((0, length))
plt.ylim((min(sdata) * 0.9, max(sdata) * 1.1))
plt.show()