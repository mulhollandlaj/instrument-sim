import numpy as np
import matplotlib.pyplot as plt
import struct


plt.style.use('_mpl-gallery')

pdata = np.zeros((1000,128), np.float32)

with open("src/sim/pdata.bin", "rb") as file:
    for i in range(1000):
        for j in range(128):
            pdata[i][j] = (np.float32) (struct.unpack("f", file.read(4)))


# plot
fig, ax = plt.subplots()
ax.imshow(pdata.transpose(), origin='lower')


plt.show()

input()

breakpoint