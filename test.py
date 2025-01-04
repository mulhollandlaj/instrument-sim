import numpy as np
import matplotlib.pyplot as plt
import struct

iterations = 1000
width = 512

plt.style.use('_mpl-gallery')

pdata = np.zeros((iterations,width), np.float32)

with open("output/pdata.bin", "rb") as file:
    for i in range(iterations):
        for j in range(width):
            pdata[i][j] = (np.float32) (struct.unpack("f", file.read(4)))


# plot
fig, ax = plt.subplots()
ax.imshow(pdata.transpose(), origin='lower')


plt.savefig("output/figure")