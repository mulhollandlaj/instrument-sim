from matplotlib import pyplot as plt 
import numpy as np
import csv
import struct

data = []

fname_short = "delusse_oboe"

with open("input\\"+fname_short+".csv", "r") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        data.append((float(row[0]), float(row[1]), float(row[2])))

npdata = np.array(data, dtype=np.float32).transpose()

nx = 128
length = npdata[0][-1]

h = length / float(nx+1)


with open("input\\"+fname_short+".bin", "wb") as outfile:
    outfile.write(struct.pack('<i',len(npdata[0])))
    outfile.write(struct.pack('<f',npdata[2][-1]))
    tdata = tuple(np.concat((npdata[0], npdata[1], npdata[2])))
    for f in tdata:
        outfile.write(struct.pack('<f', f))


# plt.plot(npdata[1], npdata[0],'r--')
# plt.plot(npdata[2], npdata[0],'k:')
# plt.grid(which='both')
# plt.xlabel("Displacement from top of instrument / mm")
# plt.ylabel("Bore diameter / mm")
# plt.xlim((0, max(npdata[1])))
# plt.ylim((0, max(npdata[0])))


# plt.show()