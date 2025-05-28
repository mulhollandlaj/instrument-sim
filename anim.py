
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation  
import struct
import os

def uniquify(path):
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + " (" + str(counter) + ")" + extension
        counter += 1

    return path

frames = 300
iterations = 1000
width = 128
step = (int) (iterations / frames)
xmin = -2
xmax = 2

# initializing a figure in  
# which the graph will be plotted 
fig, (axp, axv) = plt.subplots(
    nrows=2,
    sharex=True,
    figsize=(4, 4),
)


maxp = -1e6
minp = 1e6
with open("output/pdata.bin", "rb") as file:
    nt = int(struct.unpack("i", file.read(4))[0])
    width = int(struct.unpack("i", file.read(4))[0])
    
    pdata = np.zeros((iterations,width), np.float32)
    vdata = np.zeros((iterations,width), np.float32)

    for i in range(iterations):
        for j in range(width):
            pdata[i][j] = (np.float32) (struct.unpack("f", file.read(4)))[0]
            maxp = max(pdata[i][j], maxp)
            minp = min(pdata[i][j], minp)

maxv = -1e6
minv = 1e6
with open("output/vdata.bin", "rb") as file:
    file.read(8)
    for i in range(iterations):
        for j in range(width):
            vdata[i][j] =  (np.float32) (struct.unpack("f", file.read(4)))[0]
            maxv = max(vdata[i][j], maxv)
            minv = min(vdata[i][j], minv)

assert type(axp) == plt.Axes and type(axv) == plt.Axes
axp.set_xlim(xmin - (xmax-xmin)/width, xmax + (xmax-xmin)/width)
axp.set_ylim(-2e4, 2e4)
axv.set_ylim(-2e2, 2e2)

x = [xmin + (xmax-xmin)*(i-1.0)/width for i in range(width)]
xi = [xmin + (xmax-xmin)*(i-0.5)/width for i in range(width)]

# initializing a line variable 
linep, = axp.plot([], [], lw=1) 
linev, = axv.plot([], [], lw=1)
   
# data which the line will
# contain (x, y) 
def init():
    linep.set_data([], [])
    linev.set_data([], [])
    return linep, linev,
   
def animate(i): 
    linep.set_data(x, pdata[i*step])
    linev.set_data(xi, vdata[i*step])
    
    return linep, linev, 
    
anim = FuncAnimation(fig, animate, init_func = init, 
                     frames = iterations // step, interval = 20, blit = False)

plt.show()

if (input("Do you want to save this plot? y/n ") == "y"):
    animpath = uniquify('output/anim.mp4')
    anim.save(animpath, writer = 'ffmpeg', fps = 30) 
