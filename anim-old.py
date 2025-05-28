
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

iterations = (int) (4e4 // 1e2)
width = 128
step = 1
xmin = -2
xmax = 2

# initializing a figure in  
# which the graph will be plotted 
fig, (axp, axv) = plt.subplots(
    nrows=2,
    sharex=True,
    figsize=(4, 4),
)

pdata = np.zeros((iterations,width+2), np.float32)
vdata = np.zeros((iterations,width+1), np.float32)

with open("output/pdata.bin", "rb") as file:
    for i in range(iterations):
        for j in range(width+2):
            pdata[i][j] =  (np.float32) (struct.unpack("f", file.read(4)))[0]

with open("output/vdata.bin", "rb") as file:
    for i in range(iterations):
        for j in range(width+1):
            vdata[i][j] =  (np.float32) (struct.unpack("f", file.read(4)))[0]

assert type(axp) == plt.Axes and type(axv) == plt.Axes
axp.set_ylim(-1e4, 1e4)
axp.set_xlim(xmin - (xmax-xmin)/width, xmax + (xmax-xmin)/width)
axp.set_xlabel("x / m")
axp.set_ylabel("p' / Pa")
axp.grid()

axv.set_ylim(-2e1, 2e1)
axv.set_xlim(xmin - (xmax-xmin)/width, xmax + (xmax-xmin)/width)
axv.set_xlabel("x / m")
axv.set_ylabel("v / m/s")
axv.grid()

plt.subplots_adjust(left=0.25)

x = [xmin + (xmax-xmin)*(i-1.0)/width for i in range(width+2)]
xi = [xmin + (xmax-xmin)*(i-0.5)/width for i in range(width+1)]


# initializing a line variable 
linep, = axp.plot([], [], lw=2) 
linev, = axv.plot([], [], lw=2)
   
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
