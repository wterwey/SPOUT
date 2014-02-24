import math, os, pickle
import numpy as np
import matplotlib.pyplot as plt
import updraftClasses as up
import modelReadIn as modelread
import globalVariables as GV
import found2Dupdraft as f2D
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

global modelParams, TrackerParams, modelData

inputfile = open('finalUpdrafts.dat', 'r')
updrafts  = []

print 'Loading updrafts...'

try:
    while True:
        data = pickle.load(inputfile)
        updrafts.append(data)
except EOFError:
    print 'Done reading in updrafts...'

inputfile.close()

for t in range(240):
    fig2D = plt.figure()
    fig3D = plt.figure()
    ax2D = fig2D.add_subplot(111)
    ax   = fig3D.add_subplot(111, projection='3d')
    ax.set_xlim3d(0, 250)
    ax.set_ylim3d(0, 250)
    ax.set_zlim3d(0, 30)
    for temp4D in updrafts:
        temp3D = temp4D.findTime(t)
        if (temp3D != -1):
            x = []
            y = []
            z = []
            for temp2D in temp3D.linked2DUpdrafts:
                x.append(temp2D.xPos)
                y.append(temp2D.yPos)
                z.append(temp2D.zPos)

            ax.plot(x, y, z, color="black")
            ax.text(x[0], y[0], z[0]+1, str(temp4D.ident))

            xAvg = np.mean(x)
            yAvg = np.mean(y)
            ax2D.set_xlim(0, 250)
            ax2D.set_ylim(0, 250)
            ax2D.scatter(xAvg, yAvg, s=30, color="black")
            ax2D.text(xAvg, yAvg+5, str(temp4D.ident), fontsize=6)

    tStr = str(t)
    figName = 'currentUpdrafts' + tStr + '.png'
    fig3D.savefig(figName)
    figName = 'currentUpdrafts2D' + tStr + '.png'
    fig2D.savefig(figName)


