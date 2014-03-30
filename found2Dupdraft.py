import math
import numpy as np
import globalVariables as GV
global modelParams, TrackerParams, modelData

def boxBounds(startValue, minimum, maximum, boxSize):
    if (startValue < minimum + boxSize):
        minValue = minimum
    else:
        minValue = startValue - boxSize
    if (startValue > maximum - boxSize):
        maxValue = maximum
    else:
        maxValue = startValue + boxSize
    return [minValue, maxValue]

#======================================================================
# found2Dupdraft
#
# Function that takes the point of a found 2D updraft and returns a 
# dictionary of calculations pertaining to that point.
#
# This function should be updated to reflect the user's interests

def found2Dupdraft(x, y, z):
    data = {}

    nearUMean = (GV.U[x-1, y, z] + GV.U[x+1, y, z] + GV.U[x, y-1, z] + GV.U[x, y+1, z] + \
                 GV.U[x-1, y-1, z] + GV.U[x-1, y+1, z] + GV.U[x+1, y-1, z] + GV.U[x+1, y+1, z]) / 8.0
    data['nearUMean'] = nearUMean

    nearVMean = (GV.V[x-1, y, z] + GV.V[x+1, y, z] + GV.V[x, y-1, z] + GV.V[x, y+1, z] + \
                 GV.V[x-1, y-1, z] + GV.V[x-1, y+1, z] + GV.V[x+1, y-1, z] + GV.V[x+1, y+1, z]) / 8.0
    data['nearVMean'] = nearVMean

    data['WInside'] = GV.W[x, y, z]

    #------------------------------------------
    # 3-box calculation
    #------------------------------------------
    xRange = boxBounds(x, 0, GV.modelParams['NX']-1, 3)
    minX = xRange[0]
    maxX = xRange[1]

    yRange = boxBounds(y, 0, GV.modelParams['NY']-1, 3)
    minY = yRange[0]
    maxY = yRange[1]

    totalPoints = 2 * (maxX - minX) + 2 * (maxY - minY)

    data['ThetaBox3'] =  (sum(GV.THETA[minX:maxX+1, minY, z]) + \
                          sum(GV.THETA[minX:maxX+1, maxY, z]) + \
                          sum(GV.THETA[minX, minY:maxY+1, z]) + \
                          sum(GV.THETA[maxX, minY:maxY+1, z])) / totalPoints

    #---------------------------------------
    # Other data calculations
    #---------------------------------------

    data['NearCoreIntQV'] = np.mean(GV.QV[minX:maxX+1, minY:maxY+1, z])


    return data
