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

    data['ThetaEIn'] = GV.THETAE[x, y, z]

    #------------------------------------------
    # 2-box calculations
    #------------------------------------------
    xRange = boxBounds(x, 0, GV.modelParams['NX']-1, 2)
    minX = xRange[0]
    maxX = xRange[1]

    yRange = boxBounds(y, 0, GV.modelParams['NY']-1, 2)
    minY = yRange[0]
    maxY = yRange[1]

    totalPoints = 2 * (maxX - minX) + 2 * (maxY - minY)

    data['ThetaEBox2'] = (sum(GV.THETAE[minX:maxX+1, minY, z]) + \
                          sum(GV.THETAE[minX:maxX+1, maxY, z]) + \
                          sum(GV.THETAE[minX, minY:maxY+1, z]) + \
                          sum(GV.THETAE[maxX, minY:maxY+1, z])) / totalPoints

    #------------------------------------------
    # 4-box calculations
    #------------------------------------------
    xRange = boxBounds(x, 0, GV.modelParams['NX']-1, 4)
    minX = xRange[0]
    maxX = xRange[1]

    yRange = boxBounds(y, 0, GV.modelParams['NY']-1, 4)
    minY = yRange[0]
    maxY = yRange[1]

    totalPoints = 2 * (maxX - minX) + 2 * (maxY - minY)

    data['ThetaEBox4'] = (sum(GV.THETAE[minX:maxX+1, minY, z]) + \
                          sum(GV.THETAE[minX:maxX+1, maxY, z]) + \
                          sum(GV.THETAE[minX, minY:maxY+1, z]) + \
                          sum(GV.THETAE[maxX, minY:maxY+1, z])) / totalPoints

    #------------------------------------------
    # 3-box calculations
    #------------------------------------------
    xRange = boxBounds(x, 0, GV.modelParams['NX']-1, 3)
    minX = xRange[0]
    maxX = xRange[1]

    yRange = boxBounds(y, 0, GV.modelParams['NY']-1, 3)
    minY = yRange[0]
    maxY = yRange[1]

    totalPoints = 2 * (maxX - minX) + 2 * (maxY - minY)

    data['ThetaEBox3'] = (sum(GV.THETAE[minX:maxX+1, minY, z]) + \
                          sum(GV.THETAE[minX:maxX+1, maxY, z]) + \
                          sum(GV.THETAE[minX, minY:maxY+1, z]) + \
                          sum(GV.THETAE[maxX, minY:maxY+1, z])) / totalPoints

    data['RHBox3'] =     (sum(GV.RH[minX:maxX+1, minY, z]) + \
                          sum(GV.RH[minX:maxX+1, maxY, z]) + \
                          sum(GV.RH[minX, minY:maxY+1, z]) + \
                          sum(GV.RH[maxX, minY:maxY+1, z])) / totalPoints

    data['TempBox3'] =   (sum(GV.TEMPER[minX:maxX+1, minY, z]) + \
                          sum(GV.TEMPER[minX:maxX+1, maxY, z]) + \
                          sum(GV.TEMPER[minX, minY:maxY+1, z]) + \
                          sum(GV.TEMPER[maxX, minY:maxY+1, z])) / totalPoints

    #---------------------------------------
    # Other data calculations
    #---------------------------------------

    data['NearCoreIntQl'] = np.mean(GV.QL[minX:maxX+1, minY:maxY+1, z])


    return data
