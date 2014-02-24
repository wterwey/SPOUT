import math, os, pickle
import numpy as np
import matplotlib.pyplot as plt
import updraftClasses as up
import modelReadIn as modelread
import globalVariables as GV
import found2Dupdraft as f2D
from matplotlib import cm

global modelParams, TrackerParams, modelData

updrafts = []

fileList = os.listdir(GV.modelParams['FILEROOT'])
fileList = sorted(fileList)
t = -1

outputfile = open('finalUpdrafts.dat', 'w')

idNum = 0

for currFile in fileList:
    t = t + 1
    modelread.modelData_ReadInFromFile(GV.modelParams['FILEROOT']+'/'+currFile)

    print 'Establishing current time updrafts...'
    
    curr3D = []

    #-----------------------------------------------------
    # Begin searching through the W array for 2D updrafts,
    # then look for connections to make the collection of
    # 3D updrafts.
    #-----------------------------------------------------

    for m in reversed(GV.TrackerParams['checkLevels']):
        testArray = GV.W[:, :, m]
        mVal = m - GV.TrackerParams['MINUPHEIGHT'] + 1
        checkVal = GV.TrackerParams['checkValues'][mVal]
        wTh = GV.TrackerParams['WTHRES']

        for i in range(2, GV.modelParams['NX']-2):
            for j in range(2, GV.modelParams['NY']-2):
                if (testArray[i, j] >= checkVal):
                    tempTotal = 0
                    for x in range(-1, 2):
                        for y in range(-1, 2):
                            if (testArray[i+x, j+y] > checkVal):
                                tempTotal = tempTotal + 1
                    if (tempTotal >= 5):
                        if ((testArray[i, j] >= testArray[i-2, j] + wTh) and \
                            (testArray[i, j] >= testArray[i-1, j] + wTh) and \
                            (testArray[i, j] >= testArray[i+1, j] + wTh) and \
                            (testArray[i, j] >= testArray[i+2, j] + wTh) and \
                            (testArray[i, j] >= testArray[i-1, j-1] + wTh) and \
                            (testArray[i, j] >= testArray[i, j-1] + wTh) and \
                            (testArray[i, j] >= testArray[i+1, j-1] + wTh) and \
                            (testArray[i, j] >= testArray[i, j-2] + wTh) and \
                            (testArray[i, j] >= testArray[i-1, j+1] + wTh) and \
                            (testArray[i, j] >= testArray[i, j+1] + wTh) and \
                            (testArray[i, j] >= testArray[i+1, j+1] + wTh) and \
                            (testArray[i, j] >= testArray[i, j+2] + wTh)):
                            
                            tempData = f2D.found2Dupdraft(i, j, m)
                            tempDraft = up.updraft2D(i, j, m, t, tempData)
                            
                            if (m == max(GV.TrackerParams['checkLevels'])):
                                curr3D.append(up.updraft3D([tempDraft], m))
                            # If at the top level, start new 3D updrafts; if
                            # not, look for connections above.
                            else:
                                connection = 0
                                tempRank   = 1000.0
                                for upTemp3 in curr3D:
                                    if (((upTemp3.vertRange())[0] == m + 1)):
                                        for upTemp2 in upTemp3.linked2DUpdrafts:
                                            if (upTemp2.zPos == m+1):
                                                tempRangeX = abs(i - upTemp2.xPos)
                                                tempRangeY = abs(j - upTemp2.yPos)
                                                if ((tempRangeX <= GV.TrackerParams['SPACETHRES']) and \
                                                    (tempRangeY <= GV.TrackerParams['SPACETHRES'])):
                                                    
                                                    currRank = (tempRangeX ** 2.0 + tempRangeY ** 2.0) ** 0.5
                                                    if (currRank < tempRank):
                                                        tempRank = currRank
                                                        temp3Updraft = upTemp3
                                                        connection = 1

                                                # if within SPACETHRES, add this point to the 3D updraft, if
                                                # none already matched up.
                                if (connection == 0):
                                    curr3D.append(up.updraft3D([tempDraft], m))
                                # If not connected, start a new 3D updraft
                                else:
                                    temp3Updraft.linked2DUpdrafts.append(tempDraft)
                                    temp3Updraft.bottomHeight = m

    # Complete the 3D updraft finder for the current time.
    newCurr3D = []
    for upTemp3 in curr3D:
        tempRange = upTemp3.vertRange()
        if (tempRange[1]-tempRange[0]+1 >= GV.TrackerParams['MINUPHEIGHT']):
            newCurr3D.append(upTemp3)
    
    if (t == 0):
        for temp3D in newCurr3D:
            updrafts.append(up.updraft4D([temp3D], t, idNum))
            idNum = idNum + 1
    else:
        for temp3D in newCurr3D:
            minRank1 = 1000.0
            minRank2 = 1001.0
            minRank3 = 1002.0
            minUp1   = []
            minUp2   = []
            minUp3   = []
            minPts1  = 0
            minPts2  = 0
            minPts3  = 0
            for temp4D in updrafts:
                up3Dprior = temp4D.findTime(t-1)
                up3Dcurr  = temp4D.findTime(t)

                if ((up3Dcurr == -1) and (up3Dprior != -1)):
                    priorRange = up3Dprior.vertRange()
                    totalUadv = 0.0
                    totalVadv = 0.0
                    tempCount = 0
                    for temp2D in up3Dprior.linked2DUpdrafts:
                        totalUadv = totalUadv + temp2D.savedData['nearUMean']
                        totalVadv = totalVadv + temp2D.savedData['nearVMean']
                        tempCount = tempCount + 1
                    avgUadv = totalUadv / tempCount
                    avgVadv = totalVadv / tempCount

                    totalXadv = 0.0
                    totalYadv = 0.0
                    tempCount = 0
                    gridCorrX = GV.modelParams['DT'] * 3600.0 / GV.modelParams['DX']
                    gridCorrY = GV.modelParams['DT'] * 3600.0 / GV.modelParams['DY']
                    for temp2D in up3Dprior.linked2DUpdrafts:
                        totalXadv = (temp2D.xPos + GV.TrackerParams['velCorrect'] * avgUadv * gridCorrX) + \
                                    totalXadv
                        totalYadv = (temp2D.yPos + GV.TrackerParams['velCorrect'] * avgVadv * gridCorrY) + \
                                    totalYadv
                        tempCount = tempCount + 1
                    advectedX = totalXadv / tempCount
                    advectedY = totalYadv / tempCount

                    totalXcur = 0.0
                    totalYcur = 0.0
                    tempCount = 0
                    for temp2D in temp3D.linked2DUpdrafts:
                        totalXcur = temp2D.xPos + totalXcur
                        totalYcur = temp2D.yPos + totalYcur
                        tempCount = tempCount + 1
                    currentX = totalXcur / tempCount
                    currentY = totalYcur / tempCount

                    checkRank = ((currentX - advectedX) ** 2 + (currentY - advectedY) ** 2) ** 0.5

                    currentRange = temp3D.vertRange()
                    numRankPts = min([currentRange[1], priorRange[1]]) - \
                                 max([currentRange[0], priorRange[0]]) + 1
                    if (numRankPts <= 0):
                        numRankPts = 0

                    if (checkRank < minRank1):
                        minRank3 = minRank2
                        minRank2 = minRank1
                        minRank1 = checkRank
                        minUp3   = minUp2
                        minUp2   = minUp1
                        minUp1   = temp4D
                        minPts3  = minPts2
                        minPts2  = minPts1
                        minPts1  = numRankPts
                    elif (checkRank < minRank2):
                        minRank3 = minRank2
                        minRank2 = checkRank
                        minUp3   = minUp2
                        minUp2   = temp4D
                        minPts3  = minPts2
                        minPts2  = numRankPts
                    elif (checkRank < minRank3):
                        minRank3 = checkRank
                        minUp3   = temp4D
                        minPts3  = numRankPts

            #print (minRank1, minPts1), (minRank2, minPts2), (minRank3, minPts3)
            if ((minRank1 <= GV.TrackerParams['RANKTHRES']) and \
                (minPts1 >= GV.TrackerParams['POINTTHRES'])):
                minUp1.linked3DUpdrafts.append(temp3D)
                minUp1.endTime = t
            elif ((minRank2 <= GV.TrackerParams['RANKTHRES']) and \
                  (minPts2 >= GV.TrackerParams['POINTTHRES'])):
                minUp2.linked3DUpdrafts.append(temp3D)
                minUp2.endTime = t
            elif ((minRank3 <= GV.TrackerParams['RANKTHRES']) and \
                  (minPts3 >= GV.TrackerParams['POINTTHRES'])):
                minUp3.linked3DUpdrafts.append(temp3D)
                minUp3.endTime = t
            else:
                updrafts.append(up.updraft4D([temp3D], t, idNum))
                idNum = idNum + 1

    currUpdrafts = []
    for temp4D in updrafts:
        if (temp4D.endTime != t):
            if (temp4D.endTime - temp4D.beginTime + 1 >= GV.TrackerParams['OUTPUTTIMETHRES']):
                pickle.dump(temp4D, outputfile)
        else:
            currUpdrafts.append(temp4D)
    updrafts = currUpdrafts

#-----------------------------------------------
# End time/file loop...
#-----------------------------------------------

for temp4D in updrafts:
    if(temp4D.endTime - temp4D.beginTime + 1 >= GV.TrackerParams['OUTPUTTIMETHRES']):
        pickle.dump(temp4D, outputfile)

outputfile.close()
