import math
import numpy as np
import globalVariables as GV
global modelParams, TrackerParams, modelData

def modelData_ReadInFromFile(filename):
    f = open(filename, 'r')

    ZARRU = [   147.64,    457.24,    786.96,   1138.11,   1512.08,   1910.37,   
               2334.54,   2786.29,   3267.40,   3779.78,   4325.47,   4906.62,   
               5525.55,   6184.71,   6886.72,   7634.36,   8430.59,   9278.58,  
              10181.69,  11143.50,  12167.83,  13258.74,  14420.56,  15657.90,  
              16975.66,  18379.08,  19873.72,  21465.52,  23160.78,  24966.23]
    ZARRW = [   300.00,    619.50,    959.77,   1322.15,   1708.09,   2119.12,
               2556.86,   3023.06,   3519.56,   4048.33,   4611.47,   5211.21,
               5849.94,   6530.19,   7254.65,   8026.21,   8847.91,   9723.03,
              10655.02,  11647.60,  12704.70,  13830.50,  15029.49,  16306.41,
              17666.32,  19114.64,  20657.09,  22299.80,  24049.29,  25912.50]
    GV.ZARR = ZARRU
    
    NX = GV.NX
    NY = GV.NY
    NZ = GV.NZ
    DX = GV.modelParams['DX']
    DY = GV.modelParams['DY']
    
    URAW = np.zeros((NX, NY, NZ), float)
    VRAW = np.zeros((NX, NY, NZ), float)
    WRAW = np.zeros((NX, NY, NZ), float)
    
    for k in range(NZ):
        for i in range(NX):
            for j in range(NY):
                tempstr = f.readline()
                temparr = map(float, tempstr.split())
                GV.P[i, j, k] = temparr[0] * 100.0
                URAW[i, j, k] = temparr[1]
                VRAW[i, j, k] = temparr[2]
                WRAW[i, j, k] = temparr[3]
                tempstr = f.readline()
                temparr = map(float, tempstr.split())
                GV.THETA[i, j, k] = temparr[0]
                GV.QV[i, j, k] = temparr[1]
                GV.QL[i, j, k] = temparr[2]
    
    print 'Finished reading in ' + filename
    print 'Calculating other variables...'
    
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
                if ((i == 0) and (j != 0)):
                    GV.U[i, j, k] = URAW[i, j, k] - ((URAW[i+1, j, k] - URAW[i, j, k]) / 2.0)
                    GV.V[i, j, k] = (VRAW[i, j, k] + VRAW[i, j-1, k]) / 2.0
                elif ((i != 0) and (j == 0)):
                    GV.U[i, j, k] = (URAW[i, j, k] + URAW[i-1, j, k]) / 2.0
                    GV.V[i, j, k] = VRAW[i, j, k] - ((VRAW[i, j+1, k] - VRAW[i, j, k]) / 2.0)
                elif ((i == 0) and (j == 0)):
                    GV.U[i, j, k] = URAW[i, j, k] - ((URAW[i+1, j, k] - URAW[i, j, k]) / 2.0)
                    GV.V[i, j, k] = VRAW[i, j, k] - ((VRAW[i, j+1, k] - VRAW[i, j, k]) / 2.0)
                else:
                    GV.U[i, j, k] = (URAW[i, j, k] + URAW[i-1, j, k]) / 2.0
                    GV.V[i, j, k] = (VRAW[i, j, k] + VRAW[i, j-1, k]) / 2.0

    for k in range(NZ - 1):
        for i in range(NX):
            for j in range(NY):
                RELPOS = ZARRU[k] - ZARRW[k]
                WSLOPE = (WRAW[i, j, k+1] - WRAW[i, j, k]) / (ZARRW[k+1] - ZARRW[k])
                GV.W[i, j, k] = WRAW[i, j, k] + WSLOPE * RELPOS

    for i in range(NX):
        for j in range(NY):
            RELPOS = ZARRU[NZ-1] - ZARRW[NZ-1]
            WSLOPE = (WRAW[i, j, NZ-1] - WRAW[i, j, NZ-2]) / (ZARRW[NZ-1] - ZARRW[NZ-2])
            GV.W[i, j, NZ-1] = WRAW[i, j, NZ-1] + WSLOPE * RELPOS
#===============================================
# End of file read in
#===============================================

    for j in range(NY):
        for i in range(NX):
            for k in range(NZ):
            
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # Calculate derivatives... use raw variables when 
                # appropriate.
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (i == 0):
                    dudx = (URAW[i+1, j, k] - URAW[i, j, k]) / DX
                    dvdx = (GV.V[i+1, j, k] - GV.V[i, j, k]) / DX
                    dwdx = (GV.W[i+1, j, k] - GV.W[i, j, k]) / DX
                    dtdx = (GV.THETA[i+1, j, k] - GV.THETA[i, j, k]) / DX
                elif (i == NX-1):
                    dudx = (URAW[i, j, k] - URAW[i-1, j, k]) / DX
                    dvdx = (GV.V[i, j, k] - GV.V[i-1, j, k]) / DX
                    dwdx = (GV.W[i, j, k] - GV.W[i-1, j, k]) / DX
                    dtdx = (GV.THETA[i, j, k] - GV.THETA[i-1, j, k]) / DX
                else:
                    dudx = (URAW[i, j, k] - URAW[i-1, j, k]) / DX
                    dvdx = (GV.V[i+1, j, k] - GV.V[i-1, j, k]) / (2.0 * DX)
                    dwdx = (GV.W[i+1, j, k] - GV.W[i-1, j, k]) / (2.0 * DX)
                    dtdx = (GV.THETA[i+1, j, k] - GV.THETA[i-1, j, k]) / (2.0 * DX)
            
                if (j == 0):
                    dudy = (GV.U[i, j+1, k] - GV.U[i, j, k]) / DY
                    dvdy = (VRAW[i, j+1, k] - VRAW[i, j, k]) / DY
                    dwdy = (GV.W[i, j+1, k] - GV.W[i, j, k]) / DY
                    dtdy = (GV.THETA[i, j+1, k] - GV.THETA[i, j, k]) / DY
                elif (j == NY-1):
                    dudy = (GV.U[i, j, k] - GV.U[i, j-1, k]) / DY
                    dvdy = (VRAW[i, j, k] - VRAW[i, j-1, k]) / DY
                    dwdy = (GV.W[i, j, k] - GV.W[i, j-1, k]) / DY
                    dtdy = (GV.THETA[i, j, k] - GV.THETA[i, j-1, k]) / DY
                else:
                    dudy = (GV.U[i, j+1, k] - GV.U[i, j-1, k]) / (2.0 * DY)
                    dvdy = (VRAW[i, j+1, k] - VRAW[i, j, k]) / DY
                    dwdy = (GV.W[i, j+1, k] - GV.W[i, j-1, k]) / (2.0 * DY)
                    dtdy = (GV.THETA[i, j+1, k] - GV.THETA[i, j-1, k]) / (2.0 * DY)
            
                if (k == 0):
                    dudz = (GV.U[i, j, k+1] - GV.U[i, j, k]) / (ZARRU[k+1] - ZARRU[k])
                    dvdz = (GV.V[i, j, k+1] - GV.V[i, j, k]) / (ZARRU[k+1] - ZARRU[k])
                    #dwdz = (WRAW[i, j, k+1] - WRAW[i, j, k]) / (ZARRW[k+1] - ZARRW[k])
                    dtdz = (GV.THETA[i, j, k+1] - GV.THETA[i, j, k]) / (ZARRU[k+1] - ZARRU[k])
                elif (k == NZ-1):
                    dudz = (GV.U[i, j, k] - GV.U[i, j, k-1]) / (ZARRU[k] - ZARRU[k-1])
                    dvdz = (GV.V[i, j, k] - GV.V[i, j, k-1]) / (ZARRU[k] - ZARRU[k-1])
                    #dwdz = (WRAW[i, j, k] - WRAW[i, j, k-1]) / (ZARRW[k] - ZARRW[k-1])
                    dtdz = (GV.THETA[i, j, k] - GV.THETA[i, j, k-1]) / (ZARRU[k] - ZARRU[k-1])
                else:
                    dudz = (GV.U[i, j, k+1] - GV.U[i, j, k-1]) / (ZARRU[k+1] - ZARRU[k-1])
                    dvdz = (GV.V[i, j, k+1] - GV.V[i, j, k-1]) / (ZARRU[k+1] - ZARRU[k-1])
                    #dwdz = (WRAW[i, j, k+1] - WRAW[i, j, k]) / (ZARRW[k+1] - ZARRW[k])
                    dtdz = (GV.THETA[i, j, k+1] - GV.THETA[i, j, k-1]) / (ZARRU[k+1] - ZARRU[k-1])
            
                GV.DIVERG[i, j, k] = dudx + dvdy
            
                GV.VORT[i, j, k] = dvdx - dudy
            
                GV.AVORT[i, j, k] = GV.VORT[i, j, k] + GV.modelParams['f0']
            
                GV.SPEED[i, j, k] = ((GV.U[i, j, k] ** 2.0) + (GV.V[i, j, k] ** 2.0)) ** 0.5
            
                GV.TEMPER[i, j, k] = GV.THETA[i, j, k] * ((GV.P[i, j, k] / 100000.0) ** 0.286)
            
                GV.TEMPVT[i, j, k] = GV.TEMPER[i, j, k] * (1.0 + GV.QV[i, j, k] * 0.607717e-3)
            
                GV.RHO[i, j, k] = GV.P[i, j, k] / (287 * GV.TEMPER[i, j, k])
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate theta-e, Bolton (1980)
#
# mixratio = [g/g]
# eVapor   = [mb]

                dQV = GV.QV[i, j, k]
                dP  = GV.P[i, j, k]
                dTK = GV.TEMPER[i, j, k]
            
                mixratio = (dQV / 1.0E3) / (1.0E0 - (dQV / 1.0E3))
                if (mixratio == 0):
                    mixratio = 1.0E-11

                eVapor   = (dP / 1.0E2) * mixratio / (mixratio + 0.622E0)
            
                tempf1   = 0.2854E0 * (1.0E0 - 0.28E0 * mixratio)
          
                tempf2   = 2840.0E0 / (3.5E0 * math.log(dTK) - math.log(eVapor) - 4.805E0) + 55.0E0
            
                tempf3   = math.exp(1.0E3 * mixratio * (1.0E0 - 0.81E0 * mixratio) * \
                           ((3.376E0 / tempf2) - 0.00254E0))
            
                GV.THETAE[i, j, k] = dTK * tempf3 * (1.0E5 / dP) ** tempf1
            

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Calculate Ertel's PV, full version.
            
                GV.PV[i, j, k] = ((dwdy - dvdz) * dtdx + (dudz - dwdx) * dtdy + \
                                  (dvdx - dudy + GV.modelParams['f0']) * dtdz) / GV.RHO[i, j, k]
            
                tempf1      = 6.112 * math.exp(17.67 * (GV.TEMPER[i, j, k] - 273.15) / (GV.TEMPER[i, j, k] - 29.65))
                GV.RH[i, j, k] = (GV.QV[i, j, k] / 1000.0) / (0.622 * tempf1 / ((GV.P[i, j, k] / 100.0) - tempf1))
                if (GV.RH[i, j, k] > 1.0):
                    GV.RH[i, j, k] = 1.0
            
    print 'Completed calculations on file read-in...'
    
