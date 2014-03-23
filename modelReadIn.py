import netCDF4
import math
import numpy as np
import globalVariables as GV
global modelParams, TrackerParams, modelData

#==================================
# Basic interpolation function
#==================================
def interp2hgt(zlo, zhi, varlo, varhi, zcur):
    return (varhi - varlo) * (zcur - zlo) / (zhi - zlo) + varlo


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For users of WRF, a standard WRF read-in function, assuming netCDF files is
# included below.
#
# SPOUT assumes everything on an Arakawa A-grid in height coordinates, so 
# interpolation to this grid is done here.  The zLevels array specifies the 
# vertical level interpolation points.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def modelData_ReadInFromFile_WRF(filename):
    zgrid = [   50.0, 107.307, 205.706, 334.344, 496.915, 
	     701.513, 953.557, 1254.76, 1617.94, 2035.14, 
	     2504.95, 3027.17, 3605.44, 4247.20, 4941.06, 
	     5715.49, 6492.84, 7305.57, 8122.98, 8914.02, 
	     9741.23, 10578.8, 11418.7, 12260.2, 13107.5, 
	     13977.3, 14843.9, 15830.8, 16832.0, 17951.9, 
	     19213.7, 20862.0, 22934.2, 24000.0]
    nz2 = len(zgrid)
    
    fileRead = netCDF4.Dataset(filename, 'r')

    nx     = fileRead.variables['west_east'][:]
    ny     = fileRead.variables['south_north'][:]
    nz     = fileRead.variables['bottom_top'][:]
    nxp1   = fileRead.variables['west_east_stag'][:]
    nyp1   = fileRead.variables['south_north_stag'][:]
    nzp1   = fileRead.variables['bottom_top_stag'][:]
    
    print nx, ny, nz, nyp1, nzp1

    ph     = fileRead.variables['PH'][:]
    phb    = fileRead.variables['PHB'][:]
    theta  = fileRead.variables['T'][:]
    p      = fileRead.variables['P'][:]
    pb     = fileRead.variables['PB'][:]
    qvapor = fileRead.variables['QVAPOR'][:]
    u      = fileRead.variables['U'][:]
    v      = fileRead.variables['V'][:]
    w      = fileRead.variables['W'][:]
    znu    = fileRead.variables['ZNU'][:]
    znw    = fileRead.variables['ZNW'][:]
    mu     = fileRead.variables['MU'][:]
    p_top  = fileRead.variables['P_TOP'][:]
    
    p = p + pb
    ph = ph + phb
    del pb
    del phb
    
    p0 = 1e+05
    g  = 9.81e+00
    cp = 1.0046e+03
    rd = 2.87e+02
    kappa = rd / cp

    theta = theta + 300
    rho = (p0 ** kappa /rd) * p ** (1.0e+00 - kappa) / (theta * \
	(1.0e+00 + 1.61e+00 * qvapor) )
    temp = theta * (p / p0) ** kappa


    print 'Calculating geopotential heights at regular grid points'
    phup = np.zeros((nx, ny, nz), float)
    phdown = np.zeros((nx, ny, nz), float)
    phit = np.zeros((nx, ny, nz), float)

    phup[:, :, 0] = ph[:, :, 0] - \
	rd * (3.0e+00 * temp[:, :, 0] - \
	1.5e+00 * (temp[:, :, 0] + temp[:, :, 1]) + \
	1.0e+00 * temp[:, :, 1] ) * \
	math.log( (mu * znu[0] + p_top) / (mu * znw[0] + p_top) )
    phdown[:, :, 0] = ph[:, :, 1] - \
	rd * 0.5 * (temp[:, :, 0] + temp[:, :, 1]) * \
	math.log( (mu * znu[0] + p_top) / (mu * znw[1] + p_top) )
    phit[:, :, 0] = 0.5 * (phup[:, :, 0] + phdown[:, :, 0])
    
    for k in np.arange(1, nz-1, 1):
        phup[:, :, k] = ph[:, :, k] - \
		rd * 0.5 * (temp[:, :, k] + temp[:, :, k-1]) * \
		math.log( (mu * znu[k] + p_top) / (mu * znw[k] + p_top) )
	phdown[:, :, k] = ph[:, :, k+1] - \
        	rd * 0.5 * (temp[:, :, k] + temp[:, :, k+1]) * \
		math.log( (mu * znu[k] + p_top) / (mu * znw[k+1] + p_top) )
	phit[:, :, k] = 0.5 * (phup[:, :, k] + phdown[:, :, k])

    phup[:, :, nz-1] = ph[:, :, nz-1] - \
	rd * 0.5 * (temp[:, :, nz-1] + temp[:, :, nz-2]) * \
	math.log( (mu * znu[nz-1] + p_top) / (mu * znw[nz-1] + p_top) )
    phdown[:, :, nz-1] = ph[:, :, nzp1-1] - \
	rd * (3.0e+00 * temp[:, :, nz-1] - \
	1.5e+00 * (temp[:, :, nz-1] + temp[:, :, nz-2]) + \
	1.0e+00 * temp[:, :, nz-2] ) * \
	math.log( (mu * znu[nz-1] + p_top) / (mu * znw[nzp1-1] + p_top) )
    phit[:, :, nz-1] = 0.5 * (phup[:, :, nz-1] + phdown[:, :, nz-1])

    del phup, phdown
    

    print "Interpolating velocity components..."
    u2 = np.zeros((nx, ny, nz), float)
    v2 = np.zeros((nx, ny, nz), float)
    w2 = np.zeros((nx, ny, nz), float)

    for i in np.arange(0, nx-1, 1):
        u2[i, :, :] = 0.5 * (u[i, :, :] + u[i+1, :, :])
        v2[:, i, :] = 0.5 * (v[:, i, :] + v[:, i+1, :])
        
    for i in np.arange(0, nz-1, 1):
        w2[:, :, i] = 0.5 * (w[:, :, i] + w[:, :, i+1])


    print 'Interpolating variables to physical height grid'

    ugd = np.zeros((nx, ny, nz2-2), float)
    

    for k in np.arange(1, nz2-2, 1):
        for i in np.arange(0, nx, 1):
            for j in np.arange(0, ny, 1):
                dzhi = 9.9e+09
                dzlo = 9.9e+09
                zhi  = -9.9e+09
                zlo  = -9.9e+09
                for l in np.arange(0, nz, 1):
                    zlev = phit[i, j, l] / g
                    if ((zgrid[k] - zlev) >= 0.0 and (zgrid[k] - zlev) < dzlo):
                        dzlo = zgrid[k] - zlev
                        zlo = zlev
                        klo = 1
                    if ((zlev - zgrid[k]) > 0.0 and (zlev - zgrid[k]) < dzhi):
                        dzhi = zgrid[k] - zlev
                        zhi = zlev
                        khi = 1
                
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read in the model data from plain text (ASCII) files
#
# This code will need to be changed to fit with the user's dataset and how
# that data is set up in their files.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def modelData_ReadInFromFile_PlainText(filename):
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
    
