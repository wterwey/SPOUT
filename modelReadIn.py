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

#    nx     = fileRead.dimensions['west_east'][:]
#    ny     = fileRead.dimensions['south_north'][:]
#    nz     = fileRead.dimensions['bottom_top'][:]
#    nxp1   = fileRead.dimensions['west_east_stag'][:]
#    nyp1   = fileRead.dimensions['south_north_stag'][:]
#    nzp1   = fileRead.dimensions['bottom_top_stag'][:]
    dx     = 1.0 / fileRead.variables['RDX'][:][0]
    dy     = 1.0 / fileRead.variables['RDY'][:][0]

    ph     = fileRead.variables['PH'][:]
    f0     = fileRead.variables['F'][:]

    nz = ph.shape[1]
    nx = ph.shape[2]
    ny = ph.shape[3]
    nzp1 = nz + 1
    
    ph     = fileRead.variables['PH'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    phb    = fileRead.variables['PHB'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    theta  = fileRead.variables['T'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    p      = fileRead.variables['P'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    pb     = fileRead.variables['PB'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    qvapor = fileRead.variables['QVAPOR'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    u      = fileRead.variables['U'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    v      = fileRead.variables['V'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    w      = fileRead.variables['W'][:].transpose((2, 3, 1, 0))[:, :, :, 0]
    znu    = fileRead.variables['ZNU'][:][0,:]
    znw    = fileRead.variables['ZNW'][:][0,:]
    mu     = fileRead.variables['MU'][:][0,:,:]
    mub    = fileRead.variables['MUB'][:][0,:,:]
    p_top  = fileRead.variables['P_TOP'][:][0]
    
    p = p + pb
    ph = ph + phb
    mu = mu + mub
    del pb
    del phb
    del mub
    
    p0 = 1e+05
    g  = 9.81e+00
    cp = 1.0046e+03
    rd = 2.87e+02
    kappa = rd / cp

    theta = theta + 300
    rho = (p0 ** kappa /rd) * p ** (1.0e+00 - kappa) / (theta * \
	(1.0e+00 + 1.61e+00 * qvapor) )
    temp = theta * np.power((p / p0), kappa)


    print 'Calculating geopotential heights at regular grid points'
    phup = np.zeros((nx, ny, nz), float)
    phdown = np.zeros((nx, ny, nz), float)
    phit = np.zeros((nx, ny, nz), float)

    for i in range(nx):
        for j in range(ny):
            phup[i, j, 0] = ph[i, j, 0] - \
   	       rd * (3.0e+00 * temp[i, j, 0] - \
   	       1.5e+00 * (temp[i, j, 0] + temp[i, j, 1]) + \
	       1.0e+00 * temp[i, j, 1] ) * \
	       math.log( (mu[i, j] * znu[0] + p_top) / (mu[i, j] * znw[0] + p_top) )
            phdown[i, j, 0] = ph[i, j, 1] - \
	       rd * 0.5 * (temp[i, j, 0] + temp[i, j, 1]) * \
	       math.log( (mu[i, j] * znu[0] + p_top) / (mu[i, j] * znw[1] + p_top) )
            phit[i, j, 0] = 0.5 * (phup[i, j, 0] + phdown[i, j, 0])
    
            for k in np.arange(1, nz-3, 1):
                phup[i, j, k] = ph[i, j, k] - \
	           rd * 0.5 * (temp[i, j, k] + temp[i, j, k-1]) * \
		   math.log( (mu[i, j] * znu[k] + p_top) / (mu[i, j] * znw[k] + p_top) )
	        phdown[i, j, k] = ph[i, j, k+1] - \
        	   rd * 0.5 * (temp[i, j, k] + temp[i, j, k+1]) * \
		   math.log( (mu[i, j] * znu[k] + p_top) / (mu[i, j] * znw[k+1] + p_top) )
	        phit[i, j, k] = 0.5 * (phup[i, j, k] + phdown[i, j, k])

            phup[i, j, nz-2] = ph[i, j, nz-2] - \
	       rd * 0.5 * (temp[i, j, nz-2] + temp[i, j, nz-2]) * \
	       math.log( (mu[i, j] * znu[nz-2] + p_top) / (mu[i, j] * znw[nz-2] + p_top) )
            phdown[i, j, nz-2] = ph[i, j, nzp1-2] - \
	       rd * (3.0e+00 * temp[i, j, nz-2] - \
	       1.5e+00 * (temp[i, j, nz-2] + temp[i, j, nz-2]) + \
	       1.0e+00 * temp[i, j, nz-2] ) * \
	       math.log( (mu[i, j] * znu[nz-2] + p_top) / (mu[i, j] * znw[nzp1-2] + p_top) )
            phit[i, j, nz-2] = 0.5 * (phup[i, j, nz-2] + phdown[i, j, nz-2])

    del phup, phdown
    

    print "Interpolating velocity components..."
    u2 = np.zeros((nx, ny, nz), float)
    v2 = np.zeros((nx, ny, nz), float)
    w2 = np.zeros((nx, ny, nz), float)

    for i in np.arange(0, nx-1, 1):
        for j in np.arange(0, ny, 1):
            for k in np.arange(0, nz-1, 1):
                u2[i, j, k] = 0.5 * (u[i, j, k] + u[i+1, j, k])
                v2[j, i, k] = 0.5 * (v[j, i, k] + v[j, i+1, k])
        
    for i in np.arange(0, nz-1, 1):
        for j in np.arange(0, nx, 1):
            for k in np.arange(0, ny, 1):
                w2[j, k, i] = 0.5 * (w[j, k, i] + w[j, k, i+1])


    print 'Interpolating variables to physical height grid'

    pgd   = np.zeros((nx, ny, nz2-2), float)
    ugd   = np.zeros((nx, ny, nz2-2), float)
    vgd   = np.zeros((nx, ny, nz2-2), float)
    wgd   = np.zeros((nx, ny, nz2-2), float)
    thgd  = np.zeros((nx, ny, nz2-2), float)
    rhogd = np.zeros((nx, ny, nz2-2), float)
    qvgd  = np.zeros((nx, ny, nz2-2), float)

    klo = 0
    khi = 0
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
                pgd[i, j, k-1]   = interp2hgt(zlo, zhi, p[i, j, klo], p[i, j, khi], zgrid[k])
                ugd[i, j, k-1]   = interp2hgt(zlo, zhi, u2[i, j, klo], u2[i, j, khi], zgrid[k])
                vgd[i, j, k-1]   = interp2hgt(zlo, zhi, v2[i, j, klo], v2[i, j, khi], zgrid[k])
                wgd[i, j, k-1]   = interp2hgt(zlo, zhi, w2[i, j, klo], w2[i, j, khi], zgrid[k])
                thgd[i, j, k-1]  = interp2hgt(zlo, zhi, theta[i, j, klo], theta[i, j, khi], zgrid[k])
                rhogd[i, j, k-1] = interp2hgt(zlo, zhi, rho[i, j, klo], rho[i, j, khi], zgrid[k])
                qvgd[i, j, k-1]  = interp2hgt(zlo, zhi, qvapor[i, j, klo], qvapor[i, j, khi], zgrid[k])
                

    
    GV.modelParams['NX'] = nx
    GV.modelParams['NY'] = ny
    GV.modelParams['NZ'] = nz-2
                    
    GV.modelParams['DX']  = dx
    GV.modelParams['DY']  = dy
    GV.modelParams['X0']  = -1.0 * (GV.modelParams['NX'] - 1) * GV.modelParams['DX'] / 2.0
    GV.modelParams['Y0']  = -1.0 * (GV.modelParams['NY'] - 1) * GV.modelParams['DY'] / 2.0
    GV.modelParams['T0']  = 0.0
    GV.modelParams['F0']  = f0[0,:,0]
    
    GV.modelParams['P'] = pgd
    GV.modelParams['U'] = ugd
    GV.modelParams['V'] = vgd
    GV.modelParams['W'] = wgd
    GV.modelParams['THETA'] = thgd
    GV.modelParams['RHO'] = rhogd
    GV.modelParams['QV'] = qvgd

    GV.modelParams['ZARR'] = zgrid[1:-2]

                
                
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read in the model data from plain text (ASCII) files
#
# This code will need to be changed to fit with the user's dataset and how
# that data is set up in their files.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def modelData_ReadInFromFile_PlainText(filename):
    GV.modelParams['f0']  = 2 * 7.292E-5 * math.sin(15.0 * math.pi / 180.0)
    GV.modelParams['NX']  = 251
    GV.modelParams['NY']  = 251
    GV.modelParams['NZ']  = 30
    GV.modelParams['NT']  = 240
    GV.modelParams['DX']  = 2000.0
    GV.modelParams['DY']  = 2000.0
    GV.modelParams['DT']  = 0.1
    GV.modelParams['X0']  = -1.0 * (GV.modelParams['NX'] - 1) * GV.modelParams['DX'] / 2.0
    GV.modelParams['Y0']  = -1.0 * (GV.modelParams['NY'] - 1) * GV.modelParams['DY'] / 2.0
    GV.modelParams['T0']  = 0.0

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
    GV.modelParams['ZARR'] = ZARRU
    
    NX = GV.NX
    NY = GV.NY
    NZ = GV.NZ
    DX = GV.modelParams['DX']
    DY = GV.modelParams['DY']
    
    GV.modelData['P']      = np.zeros((NX, NY, NZ), float)
    GV.modelData['U']      = np.zeros((NX, NY, NZ), float)
    GV.modelData['V']      = np.zeros((NX, NY, NZ), float)
    GV.modelData['W']      = np.zeros((NX, NY, NZ), float)
    GV.modelData['THETA']  = np.zeros((NX, NY, NZ), float)
    GV.modelData['QV']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['QL']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['DIVERG'] = np.zeros((NX, NY, NZ), float)
    GV.modelData['VORT']   = np.zeros((NX, NY, NZ), float)
    GV.modelData['AVORT']  = np.zeros((NX, NY, NZ), float)
    GV.modelData['PV']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['SPEED']  = np.zeros((NX, NY, NZ), float)
    GV.modelData['VT']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['VR']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['TEMPER'] = np.zeros((NX, NY, NZ), float)
    GV.modelData['THETAE'] = np.zeros((NX, NY, NZ), float)
    GV.modelData['TEMPVT'] = np.zeros((NX, NY, NZ), float)
    GV.modelData['RHO']    = np.zeros((NX, NY, NZ), float)
    GV.modelData['RH']     = np.zeros((NX, NY, NZ), float)
    GV.modelData['ZARR']   = np.zeros(NZ, float)

    
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
    
