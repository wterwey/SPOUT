import math
import numpy as np
global modelParams, TrackerParams, modelData

#----------------------------------------------------------------------
# Set up common arrays and constants.
#----------------------------------------------------------------------
modelParams = {}

modelParams['f0']  = 2 * 7.292E-5 * math.sin(15.0 * math.pi / 180.0)
modelParams['NX']  = 251
modelParams['NY']  = 251
modelParams['NZ']  = 30
modelParams['NT']  = 240
modelParams['DX']  = 2000.0
modelParams['DY']  = 2000.0
modelParams['DT']  = 0.1
modelParams['X0']  = -1.0 * (modelParams['NX'] - 1) * modelParams['DX'] / 2.0
modelParams['Y0']  = -1.0 * (modelParams['NY'] - 1) * modelParams['DY'] / 2.0
modelParams['T0']  = 0.0

modelParams['ZARR'] = []
modelParams['FILEROOT'] = './AnalysisFiles/ascii'
modelParams['g'] = 9.8

#----------------------------------------------------------------------
# Establish updraft arrays and parameters.
#
# BADVALUE is a tracking value that can be used for checking.  Choose
# a value that NONE of your potential tracked variables could EVER
# take on.
#
# NVAR is the number of variables being tracked.  This needs to be a
# number greater than or equal to 5.  The first five variables tracked
# are X-position (0), Y-position (1), U-environment (2), V-environment
# (3), and W-insitu (4).  Any and all other variables can be modified
# by the user in the appropriate section of code.
#
# MAXUP is the maximum number of updrafts that can be tracked at
# once.
#
# WTHRES is a threshold value for checking by how much a local maximum
# must exceed its neighbors.
#
# MINUPHEIGHT is a value stating the minimum updraft height to be
# considered (in vertical gridpoints).
#
# NUMLEVELS states through how many vertical gridpoints the code will
# be searching for updrafts.  It is advisable to count the number of
# gridpoints between the boundary layer and the tropopause and use a
# number close to this.
#
# SPACETHRES is the spatial threshold (in horizontal gridpoints) that
# must be satisfied for the vertical linking of updrafts.  This is
# used for both horizontal dimensions independently, not as a total
# distance.
#
# RANKTHRES is the spatial threshold (in horizontal gridpoint
# distance) that must be satisfied for the temporal linking of
# updrafts.  Using the previous time step's environmental U and
# V winds (and the velFudge factor), a prior updraft is "advected" to
# a new point.  For it to be (potentially) temporally linked, it must
# be within RANKTHRES distance of this advected point.
# 
# POINTTHRES is the vertial spatial threshold that must be satisfied
# for the temporal linking of updrafts.  Prior and current updrafts
# must have a number of updraft levels within POINTTHRES number of
# vertical points to be temporally linked.
#
# OUTPUTTIMETHRES is a threshold value for the minimum number of time
# steps that a tracked updraft must be to be output to the final data
# file.
#
# velFudge is the fudge-factor for the advective process in ranking
# potential temporal linking.  This was included since it is
# well-known that thunderstorms can often move at fractions of the
# imposed advective velocities.
#
# checkLevels is the array that states which set of consecutive levels
# (in vertical gridpoint space) are being checked.
#
# checkValues is the array that states what the minimal vertical
# velocity is for each level's maxima to be qualified for
# updraft status.
# 
# updraftData is the actual array that tracks the information for the
# updrafts currently being tracked.  It is a 4-D array:
#   1.) ID of the variable being tracked
#   2.) ID number of the updraft
#   3.) Vertical profile of each updraft
#   4.) Time for each updraft.
# This array is filled usually with BADVALUE. 
#----------------------------------------------------------------------

TrackerParams = {}
TrackerParams['BADVALUE']        = -999.9
TrackerParams['WTHRES']          = 0.0
TrackerParams['MINUPLEVEL']      = 3
TrackerParams['MINUPHEIGHT']     = 4
TrackerParams['NUMLEVELS']       = 20
TrackerParams['SPACETHRES']      = 2
TrackerParams['RANKTHRES']       = 5.0
TrackerParams['POINTTHRES']      = 6
TrackerParams['OUTPUTTIMETHRES'] = 3
TrackerParams['velFudge']        = 1.0

tempar = np.array(range(TrackerParams['NUMLEVELS']))
tempar = tempar + TrackerParams['MINUPLEVEL']
TrackerParams['checkLevels'] = tempar

tempar = np.array(range(TrackerParams['NUMLEVELS']), float)
tempar = 0.2 * tempar + 0.6
TrackerParams['checkValues'] = tempar

updraftData = []


#######################################################################
# File read-in...
# 
# This portion of code will need to be changed if one is using a
# different model result.  This portion of code can be expanded to
# include any and all read-in, but the important data arrays (U, V, W,
# etc.) need to be filled appropriately.

#----------------------------------------------------------------------
# Set up data arrays
#----------------------------------------------------------------------
NX = modelParams['NX']
NY = modelParams['NY']
NZ = modelParams['NZ']

P      = np.zeros((NX, NY, NZ), float)
U      = np.zeros((NX, NY, NZ), float)
V      = np.zeros((NX, NY, NZ), float)
W      = np.zeros((NX, NY, NZ), float)
THETA  = np.zeros((NX, NY, NZ), float)
QV     = np.zeros((NX, NY, NZ), float)
QL     = np.zeros((NX, NY, NZ), float)
DIVERG = np.zeros((NX, NY, NZ), float)
VORT   = np.zeros((NX, NY, NZ), float)
AVORT  = np.zeros((NX, NY, NZ), float)
PV     = np.zeros((NX, NY, NZ), float)
SPEED  = np.zeros((NX, NY, NZ), float)
VT     = np.zeros((NX, NY, NZ), float)
VR     = np.zeros((NX, NY, NZ), float)
TEMPER = np.zeros((NX, NY, NZ), float)
THETAE = np.zeros((NX, NY, NZ), float)
TEMPVT = np.zeros((NX, NY, NZ), float)
RHO    = np.zeros((NX, NY, NZ), float)
RH     = np.zeros((NX, NY, NZ), float)
ZARR   = np.zeros(NZ, float)

modelData = {'P' : P, 'U' : U, 'V' : V, 'W' : W, 'THETA' : THETA, 
             'QV' : QV, 'QL' : QL, 'DIVERG' : DIVERG, 'VORT' : VORT,
             'AVORT' : AVORT, 'PV' : PV, 'SPEED' : SPEED, 'TEMPER' : TEMPER, 
             'THETAE' : THETAE, 'TEMPVT' : TEMPVT, 'RHO' : RHO, 'RH' : RH,
             'ZARR' : ZARR}
