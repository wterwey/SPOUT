import math
import numpy as np
global modelParams, TrackerParams, modelData

#----------------------------------------------------------------------
# Set up common arrays and constants.
#----------------------------------------------------------------------
modelParams = {}
modelParams['INPUTROOT']  = 'C:\Users\wterwey\Desktop\WRFFiles'
modelParams['OUTPUTROOT'] = 'C:\Users\wterwey\Desktop\WRFTest'
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
# velCorrect is the correction factor for the advective process in ranking
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
TrackerParams['velCorrect']        = 1.0

tempar = np.zeros(TrackerParams['NUMLEVELS']) + TrackerParams['MINUPLEVEL']
TrackerParams['checkLevels'] = tempar

tempar = np.zeros(TrackerParams['NUMLEVELS'], float) + 1.5
TrackerParams['checkValues'] = tempar

updraftData = []


#######################################################################
#----------------------------------------------------------------------
# Set up data array dictionary
#
# This dictionary will be populated by the modelReadIn modules.
#----------------------------------------------------------------------

modelData = {}
