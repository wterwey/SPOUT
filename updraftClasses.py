class updraft2D:
    def __init__(self, x, y, z, t, data):
        self.xPos = x
        self.yPos = y
        self.zPos = z
        self.time = t
        self.savedData = data

class updraft3D:
    def __init__(self, data, startHeight):
        try:
            (x for x in data)
            self.linked2DUpdrafts = data
        except TypeError:
            self.linked2DUpdrafts = [data]
        self.bottomHeight = startHeight
        self.topHeight    = startHeight

#    def __init__(self):
#        self.linked2DUpdrafts = []

    def vertRange(self):
        values = [self.bottomHeight, self.topHeight]
        return values

class updraft4D:
    def __init__(self, data, startTime, identNum):
        try:
            (x for x in data)
            self.linked3DUpdrafts = data
        except TypeError:
            self.linked3DUpdrafts = [data]
        self.beginTime = startTime
        self.endTime   = startTime
        self.ident     = identNum

    def findTime(self, tVal):
        if (tVal <= self.endTime and tVal >= self.beginTime):
            for temp3D in self.linked3DUpdrafts:
                if (temp3D.linked2DUpdrafts[0].time == tVal):
                    return temp3D
        else:
            return -1


#    def __init__(self):
#        self.linked3DUpdrafts = []



