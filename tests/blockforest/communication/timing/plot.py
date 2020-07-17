import numpy as np
import matplotlib
import copy
import matplotlib.pyplot as plt

matplotlib.use('Qt4Agg')


# ------------------- Timing -------------------------------------

class Timing:
    def __init__(self):
        self.min = 0
        self.max = 0
        self.avg = 0

    def get(self, value):
        if value == "min":
            return self.min
        if value == "max":
            return self.max
        if value == "avg":
            return self.avg

    def __add__(self, other):
        ret = copy.deepcopy(self)
        ret.min += other.min
        ret.max += other.max
        ret.avg += other.avg
        return ret

    def __div__(self, other):
        ret = copy.deepcopy(self)
        ret.min /= other
        ret.max /= other
        ret.avg /= other
        return ret

    def readFromArray(self, arr):
        self.min = np.double(arr[0])
        self.avg = np.double(arr[1])
        self.max = np.double(arr[2])


# ---------------Module Timing------------------------------------


class ModuleTiming:
    def __init__(self):
        self.total = Timing()
        self.pack = Timing()
        self.mpi = Timing()
        self.unpack = Timing()

    def __add__(self, other):
        ret = copy.deepcopy(self)
        ret.total += other.total
        ret.pack += other.pack
        ret.mpi += other.mpi
        ret.unpack += other.unpack
        return ret

    def __div__(self, other):
        ret = copy.deepcopy(self)
        ret.total /= other
        ret.pack /= other
        ret.mpi /= other
        ret.unpack /= other
        return ret

    def readFromArray(self, arr):
        self.pack = Timing()
        self.pack.readFromArray(arr[0:3])

        self.mpi = Timing()
        self.mpi.readFromArray(arr[3:6])

        self.unpack = Timing()
        self.unpack.readFromArray(arr[6:9])


# ---------------Timing DataSet ----------------------------------


class TimingDataSet:
    """ Represents one line of a timing file"""

    def __init__(self):
        self.functionality = ""
        self.timesteps = 0
        self.cores = 0
        self.blocksize = np.array([0, 0, 0])
        self.blocks = np.array([0, 0, 0])

        self.oldModule = ModuleTiming()
        self.newModule = ModuleTiming()

    def __add__(self, other):
        assert (self == other)
        ret = copy.deepcopy(self)
        ret.oldModule += other.oldModule
        ret.newModule += other.newModule
        return ret

    def __div__(self, other):
        ret = copy.deepcopy(self)
        ret.oldModule /= other
        ret.newModule /= other
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        if not isinstance(other, TimingDataSet):
            return False

        return (self.functionality == other.functionality) and \
               (self.timesteps == other.timesteps) and \
               (self.cores == other.cores) and \
               (self.blocksize == other.blocksize).all and \
               (self.blocks == other.blocks).all

    def normalizeUnitTimeStep(self):
        self.oldModule /= self.timesteps
        self.newModule /= self.timesteps
        self.timesteps = 1

    def readFromLine(self, line):
        splitted = line.split()
        self.functionality = splitted[0]
        self.timesteps = np.double(splitted[1])
        self.cores = int(splitted[2])
        self.blocks = np.array([splitted[3], splitted[4], splitted[5]])
        self.blocksize = np.array([splitted[6], splitted[7], splitted[8]])

        self.oldModule = ModuleTiming()
        self.oldModule.readFromArray(splitted[9:21])

        self.newModule = ModuleTiming()
        self.newModule.readFromArray(splitted[21:33])

    def getLabel(self):
        ret = self.functionality + "\n"
        ret += "blocks:  " + str(self.blocks) + "\n"
        ret += "bl size: " + str(self.blocksize) + "\n"
        ret += str(self.cores) + " cores"
        return ret


def parseTimingFile(fileObject):
    """Parses Timing file of following format
       functionality #timesteps #cores  #blocks #blocksize
       oldPack oldMpi oldUnpack oldTotal newPack newMpi newUnpack newTotal
       """

    timingDataSetList = []
    for line in fileObject:
        dataset = TimingDataSet()
        dataset.readFromLine(line)
        timingDataSetList.append(dataset)

    return timingDataSetList


def plotListOfDataSets(l, value="avg"):
    labels = [i.getLabel() for i in l]
    labelLoc = np.array([i * 3 for i in range(0, len(l))])

    xLocations = np.empty([len(l) * 2])
    packing = np.empty([len(l) * 2])
    mpi = np.empty([len(l) * 2])
    unpacking = np.empty([len(l) * 2])
    rest = np.empty([len(l) * 2])

    for i in range(0, len(l)):
        xLocations[2 * i] = i * 3
        xLocations[2 * i + 1] = i * 3 + 1

        packing[2 * i] = l[i].oldModule.pack.get(value)
        packing[2 * i + 1] = l[i].newModule.pack.get(value)

        mpi[2 * i] = l[i].oldModule.mpi.get(value)
        mpi[2 * i + 1] = l[i].newModule.mpi.get(value)

        unpacking[2 * i] = l[i].oldModule.unpack.get(value)
        unpacking[2 * i + 1] = l[i].newModule.unpack.get(value)

        rest[2 * i] = l[i].oldModule.total.get(value) - packing[2 * i] - unpacking[2 * i] - mpi[2 * i]
        rest[2 * i + 1] = l[i].newModule.total.get(value) - packing[2 * i + 1] - unpacking[2 * i + 1] - mpi[2 * i + 1]

    plt.bar(xLocations, packing, width=1.0, color='#00e9e1')
    plt.bar(xLocations, mpi, width=1.0, bottom=packing, color='#5af441')
    plt.bar(xLocations, unpacking, width=1.0, bottom=packing + mpi, color='#4d4dff')
    plt.xticks(labelLoc + 1, labels)
    # plt.bar(xLocations,rest,     width=1.0, bottom=packing+mpi+unpacking,color='y')
    plt.show()


def combineDataSetsWithSameInput(l):
    l = copy.deepcopy(l)
    res = []

    while (len(l) > 0):
        curDataSet = l[0]
        sameDataList = [x for x in l if x == curDataSet]
        l = [x for x in l if x != curDataSet]

        # calculate average over all timings with same simulation parameters
        sum = sameDataList[0]
        for i in range(1, len(sameDataList)):
            sum += sameDataList[i]

        res.append(sum / len(sameDataList))

    return res


# timingFile = "/home/bauer/code/walberlaGit/bin/tests/modules/communication2/timing.out"
timingFile = "/home/bauer/code/walberlaGit/timing.out"

l = parseTimingFile(open(timingFile))
l = combineDataSetsWithSameInput(l)
plotListOfDataSets(l)
