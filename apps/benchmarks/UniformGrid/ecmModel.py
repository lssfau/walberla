#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

kernels = dict()


class Kernel:

    def __init__(self, name, cyclesFirstLoop=0, cyclesSecondLoop=0, cyclesRegPerLUP=0):
        self.name = name
        if cyclesRegPerLUP <= 0:
            self.cyclesFirstLoop = cyclesFirstLoop
            self.cyclesSecondLoop = cyclesSecondLoop
            self.cyclesRegPerLUP = cyclesFirstLoop + 9 * cyclesSecondLoop
        else:
            self.cyclesRegPerLUP = cyclesRegPerLUP

        self.cyclesRegPerCacheLine = 8 * self.cyclesRegPerLUP

        self.cyclesL1L2 = 3 * 19 * 2
        self.cyclesL2L3 = 3 * 19 * 2

        self.freq = 2.7e9
        self.cyclesMem = 305
        # self.cyclesMem = 191

    def mlups(self, processes):
        singleCoreCycles = self.cyclesRegPerCacheLine + self.cyclesL1L2 + self.cyclesL2L3 + self.cyclesMem

        timeSingleCore = singleCoreCycles / self.freq

        mlups = 8 / timeSingleCore * 1e-6

        # todo
        mlupsMax = 78

        return min(processes * mlups, mlupsMax)

    def plot(self, divideByProcesses=False, processes=8, label=""):

        x = np.arange(1, processes + 1, 1)
        if divideByProcesses:
            y = np.array([self.mlups(i) / i for i in x])
        else:
            y = np.array([self.mlups(i) for i in x])

        if label == "":
            label = "ecm_" + self.name
        plt.plot(x, y, marker='^', markersize=5, label=label)


kernels = dict()

# kernels['srt_split'] = Kernel("srt_split", 46, 12 )

kernels['srt_pure'] = Kernel("srt_pure", 40, 8)
kernels['trt_split'] = Kernel("trt_split", 41, 11)
# SRTStreamCollide.h -  pgo and lto (20cycles first loop, 35 second)
kernels['srt_nonopt'] = Kernel("srt_nonopt",
                               cyclesRegPerLUP=1045)


# kernels['trt_pure_intelOpt'] = Kernel("trt_pure_intelOpt", 41/2, 10/2 )  # vectorized (v*pd)


def plotAllKernels(divideByProcesses=False):
    for kernel in kernels:
        kernel.plot(divideByProcesses)


def plot(kernelName, divideByProcesses=False, label=""):
    kernels[kernelName].plot(divideByProcesses, label=label)


if __name__ == "__main__":
    plotAllKernels()
    plt.legend()
    plt.show()
