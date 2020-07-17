try:
    from . import walberla_cpp
except ImportError:
    import walberla_cpp


def functorFromSweep(blocks, sweep):
    """Makes a functor from a sweep, by iterating over all the blocks in the provided block storage"""

    def functor():
        for b in blocks:
            sweep(b)

    return functor


class Timeloop(walberla_cpp.timeloop.ITimeloop):
    def __init__(self, nrOfTimesteps):
        super().__init__()
        self._nrOfTimesteps = nrOfTimesteps
        self._timestep = 0
        self._functors = []
        self._stopFlag = False

    def run(self, timesteps=None):
        if not timesteps:
            timesteps = self._nrOfTimesteps
        for t in range(timesteps):
            self.singleStep()
            if self._stopFlag:
                break

    def singleStep(self):
        for func in self._functors:
            func()
        self._timestep += 1

    def stop(self):
        self._stopFlag = True

    def synchronizedStop(self, stop=True):
        # syncStop = wlb.mpi.allreduceInt(int(stop), wlb.mpi.LOGICAL_OR)
        self._stopFlag = True

    def setCurrentTimeStep(self, ts):
        self._timestep = ts

    def getCurrentTimeStep(self):
        return self._timestep

    def getNrOfTimeSteps(self):
        return self._nrOfTimesteps

    def __getFunctor(self, functor, blocks):
        if blocks:  # assume that it is a sweep if blocks were given
            return functorFromSweep(blocks, functor)
        else:
            return functor

    def add(self, functor, blocks=None):
        self._functors.append(self.__getFunctor(functor, blocks))
        return len(self._functors) - 1

    def replace(self, handle, functor, blocks=None):
        self._functors[handle] = self.__getFunctor(functor, blocks)


def extend(cppTimeloopModule):
    cppTimeloopModule.Timeloop = Timeloop
