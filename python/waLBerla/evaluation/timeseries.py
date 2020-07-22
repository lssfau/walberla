import numpy as np
import time


class TimeSeries:
    """Stores timeseries of values that are occuring during a simulation.
       Can estimate runtime until a value has reached a certain threshold.
       Data is collected in a format that can be written to a sqlite database easily."""

    def __init__(self):
        self._data = dict()
        self._timestamps = []

    def __getitem__(self, dataName):
        return self._data[dataName]

    @property
    def dataDict(self):
        return self._data

    def addDatapoint(self, datapoint):
        """Add several measure quantities.
            Example:
               addDatapoint( { 'maxVelocity' : 0.04, 'dropWidth' : 8, } ) """
        if len(self._data) == 0:  # on first call
            for key in datapoint:
                self._data[key] = []

        if set(datapoint.keys()) != set(self._data.keys()):
            raise ValueError("Datapoints have to contain always the same quantities")

        self._timestamps.append(time.time())
        for key, val in datapoint.items():
            self._data[key].append(val)

    def __len__(self):
        return len(self._timestamps)

    def estimateRemainingTimeLinear(self, dataName, targetValue, historyLength=0):
        """
        :param historyLength: a linear estimate is generated using the last added value and the
                              value data[historyLength]
                              for out of range values abs(historyLength) > len(data): the first data point is chosen
        :return: estimated number of seconds until targetValue is reached
                (may be negative if distance to targetValue is increasing)
        """
        if abs(historyLength) > len(self):
            historyLength = 0
        if len(self) < 2:
            return None

        data = self._data[dataName]

        dt = self._timestamps[-1] - self._timestamps[historyLength]
        dVal = data[-1] - data[historyLength]

        distance = targetValue - data[-1]
        if dVal == 0:
            return None

        res = distance / dVal * dt
        return None if res < 0 else res

    def isStationary(self, dataName, historyLength, tolerance):
        """Views only the last historyLength values. Returns True if none of these values
           does deviate more than tolerance from their average"""
        if len(self) < historyLength:
            return False

        data = np.array(self._data[dataName][-historyLength:])
        data = np.abs(data - np.average(data))
        return np.max(data) < tolerance
