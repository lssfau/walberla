from __future__ import print_function
import os
import hashlib
import json
import time

try:
    from ..walberla_cpp import mpi
except ImportError:
    from walberla_cpp import mpi


def waitForFoldersToExist(folders, waitTime=3, timeout=10):
    """For big MPI parallel runs it can happen that mkdir returns but not all processes already
       see the newly created folders due to network filesystem issues. This function waits until
       the given folders exist."""

    def allFoldersExits():
        for f in folders:
            if not os.path.exists(f):
                return False
        return True

    waitedTime = 0
    while not allFoldersExits():
        if waitedTime > timeout:
            raise RuntimeError(
                "Creating folder timeout:  when waiting for folder creation. Network filesystem problems?" + str(
                    folders))
        time.sleep(waitTime)
        waitedTime += waitTime


class HashedScenarioFolderCreator:
    """
    Helper class to create the following directory structure for simulation results

    Example:

        baseFolder / 2fc7c861eb/ subfolderNames[0]
                               / subfolderNames[1]
                               / subfolderNames[2]
                               params.json
                  / 3413d21776/ subfolderNames[0]
                              / subfolderNames[1]
                              / subfolderNames[2]
                              params.json

    In the baseFolder several scenario folder are created. The name of this folders
    is  a hash value of the parameters of the scenarios.
    The subfolder names (in this example meshes, logs, img ) can be configured and
    are passed as folderNames in the constructor. Optionally a json file is written
    to the scenario folders that has contains the parameters in readable form.

    If a scenario was already simulated, i.e. the folder does already exists, a new folder is
    created with "_1" appended to the name.

    All filesystem operations are carried out by root process, but function have to be called
    by all processes.
    """

    def __init__(self, baseFolder, subfolderNames, hashFunc=hashlib.sha1, hashLength=10):
        self._baseFolder = baseFolder
        self._subfolderNames = subfolderNames
        self._hashLength = hashLength
        self._hashFunc = hashFunc
        # create base folder if it does not exist yet
        if mpi.worldRank() == 0 and not os.path.exists(self._baseFolder):
            os.makedirs(self._baseFolder)

    def _scenarioBaseName(self, paramsDict, postfix=0):
        """Scenario name without numbering postfix"""
        # Hash function needs encoded string -> encode the dictionary using json
        encodedJSON = json.dumps(paramsDict, sort_keys=True).encode('utf-8')
        hashStr = self._hashFunc(encodedJSON).hexdigest()[:self._hashLength]
        if postfix > 0:
            return "%s_%d" % (hashStr, postfix)
        else:
            return hashStr

    def scenarioExists(self, paramsDict, postfix=0):
        postfixStr = "_%d" % (postfix,) if postfix > 0 else ""
        path = os.path.join(self._baseFolder, self._scenarioBaseName(paramsDict) + postfixStr)
        return os.path.exists(path)

    def create(self, paramsDict, writeParamsJSON=True):
        """Creates a new scenario based on the given parameters, waits until all
           subfolders are created on all processes.
           Returns path to scenario folder"""
        if mpi.worldRank() == 0:
            postfix = 0
            while self.scenarioExists(paramsDict, postfix):
                postfix += 1
            scenarioName = self._scenarioBaseName(paramsDict, postfix)
            scenarioName = mpi.broadcastString(scenarioName)
        if mpi.worldRank() != 0:
            scenarioName = mpi.broadcastString("")

        folders = [os.path.join(self._baseFolder, scenarioName, subfolder) for subfolder in [""] + self._subfolderNames]
        scenarioFolder = folders[0]

        if mpi.worldRank() == 0:
            for f in folders:
                os.mkdir(f)
            if writeParamsJSON:
                paramsDict['scenarioName'] = scenarioName
                with open(os.path.join(scenarioFolder, "params.json"), 'w') as f:
                    json.dump(paramsDict, f, indent=4, sort_keys=True)

        mpi.worldBarrier()
        waitForFoldersToExist(folders)
        mpi.worldBarrier()

        return scenarioName, scenarioFolder
