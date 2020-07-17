import re
import os

getAppDirRegexp = re.compile(r'(src|test|apps)')


def getAppDir(filename):
    """Returns the first part of the path without (src|tests|app)/someModule/someFile.cpp
      Example: filename = ~/devel/walberla/src/core/timing/Timer.cpp
               returns    ~/devel/walberla           """
    filename = os.path.abspath(filename)
    matchObj = getAppDirRegexp.search(filename)
    return filename[:matchObj.start() - 1]


def getAppModules(appDir):
    """Returns all subfolders of appDir/src  i.e. all modules given the appDir"""
    searchDir = appDir + "/src"
    return [d for d in os.listdir(searchDir) if os.path.isdir(os.path.join(searchDir, d))]


getModuleRegexp = re.compile(r'(src|tests)/([^/]*)')


def getModule(filename):
    """Returns the  module of a file
       Examples:   somepath/walberla/src/core/something.cpp       returns core
                   somePrefix/yourApp/tests/greatModule/file.cpp  returns greatModule"""
    r = getModuleRegexp.search(filename)
    if r and '.' not in r.group(2):
        return r.group(2)
    else:
        return ""


def isDeprecated(filename):
    file = open(filename)
    content = file.read()
    return "#pragma message" in content and "deprecated" in content
