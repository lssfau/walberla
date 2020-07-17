#!/usr/bin/python
import re
import numpy as np

# -------------------------------------------------------------------------------------------------
# --------------------    Definition of all Stencils    -------------------------------------------
# -------------------------------------------------------------------------------------------------


# Template input file
templateFile = "Stencil.in.h"

# Directions have to be in the same order as they are defined in the Directions.h
# they also have to have same names

directions = ['C', 'N', 'S', 'W', 'E', 'T', 'B',
              'NW', 'NE', 'SW', 'SE', 'TN', 'TS',
              'TW', 'TE', 'BN', 'BS', 'BW', 'BE',
              'TNE', 'TNW', 'TSE', 'TSW', 'BNE', 'BNW', 'BSE', 'BSW']

directionCoords = [[0, 0, 0],  # C
                   [0, 1, 0],  # N
                   [0, -1, 0],  # S
                   [-1, 0, 0],  # W
                   [1, 0, 0],  # E
                   [0, 0, 1],  # T
                   [0, 0, -1],  # B
                   ]

# List of all stencils
# Edit this to add new stencils (make sure name is unique)
stencils = [
    {'name': 'D2Q4', 'dim': 2, 'dirs': directions[1: 5]},
    {'name': 'D2Q5', 'dim': 2, 'dirs': directions[: 5]},
    {'name': 'D2Q9', 'dim': 2, 'dirs': ['C', 'N', 'S', 'W', 'E', 'NW', 'NE', 'SW', 'SE']},
    {'name': 'D2CornerStencil', 'dim': 2, 'dirs': ['NW', 'NE', 'SW', 'SE']},
    {'name': 'D3Q6', 'dim': 3, 'dirs': directions[1: 7]},
    {'name': 'D3Q7', 'dim': 3, 'dirs': directions[: 7]},
    {'name': 'D3Q19', 'dim': 3, 'dirs': directions[:19]},
    {'name': 'D3Q15', 'dim': 3,
     'dirs': ['C', 'N', 'S', 'W', 'E', 'T', 'B', 'TNE', 'TNW', 'TSE', 'TSW', 'BNE', 'BNW', 'BSE', 'BSW']},
    {'name': 'D3Q27', 'dim': 3, 'dirs': directions[:27]},
    # also possible: pick directions as you need them in arbitrary order:
    {'name': 'EdgeStencil', 'dim': 3, 'dirs': ['NW', 'NE', 'SW', 'SE', 'TN', 'TS', 'TW', 'TE', 'BN', 'BS', 'BW', 'BE']},
    {'name': 'D3CornerStencil', 'dim': 3, 'dirs': directions[19:27]},
    {'name': 'D3EdgeCornerStencil', 'dim': 3, 'dirs': directions[7:27]}
]

# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# --------------------    Code for Stencil Generation -----------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

def directionToCoordinate(dir):
    newCoord = np.zeros(3, dtype='i')
    for character in dir:
        index = directions.index(character)
        newCoord = newCoord + np.array(directionCoords[index], dtype='i')
    return list(newCoord)


def coordinateToDirection(coord):
    xPart = [coord[0], 0, 0]
    yPart = [0, coord[1], 0]
    zPart = [0, 0, coord[2]]

    directionStr = ""
    if zPart in directionCoords[1:]:
        directionStr += directions[directionCoords.index(zPart)]
    if yPart in directionCoords[1:]:
        directionStr += directions[directionCoords.index(yPart)]
    if xPart in directionCoords[1:]:
        directionStr += directions[directionCoords.index(xPart)]

    if directionStr == "":
        directionStr = "C"

    return directionStr


header = """//====================================================================================================================
//  Caution: This file has been generated automatically. All manual changes are lost when file is regenerated!
//           Changes should be done in Stencil.in.h,and then all stencils classes can be generated again.
//====================================================================================================================
#ifndef DOXY_SKIP_INTERNAL
"""

footer = """
#endif // DOXY_SKIP_INTERNAL
"""


def isSubDirection(general, specific):
    """ Example: general="N", specific= "NW" -> true
                 general="W", specific= "SE" -> false
    """
    for char in general:
        if char not in specific:
            return False
    return True


def indexFromDir(dirs):
    """ Return list of indices that elements in dir array have in global
        directions array"""
    res = []
    index = 0
    for d in directions:
        if (d in dirs):
            res.append(index)
            index = index + 1
        else:
            res.append("INVALID_DIR")

    return res


def generate_d_per_d(dirs):
    """ Generate d_per_d array from directions"""
    d_per_d = []
    d_per_d_length = []

    for globalDir1 in directions:
        subdirs = []
        for localDir in dirs:
            if isSubDirection(globalDir1, localDir):
                subdirs.append(localDir)

        d_per_d.append("{" + ",".join(subdirs) + "}")
        d_per_d_length.append(len(subdirs))

    return (d_per_d, d_per_d_length)


def generate_dir_pos(dirs):
    """ Generates an array containing only half of the directions. This can be used to iterate over
        half the stencil, fetching the remaining direction using inv_dir.
        Only works for symmetrical stencils"""
    result = []
    for d in dirs:
        c = directionToCoordinate(d)
        if c[0] == 1 or (c[0] == 0 and c[1] == 1) or (c[0] == 0 and c[1] == 0 and c[2] == 1):
            result.append(d)

    # assert( len(dirs) / 2 == len(result) ) # Fails for asymmetric stencils
    return result


def getNeighborsOfDirection(startDir, walkingDirs, possibleResults):
    """ This function views a direction as a cell (i.e. the cell where the direction points to)
        The whole domain consists of the 3x3x3 grid (other cells do not exist).
        The function returns all neighboring cells ( in the 3x3x3 grid ) of a given cell

        Input: a start direction (i.e. cell where to start)
               walkingDirs: list of directions, in which the neighboring directions are searched, starting from startDir
                            in this argument, directions are interpreted as directions not as cells!
               possibleResults: if direction not in this list, then it is not included in the results
               """
    result = []
    curDirCoord = directionToCoordinate(startDir)

    for walkingDir in walkingDirs:
        neighborCoord = directionToCoordinate(walkingDir)
        sumDirCoord = [curDirCoord[0] + neighborCoord[0],
                       curDirCoord[1] + neighborCoord[1],
                       curDirCoord[2] + neighborCoord[2]]

        if sumDirCoord[0] > 1 or sumDirCoord[0] < -1:
            continue
        if sumDirCoord[1] > 1 or sumDirCoord[1] < -1:
            continue
        if sumDirCoord[2] > 1 or sumDirCoord[2] < -1:
            continue

        sumDir = coordinateToDirection(sumDirCoord)

        if sumDir not in possibleResults:
            continue

        if sumDir == startDir:
            continue

        result.append(sumDir)

    result = sorted(result, key=lambda d: directions.index(d))

    return result


def generateNeighborsOfDirection(stencilDirs):
    neighborList = []
    lengthList = []
    for d in directions:
        neighbors = getNeighborsOfDirection(d, stencilDirs, directions)
        neighborList.append("{" + ",".join(neighbors) + "}")
        lengthList.append(str(len(neighbors)))

    return (neighborList, lengthList)


# ---------------------------------------------------------------------------------------------------------------------


tmplFile = open(templateFile, "r")

# Pythons string.format() function takes a string where
# the parts that should be substituted# are marked with "{varname}"
# Every curly bracket that should be left unaffected by the substitution has  be expressed as double curly bracket
# Since we want to do the substitution in a C++ Header file, this syntax is not very comfortable.

# So we define a new substitution rule: instead of writing {varname} we write $varname
# The first step is not, to prepare the read template file

# Because of Python's substitution mechanism, all curly brackets have to be replaced by double
# curly brackets
tmpl = tmplFile.read().replace("{", "{{")
tmpl = tmpl.replace("}", "}}")

# Then we replace all "$var" constructs by "{var}"
tmpl = re.sub("\$[A-Za-z_]*", (lambda x: "{" + x.group(0)[1:] + "}"), tmpl)  # noqa: W605

"""Reads stencil dict and generates header files"""
for stencil in stencils:
    # Make sure that all directions are sorted as in global direction array
    # and implicitly check that all entries of dirs are in directions array
    dirs = sorted(stencil['dirs'], key=lambda d: directions.index(d))
    name = stencil['name']

    # Build up a dictionary of replacements
    vals = dict()
    vals['name'] = name
    vals['capitalName'] = name.upper()
    vals['dirs'] = ",".join(dirs)
    vals['indexFromDir'] = ",".join(str(i) for i in indexFromDir(dirs))
    vals['D'] = stencil['dim']
    vals['Q'] = len(dirs)

    (d_per_d, d_per_d_length) = generate_d_per_d(dirs)
    vals['d_per_d'] = ",\n\t\t\t\t\t\t\t\t".join(d_per_d)
    vals['d_per_d_length'] = ",".join(str(i) for i in d_per_d_length)
    vals['containsCenter'] = "true" if ('C' in dirs) else "false"
    vals['noCenterFirstIndex'] = "1" if ('C' in dirs) else '0'

    (dir_neighbors, dir_neighbors_length) = generateNeighborsOfDirection(dirs)
    vals['dir_neighbors'] = ",\n\t\t\t\t\t\t\t\t".join(dir_neighbors)
    vals['dir_neighbors_length'] = ",".join(str(i) for i in dir_neighbors_length)

    vals['dir_pos'] = ",".join(generate_dir_pos(dirs))

    content = tmpl.format(**vals)

    # Open file for writing
    out = open(name + ".h", 'w')
    out.write(header)
    out.write(content)
    out.write(footer)
    out.close()
