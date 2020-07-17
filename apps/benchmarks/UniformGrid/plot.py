#!/usr/bin/python

import sqlite3
import sys
import re
import math
import matplotlib.pyplot as plt

from optparse import OptionParser

import ecmModel

# Defining graph types


graphs = dict(nonOpt=["cores,MAX(MLUPS)",
                      "SPLIT IS NULL AND PURE IS NULL AND CELL_OP IS NULL AND "
                      + "COMPRESSED IS NULL AND D3Q19_OPT IS NULL"],
              split=["cores,MAX(MLUPS)", "SPLIT=1"],
              pure=["cores,MAX(MLUPS)", "PURE=1"],
              normal=["cores,MAX(MLUPS)", "CELL_OP=1", ],
              compressed=["cores,MAX(MLUPS)", "COMPRESSED=1"],
              opt=["cores,MAX(MLUPS)", "D3Q19_OPT=1"],
              split_lbm=["cores,MAX(percentage)", "SPLIT=1  AND name=\"Timeloop\" AND sweep LIKE \"LBM%\""],
              split_comm=["cores,MAX(percentage)", "SPLIT=1  AND name=\"Communication\" AND sweep LIKE \"%MPI%\""],
              split_time=["cores,MAX(average)  ", "SPLIT=1  AND name=\"Timeloop\" AND sweep LIKE \"LBM%\""]
              )

# Defining all options


def addCommonOptions(parser):
    parser.add_option("", "--smIntel", action="store_true", default=False, help="On SuperMUC using IntelMPI")
    parser.add_option("", "--smIbm", action="store_true", default=False, help="On SuperMUC using IBM MPI")
    parser.add_option("", "--zyxf", action="store_true", default=False, help="Only Layout zyxf")
    parser.add_option("", "--fzyx", action="store_true", default=False, help="Only Layout fzyx")
    parser.add_option("", "--trt", action="store_true", default=False, help="Only TRT")
    parser.add_option("", "--srt", action="store_true", default=False, help="Only SRT")

    parser.add_option("", "--lto", action="store_true", default=False, help="Only where link time optimization was on")
    parser.add_option("", "--noLto", action="store_true", default=False,
                      help="Only where link time optimization was off")

    parser.add_option("", "--pgo", action="store_true", default=False,
                      help="Only where profile guided optimization was on")
    parser.add_option("", "--noPgo", action="store_true", default=False,
                      help="Only where profile guided optimization was off")

    parser.add_option("", "--intelOpt", action="store_true", default=False, help="Only where intel pragmas were used")
    parser.add_option("", "--noIntelOpt", action="store_true", default=False, help="Only where intel pragmas not used")

    parser.add_option("", "--pinCore", action="store_true", default=False, help="First fill up socket.")
    parser.add_option("", "--pinMCM", action="store_true", default=False, help="Socket round robin")

    parser.add_option("-t", "--totalMLUPS", action="store_true", default=False,
                      help="Show total MLUPS instead of MLUPS/cores")

    parser.add_option("-w", "--where", help="SQL Where or Group By clause", default="")

    parser.add_option("-p", "--printDataset", action="store_true", default=False,
                      help="Prints the dataset that is plotted")

    parser.add_option("", "--legend", default="", help="Legend entry for the graph")


def addOuterOptions(parser):
    parser.add_option("-f", "--file", default="timing.sqlite", help="Sqlite3 File with timing data")
    parser.add_option("", "--title", default="", help="Title of Graph")
    parser.add_option("-s", "--save", default="", help="Print graph to file")
    parser.add_option("", "--xlabel", default="", help="Label of x Axis")
    parser.add_option("", "--ylabel", default="", help="Label of y Axis")
    parser.add_option("", "--figureWidth", default="", help="When using --save, the width has to be specified in pts")
    parser.add_option("", "--figureHeight", default="", help="When using --save, the height has to be specified in pts")
    parser.add_option("", "--legendPos", default="", help="Position of legend (matplotlib syntax)")


def addInnerOptions(parser):
    parser.add_option("-g", "--graph", default="", help="Graph from database (measured)")
    parser.add_option("-e", "--ecm", default="", help="Plot ECM model for given kernel name")


def stringFromOptions(opt):
    c = ""
    if opt.zyxf:
        c += "FZYX "
    elif opt.fzyx:
        c += "ZYXF"
    else:
        c += "ALL"

    c += "|"

    if opt.trt:
        c += "TRT"
    elif opt.srt:
        c += "SRT"
    else:
        c += "ALL"

    c += "|"

    if opt.smIntel:
        c += "SmIntel|"
    if opt.smIbm:
        c += "SmIBM"

    if opt.where and len(opt.where) > 0:
        c += opt.where + "|"

    if opt.pinCore:
        c += "pinCore|"
    if opt.pinMCM:
        c += "pinMCM|"

    if opt.lto:
        c += "lto|"
    if opt.noLto:
        c += "noLto|"

    if opt.pgo:
        c += "pgo|"
    if opt.noPgo:
        c += "noPgo|"

    if opt.intelOpt:
        c += "intelOpt|"
    if opt.noIntelOpt:
        c += "noIntelOpt|"

    if hasattr(opt, 'graph') and opt.graph:
        c += opt.graph

    return c


def whereClauseFromOptions(opt):
    c = ""
    if opt.zyxf:
        c += "FZYX IS NULL AND "
    if opt.fzyx:
        c += "FZYX = 1 AND "

    if opt.trt:
        c += "TRT=1 AND "
    if opt.srt:
        c += "TRT IS NULL AND "

    if opt.smIntel:
        c += "buildMachine LIKE \"supermuc_intel%\" AND machine LIKE \"i__r__a__\" AND "

    if opt.smIbm:
        c += "buildMachine LIKE \"supermuc_ibm%\" AND machine LIKE \"i__r__a__\" AND "

    if opt.pinCore:
        c += "MP_TASK_AFFINITY='CORE' AND "
    if opt.pinMCM:
        c += "MP_TASK_AFFINITY='MCM' AND "

    if opt.lto:
        c += "compilerFlags LIKE '%-ipo%' AND "
    if opt.noLto:
        c += "compilerFlags NOT LIKE '%-ipo%' AND "

    if opt.pgo:
        c += "compilerFlags LIKE '%-prof-use%' AND "
    if opt.noPgo:
        c += "compilerFlags NOT LIKE '%-prof-use%' AND "

    if opt.intelOpt:
        c += "intelCompilerOpt=1 AND "
    if opt.noIntelOpt:
        c += "intelCompilerOpt IS NULL AND "

    if opt.where and len(opt.where) > 0:
        c += opt.where + " AND "

    return c

# Database and plotting


def getListFromDatabase(databaseFile, graphKeyword, whereClause):
    db = sqlite3.connect(databaseFile)
    db.row_factory = sqlite3.Row
    query = "SELECT " + graphs[graphKeyword][0]
    query += " FROM runs JOIN timingPool ON runs.runId = timingPool.runId "
    if not whereClause or len(whereClause) == 0:
        whereClause = " 1 "

    combinedWhere = " WHERE " + graphs[graphKeyword][1] + " AND " + whereClause + "GROUP BY cores"
    query += combinedWhere
    print(combinedWhere)
    c = db.cursor()
    c.execute(query)
    return c.fetchall()


def numericPlot(dataset, legendLabel="", divideByProcesses=True):
    assert (len(dataset) > 0)
    assert (isinstance(dataset[0][0], (int, float)))
    columns = len(dataset[0])
    columnNames = dataset[0].keys()

    xVals = [e[0] for e in dataset]

    plots = []
    for i in range(1, columns):
        yVals = [e[i] for e in dataset]
        if divideByProcesses:
            yVals = [yVals[i] / xVals[i] for i in range(0, len(yVals))]

        if legendLabel == "":
            legendLabel = columnNames[i]

        p, = plt.plot(xVals, yVals, marker='^', markersize=5, label=legendLabel)
        plots.append(p)

    plt.xlabel("Cores")
    plt.ylabel("MLUPS")
    # plt.gca().yaxis.grid(color='gray', linestyle='dashed')
    # plt.gca().xaxis.grid(color='gray', linestyle='dashed')


def printDataset(dataset):
    header = dataset[0].keys()
    row_format = "{:<30}" * (len(header))
    print(row_format.format(*header).upper())

    for row in dataset:
        print(row_format.format(*row))


def matplotlibLatexSetup(figureWidth, figureHeight):
    ieee = False
    if ieee:
        params = {
            'backend': 'GTKAgg',

            'font.family': 'serif',
            'font.serif': ['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman'],
            'font.sans-serif': ['Helvetica', 'Avant Garde', 'Computer Modern Sans serif'],
            # font.cursive       : Zapf Chancery
            # font.monospace     : Courier, Computer Modern Typewriter
            'text.usetex': True,

            'axes.labelsize': 9,
            'axes.linewidth': .75,

            'figure.figsize': (3.5, 2.5),
            'figure.subplot.left': 0.175,
            'figure.subplot.right': 0.95,
            'figure.subplot.bottom': 0.15,
            'figure.subplot.top': .95,

            'figure.dpi': 150,

            'text.fontsize': 9,
            'legend.fontsize': 8,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,

            'lines.linewidth': .75,
            'savefig.dpi': 600,
        }
    else:
        # Input comes usually directly from latex, and looks like "124.0pt"
        # first the "pt" is removed, then it is converted to a floating point number
        maxFigHeight = float(figureWidth.replace("pt", ""))
        maxFigWidth = float(figureHeight.replace("pt", ""))

        # Convert pt to inch
        inches_per_pt = 1.0 / 72.27
        maxFigHeight *= inches_per_pt
        maxFigWidth *= inches_per_pt

        golden_mean = (math.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio

        if maxFigWidth * golden_mean < maxFigHeight:
            height = maxFigWidth * golden_mean
            width = maxFigWidth
        else:
            height = maxFigHeight
            width = maxFigHeight / golden_mean

        fig_size = [width, height]
        params = {'backend': 'ps',
                  'font.family': 'serif',
                  'font.serif': ['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman'],
                  'font.sans-serif': ['Helvetica', 'Avant Garde', 'Computer Modern Sans serif'],
                  'axes.labelsize': 10,
                  'text.fontsize': 10,
                  'legend.fontsize': 8,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'text.usetex': True,
                  'figure.figsize': fig_size,
                  'savefig.dpi': 600,
                  'antialiased': True, }

    plt.rcParams.update(params)
    plt.rc('lines', aa=True)

# Option parsing


def bracketsAsExtraElements(list):
    """ Transforms the list [ '[abc', 'def]' ] to [ '[', 'abc', 'def', ']' ] """
    processedList = []
    for e in list:
        # TODO: write better
        splitted = re.split('(\[|\])', e)  # noqa: W605
        while '' in splitted:
            splitted.remove('')
        processedList.extend(splitted)

    return processedList


def parseList(list):
    list = bracketsAsExtraElements(list)
    outerList = []
    subLists = []

    curSubList = None
    for element in list:

        if element == "[":
            if curSubList is not None:
                raise ValueError("Two opening brackets")
            curSubList = []
        elif element == "]":
            if curSubList is None:
                raise ValueError("Closing bracket without opening bracket")
            subLists.append(curSubList)
            curSubList = None
        else:
            if curSubList is None:
                outerList.append(element)
            else:
                curSubList.append(element)

    if curSubList is not None:
        raise ValueError("Missing closing bracket")

    return outerList, subLists


outerParser = OptionParser()
addCommonOptions(outerParser)
addOuterOptions(outerParser)

innerParser = OptionParser()
addCommonOptions(innerParser)
addInnerOptions(innerParser)

(outerList, innerLists) = parseList(sys.argv)

(outerOptions, outerArgs) = outerParser.parse_args(outerList)
outerWhereClause = whereClauseFromOptions(outerOptions)

# Matplotlib setup
if outerOptions.save != "":
    matplotlibLatexSetup(outerOptions.figureWidth, outerOptions.figureHeight)

legend = []
for innerList in innerLists:
    (innerOptions, innerArgs) = innerParser.parse_args(innerList)

    if innerOptions.ecm != "":
        ecmModel.plot(innerOptions.ecm, not outerOptions.totalMLUPS, innerOptions.legend)
    else:
        innerWhereClause = whereClauseFromOptions(innerOptions)
        dataset = getListFromDatabase(outerOptions.file, innerOptions.graph,
                                      outerWhereClause + innerWhereClause + " 1 ")

        if len(dataset) == 0:
            print("Empty")
            continue

        if innerOptions.printDataset or outerOptions.printDataset:
            printDataset(dataset)

        if innerOptions.legend == "":
            innerOptions.legend = stringFromOptions(innerOptions)
        numericPlot(dataset, innerOptions.legend, not (outerOptions.totalMLUPS or innerOptions.totalMLUPS))

# Title
if outerOptions.title == "auto":
    outerOptions.title = stringFromOptions(outerOptions)

if outerOptions.title != "":
    plt.title(outerOptions.title)

# Axis labels
if outerOptions.xlabel != "":
    plt.xlabel(outerOptions.xlabel)
if outerOptions.ylabel != "":
    plt.ylabel(outerOptions.ylabel)

# Legend
if outerOptions.legendPos != "":
    plt.legend(loc=outerOptions.legendPos)
else:
    plt.legend(loc='upper left')

if outerOptions.save != "":
    plt.savefig(outerOptions.save + ".pdf")
else:
    plt.show()
