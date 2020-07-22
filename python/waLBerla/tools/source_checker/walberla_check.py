#!/usr/bin/env python


import os
import re
from Utils import getAppDir, getAppModules, getModule, isDeprecated
from ParsedCodeFile import OldStyleParsedCodeFile, ParsedCodeFile, ParsingException

# Settings


skipDirectories = ['extern', 'stencil', 'vof']

walberlaDir = os.environ['WALBERLA_SOURCE_DIR']

authorfile = walberlaDir + "/utilities/py_waLBerla/source_checker/authors.txt"
licensefile = walberlaDir + "/utilities/py_waLBerla/source_checker/license.txt"


class CodeStyleException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


# Separator Line Length

def getSortedIncludes(parsedCodeFile):
    filename = parsedCodeFile.filename
    includeList = parsedCodeFile.includes

    appDir = getAppDir(filename)
    appName = [x for x in appDir.split('/') if len(x) > 0][-1]
    appName = appDir.split('/')[-1]
    if appName == "walberla":
        appModules = []
    else:
        appModules = getAppModules(appDir)

    filename = filename  # name of file containing the includes
    currentModule = getModule(filename)

    localIncludes = []
    moduleIncludes = []
    appIncludes = []
    externIncludes = []

    for include in includeList:
        if '"' in include and "/" not in include:
            localIncludes.append(include)
        elif '"' in include:
            moduleOfInclude = include.split('/')[0][1:]
            if moduleOfInclude in appModules:
                appIncludes.append(include)
            elif moduleOfInclude == currentModule:
                localIncludes.append(include)
            else:
                moduleIncludes.append(include)
        elif '<' in include:
            externIncludes.append(include)
        elif len(include.strip()) == 0:
            continue
        else:
            raise Exception("Internal Error when parsing include: " + include)

    def keyFunction(include):
        """Key function for sorting module and local includes"""
        splitted = include.split('/')
        return splitted[0] + "/" + str(len(splitted)) + "/" + "/".join(splitted[1:])

    def externKeyFunction(include):
        """Key function for sorting extern includes"""
        order = ['pe', 'boost']
        if "/" in include:
            lib = include.split('/')[0]
            lib = lib[1:]  # remove '<' character
            if lib in order:
                weight = order.index(lib)
            else:
                weight = 8
        else:
            weight = 9

        return str(weight) + include

    moduleIncludes = sorted(moduleIncludes, key=keyFunction)
    appIncludes = sorted(appIncludes, key=keyFunction)
    externIncludes = sorted(externIncludes, key=externKeyFunction)

    isApp = len(appModules) > 0
    splitModuleIncludes = (len(localIncludes) + len(moduleIncludes)) > 5 and not isApp
    splitExternIncludes = (len(externIncludes)) > 5 and not isApp

    result = []

    if len(appIncludes) > 0:
        result += appIncludes
        result += [""]

    if len(localIncludes + moduleIncludes) > 0:
        if splitModuleIncludes:
            result += splitUpIncludes(localIncludes + moduleIncludes)
        else:
            result += localIncludes + moduleIncludes
        result += [""]

    if len(externIncludes) > 0:
        if splitExternIncludes:
            result += splitUpIncludes(externIncludes)
        else:
            result += externIncludes

    if len(result) > 0 and result[-1] == "":
        del result[-1]

    return result


def splitUpIncludes(includes):
    result = []

    result.append(includes[0])

    for i in range(1, len(includes)):
        cur = includes[i]
        last = includes[i - 1]
        curHasSlash = '/' in cur
        lastHasSlash = '/' in last
        if (not curHasSlash and lastHasSlash) or (curHasSlash and not lastHasSlash):
            result.append("")
        if curHasSlash and lastHasSlash:
            lastModule = last.split('/')[0]
            curModule = cur.split('/')[0]
            if lastModule != curModule:
                result.append("")

        result.append(cur)

    return result


def sortIncludes(parsedCodeFile):
    parsedCodeFile.includes = getSortedIncludes(parsedCodeFile)


def checkIncludeOrder(parsedCodeFile):
    res = getSortedIncludes(parsedCodeFile)

    if res != parsedCodeFile.includes:
        raise CodeStyleException("Includes are not sorted correctly")


# Separator Line Length

def getBodyWithCorrectSeparatorLines(parsedCodeFile, targetLineLength=120, maxLength=80):
    """ Lines that have at least 'maxLength' times the same filling char (*,-,=,/ )
        are extended or shrunk to have exactly targetLineLength  """

    def correctSeparatorLineLength(line, fillChar):

        if fillChar * maxLength in line and '"' in line:
            # print ( parsedCodeFile.filename + " Warning: Fill line with quotation mark detected: " + line)
            return line

        if fillChar * maxLength in line and len(line.rstrip()) != targetLineLength:
            line = line.rstrip()  # remove newline
            if len(line) > targetLineLength:  # remove chars
                numberCharsToDelete = len(line) - targetLineLength
                line = line.replace(fillChar * numberCharsToDelete, "", 1)
            else:  # add chars
                numberCharsToAdd = targetLineLength - len(line)
                line = line.replace(fillChar * maxLength, fillChar * (numberCharsToAdd + maxLength), 1)

            line += "\n"

        return line

    result = []
    for le in parsedCodeFile.body:
        le = correctSeparatorLineLength(le, '*')
        le = correctSeparatorLineLength(le, '-')
        le = correctSeparatorLineLength(le, '=')
        le = correctSeparatorLineLength(le, '/')
        result.append(le)

    return result


def correctSeparatorLines(parsedCodeFile):
    parsedCodeFile.body = getBodyWithCorrectSeparatorLines(parsedCodeFile)


def checkSeparatorLines(parsedCodeFile):
    res = getBodyWithCorrectSeparatorLines(parsedCodeFile)
    if res != parsedCodeFile.body:
        raise CodeStyleException("Not all separator lines have length 120")


# Doxygen Tags

def checkTags(parsedCodeFile):
    tags = parsedCodeFile.tags

    sourceFileName = os.path.basename(parsedCodeFile.filename)
    # currentModule = getModule(parsedCodeFile.filename)

    fileTagFound = False
    # authorTagFound = False
    ingroupTagFound = False
    briefTagFound = False

    for (tag, value) in tags:
        if tag == "file":
            if fileTagFound:
                raise CodeStyleException("Duplicate '\file' tag ")
            fileTagFound = True

            if not value == sourceFileName:
                raise CodeStyleException("'file' tag in header has wrong value " + value)
        elif tag == "author":
            if value not in checkTags.knownAuthors:
                raise CodeStyleException(
                    "Unknown author or email. If valid extend the authors.txt. Invalid Value was: " + value)
        elif tag == "ingroup":
            if ingroupTagFound:
                raise CodeStyleException("Duplicate ingroup tag")
        elif tag == "brief":
            if briefTagFound:
                raise CodeStyleException("Duplicate '\brief' tag")
            briefTagFound = True


checkTags.knownAuthors = []
for line in open(authorfile):
    checkTags.knownAuthors.append(line.split("|")[0].strip())


def correctTags(parsedCodeFile):
    tags = parsedCodeFile.tags
    sourceFileName = os.path.basename(parsedCodeFile.filename)
    currentModule = getModule(parsedCodeFile.filename)

    result = []
    result.append(("file", sourceFileName))
    if len(currentModule) > 0:
        result.append(("ingroup", currentModule))

    for (tag, value) in tags:
        if tag == 'author':
            authorText = re.sub(r'<[^>]*>', "", value)
            if len(authorText.split()) > 2 or "," in value:
                raise CodeStyleException("Multiple Authors in one line?")

            matchingAuthor = ""
            for (author, regexp) in correctTags.authorRegexps:
                if regexp.search(value):
                    matchingAuthor = author
                    break

            if matchingAuthor == "":
                raise CodeStyleException("Unknown Author: " + value)
            else:
                result.append(('author', matchingAuthor))

    for (tag, value) in tags:
        if tag != 'author' and tag != "file" and tag != "ingroup" and value.strip() != "":
            result.append((tag, value))

    parsedCodeFile.tags = result


# Parse authors file
correctTags.authorRegexps = []
for line in open(authorfile):
    splittedLine = line.strip().split("|")
    author = splittedLine[0].strip()
    if len(splittedLine) > 1:
        regexp = re.compile(splittedLine[1].strip(), re.IGNORECASE)
    else:
        regexp = re.compile(splittedLine[0].strip(), re.IGNORECASE)
    correctTags.authorRegexps.append((author, regexp))

# License Check

referenceLicense = [l.replace("\n", "") for l in open(licensefile)]


def checkLicense(parsedCodeFile):
    ref = [l.strip() for l in referenceLicense]
    lic = [l.strip() for l in parsedCodeFile.license]

    if not lic == ref:
        raise CodeStyleException("License is not correct")


def correctLicense(parsedCodeFile):
    parsedCodeFile.license = referenceLicense


# Write Parsed Code File


def writeParsedCodeFile(parsedCodeFile, filename):
    def prependBeforeEachLine(textToPrepend, text):
        l = [textToPrepend + line for line in text.split('\n')]
        return "\n".join(l)

    tags = "\n".join(["\\" + e[0] + " " + e[1] for e in parsedCodeFile.tags])
    doc = "\n".join(parsedCodeFile.doc)

    separatorLine = '//' + (120 - 2) * "="

    result = ""
    result += separatorLine + "\n//\n"

    result += prependBeforeEachLine("//  ", "\n".join(parsedCodeFile.license)) + "\n//\n"
    result += prependBeforeEachLine("//! ", tags)

    if doc.strip() != "":
        result += "\n//!\n"
        result += prependBeforeEachLine("//! ", doc)
        result += "\n//\n"
    else:
        result += "\n//\n"

    result += separatorLine + "\n"

    result += "\n"

    if parsedCodeFile.isHeaderFile:
        result += "#pragma once\n\n"

    for i in parsedCodeFile.includes:
        if len(i.strip()) > 0:
            result += "#include " + i + '\n'
        else:
            result += "\n"

    result += "\n\n"
    result += "".join(parsedCodeFile.body)

    file = open(filename, 'w')
    file.write(result)
    file.close()


# Main

def processFile(filename, args):
    fileExtension = os.path.splitext(filename)[1]

    if not (((fileExtension == ".cpp") or (fileExtension == ".h")) and not filename.endswith(".in.h")):
        return

    if (isDeprecated(filename)):
        return

    parsedCodeFile = None
    try:
        if args.oldStyle:
            parsedCodeFile = OldStyleParsedCodeFile(filename)
        else:
            parsedCodeFile = ParsedCodeFile(filename)

    except ParsingException as exception:
        print(filename + " " + str(exception))

    if parsedCodeFile:
        try:
            checkTags(parsedCodeFile)
            checkLicense(parsedCodeFile)
            checkSeparatorLines(parsedCodeFile)
            checkIncludeOrder(parsedCodeFile)
        except CodeStyleException as exception:
            print(filename + " " + str(exception))

    if args.autoCorrect:
        correctLicense(parsedCodeFile)
        correctSeparatorLines(parsedCodeFile)
        correctTags(parsedCodeFile)
        sortIncludes(parsedCodeFile)
        writeParsedCodeFile(parsedCodeFile, filename)


if __name__ == "__main__":
    import argparse

    argParser = argparse.ArgumentParser()
    argParser.add_argument("target", nargs='?', default=os.getcwd(), help="file or foldername to check")
    argParser.add_argument("-a", "--autoCorrect", action="store_true",
                           help="automatically correct source files if possible")
    argParser.add_argument("-o", "--oldStyle", action="store_true",
                           help="Parse old style headers")
    # argParser.add_argument( "-g", "--graph", help="write depencency graph in json format to specified file" )
    args = argParser.parse_args()

    if os.path.isdir(args.target):
        for root, dirs, files in os.walk(args.target):
            # Skip specified directories
            for d in skipDirectories:
                if d in dirs:
                    dirs.remove(d)

            for filename in files:
                filename = root + "/" + filename
                processFile(filename, args)
    else:
        processFile(os.path.abspath(args.target), args)
