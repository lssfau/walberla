import re
import os


class ParsingException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


# Reads the  comment block at the top of a waLBerla source file
#
# The top comment block has to have the following format:
#
#         //===========================================================
#         /*!
#          *  \file   MyHeader.h
#          *  \author Foo Bar <foo.bar@fau.de>
#          *  \author Another author <another.author@fau.de>
#          *  \brief  In this file the great xyz feature is implemented
#          *
#          *  Here comes some additional text
#          *  which describes the file
#          */
#         //============================================================
#
#
#
def parseTopCommentBlock(file):
    separatorRe = re.compile(r'//={20,}')  # A separator line has at least 20 "=" characters
    # Finds pairs like the following
    # \brief some description
    propertiesRe = re.compile(r'\s*\\(file|brief|author|ingroup) (.*)')

    firstLine = file.readline()
    if not separatorRe.match(firstLine):
        raise ParsingException("No Header block at beginning of file")

    tags = []  # tuples of key-value, for example ( 'file', 'MyHeader.h'), (author, "Foo Bar <foo.bar@fau.de>" )
    rest = []  # all additional lines

    # lines = []
    for line in file:
        if separatorRe.match(line):
            break

        replacer = re.compile(r'\/\*!?|\*\/|\*|//')
        line = replacer.sub("", line)

        r = propertiesRe.match(line)
        if r:
            tags.append((r.group(1).strip(), r.group(2).strip()))
        else:
            if line.strip() != "":
                rest.append(line[:-1])

    return (tags, rest)


########################################################################################################################
# Parses the file content after the top comment block
#
#  Allowed content:
#    - pragma once
#    - include directives
#    - legacy: "local includes, module includes, extern includes" comments
#
# Parameter: isHeaderFile  - if true the pragma once has to be present
#
########################################################################################################################
def parseIncludes(file, isHeaderFile, filename):
    includes = []
    content = []
    pragmaOnceFound = False
    headerStopLine = ""
    for line in file:
        cleanLine = line.strip()
        if cleanLine == "":
            continue
        if cleanLine.startswith("#include"):
            includes.append(cleanLine[8:].strip())
        elif cleanLine == "#pragma once":
            if not pragmaOnceFound:
                pragmaOnceFound = True
            else:
                raise ParsingException("Double 'pragma once' found")
        elif cleanLine.startswith("namespace") or 'cond internal' in cleanLine or 'using namespace' in cleanLine:
            headerStopLine = line
            content.append(line)
            break
        else:
            if (cleanLine.lower() != "// local includes") and \
                    (cleanLine.lower() != "// module includes") and \
                    (cleanLine.lower() != "// extern includes") and \
                    (cleanLine.lower() != "// boost includes") and \
                    (cleanLine.lower() != "//boost includes") and \
                    (cleanLine.lower() != "// core includes") and \
                    (cleanLine.lower() != "// local") and \
                    (cleanLine.lower() != "// extern") and \
                    (cleanLine.lower() != "//core") and \
                    (cleanLine.lower() != "// Modules") and \
                    (cleanLine.lower() != "// modules includes") and \
                    (cleanLine.lower() != "//STL") and \
                    (cleanLine.lower() != "// stl includes"):
                headerStopLine = line
                content.append(line)
                break
                # raise ParsingException( "Additional Line at top of file: \n" + cleanLine )

    for line in file:
        content.append(line)
        cleanLine = line.strip()
        if cleanLine.startswith("#include"):
            includedFile = cleanLine[8:].strip()
            includedFile = includedFile[1:-1]
            if not includedFile.endswith(".impl.h"):
                print(filename + " Undetected include:" + cleanLine + " Due to line " + headerStopLine)

    if isHeaderFile and not pragmaOnceFound:
        raise ParsingException("No pragma once found")

    return (includes, content)


########################################################################################################################
# Gets a list of includes and the folder of the file where they were found
# the includes are sorted such that includes from the current module come first, then other waLBerla module includes
# and external includes occur last
#
########################################################################################################################
def sortIncludes(includes, sourceFileName):
    # externRe = re.compile(r'<([^>]*)>')
    # internRe = re.compile(r'"([^"]*)"')

    modulesRe = re.compile(r'free_surface/([^/]*)/([^/]*)/')
    r = modulesRe.search(sourceFileName)
    currentModule = r.group(2)

    # List of tuples with (weight, includeFile)
    # lists with lower weight come first
    weightedIncludes = []

    for include in includes:
        if '"' in include:
            if "/" not in include:
                weightedIncludes.append((1, 0, include))  # same folder include -> weighted 1
            else:
                fileName = include.strip()[1:-1]
                splittedName = fileName.split("/")
                module = splittedName[0]
                if (module == currentModule):
                    weightedIncludes.append((2, len(splittedName), include))  # same module include -> weighted 2
                else:
                    weightedIncludes.append((3, len(splittedName), include))  # other module include -> weighted 3
        elif "<" in include:
            if "/" in include:
                fileName = include.strip()[1:-1]
                splittedName = fileName.split("/")
                weightedIncludes.append((4, len(splittedName), include))  # external include with subpath -> weighted 4
            else:
                weightedIncludes.append((5, 0, include))  # external include without subpath -> weighted 5

    sortedList = sorted(weightedIncludes, key=lambda e: str(e[0]) + str(e[1]) + e[2])
    return [e[2] for e in sortedList]


def prependBeforeEachLine(textToPrepend, text):
    l = [textToPrepend + line for line in text.split('\n')]
    return "\n".join(l)


def writeHeader(commentBlockTags, commentBlockRest, includes, walberlaBaseDir, isHeaderFile):
    sourceCheckerDir = walberlaBaseDir + "/utilities/source_checker/"
    license = open(sourceCheckerDir + "license.txt").read()

    tags = "\n".join(["\\" + e[0] + " " + e[1] for e in commentBlockTags])
    rest = "\n".join(commentBlockRest)

    separatorLine = '//================================================================================================'

    result = ""
    result += separatorLine + "\n//\n"

    result += prependBeforeEachLine("//  ", license) + "\n//\n"
    result += prependBeforeEachLine("//! ", tags) + "\n//\n"
    if rest.strip() != "":
        result += prependBeforeEachLine("//! ", rest)
        result += "\n//\n"
    result += separatorLine + "\n"

    result += "\n"

    if isHeaderFile:
        result += "#pragma once\n\n"

    externNewlineAdded = False
    for i in includes:
        if "<" in i and not externNewlineAdded:
            result += "\n"
            externNewlineAdded = True

        result += "#include " + i + '\n'

    return result


def correctSeparatorLineLength(line, fillChar, targetLineLength=120):
    line = line.rstrip()  # remove newline
    maxLength = 80
    if fillChar * maxLength in line and len(line) != targetLineLength:
        if len(line) > targetLineLength:  # remove chars
            numberCharsToDelete = len(line) - targetLineLength
            line = line.replace(fillChar * numberCharsToDelete, "", 1)
        else:  # add chars
            numberCharsToAdd = targetLineLength - len(line)
            line = line.replace(fillChar * maxLength, fillChar * (numberCharsToAdd + maxLength), 1)

    return line


def checkContent(content):
    # targetLineLength = 120
    result = []
    for line in content:
        line = correctSeparatorLineLength(line, '*')
        line = correctSeparatorLineLength(line, '-')
        line = correctSeparatorLineLength(line, '=')
        line = correctSeparatorLineLength(line, '/')
        result.append(line + "\n")

    return result


def cleanupHeaderTags(tags, sourceFileName):
    modulesRe = re.compile(r'free_surface/([^/]*)/([^/]*)/')
    r = modulesRe.search(sourceFileName)
    module = r.group(2)

    sourceFileName = os.path.basename(sourceFileName)

    result = []
    result.append(("file", sourceFileName))
    result.append(("ingroup", module))
    for tuple in tags:
        if tuple[0] == 'author':
            authorText = re.sub(r'<[^>]*>', "", tuple[1])
            if len(authorText.split()) > 2 or "," in tuple[1]:
                raise ParsingException("Multiple Authors in one line?")

            if "Bauer" in tuple[1]:
                result.append(('author', "Martin Bauer <martin.bauer@fau.de>"))
            elif "Schornbaum" in tuple[1]:
                result.append(('author', "Florian Schornbaum <florian.schornbaum@fau.de>"))
            elif "Godenschwager" in tuple[1]:
                result.append(('author', "Christian Godenschwager <christian.godenschwager@fau.de>"))
            elif "Markl" in tuple[1]:
                result.append(('author', "Matthias Markl <matthias.markl@fau.de>"))
            elif "Anderl" in tuple[1]:
                result.append(('author', "Daniela Anderl <daniela.anderl@lstm.uni-erlangen.de>"))
            elif "Staubach" in tuple[1]:
                result.append(('author', "David Staubach <david.staubach@fau.de>"))
            elif "Fattahi" in tuple[1]:
                result.append(('author', "Ehsan Fattahi <ehsan.fattahi@fau.de>"))
            elif "Bogner" in tuple[1]:
                result.append(('author', "Simon Bogner <simon.bogner@fau.de>"))
            elif "Ammer" in tuple[1]:
                result.append(('author', "Regina Ammer <regina.ammer@fau.de>"))
            elif "Feichtinger" in tuple[1]:
                result.append(('author', "Christian Feichtinger"))
            elif "Iglberger" in tuple[1]:
                result.append(('author', "Klaus Iglberger"))
            elif "Donath" in tuple[1]:
                result.append(('author', "Stefan Donath"))
            else:
                raise ParsingException("Unknown Author: " + tuple[1])

    for tuple in tags:
        if tuple[0] != 'author' and tuple[0] != "file" and tuple[1].strip() != "":
            result.append(tuple)

    return result


def isDeprecated(filename):
    file = open(filename)
    content = file.read()
    return "#pragma message" in content and "deprecated" in content


walberlaBaseDir = "/home/bauer/devel/free_surface/"

for root, dirs, files in os.walk(walberlaBaseDir + "tests/"):
    if 'extern' in dirs:
        dirs.remove('extern')
    if 'stencil' in dirs:
        dirs.remove('stencil')
    if 'vof' in dirs:
        dirs.remove('vof')

    for filename in files:
        fileExtension = os.path.splitext(filename)[1]
        isHeader = (fileExtension == ".h") and not filename.endswith(".impl.h")

        isSourceFile = ((fileExtension == ".cpp") or (fileExtension == ".h")) and not filename.endswith(".in.h")

        if not (isSourceFile):
            continue

        filename = root + "/" + filename

        if (isDeprecated(filename)):
            continue

        file = open(filename)

        try:

            (commentBlockTags, commentBlockRest) = parseTopCommentBlock(file)
            commentBlockTags = cleanupHeaderTags(commentBlockTags, filename)

            (includes, content) = parseIncludes(file, isHeader, filename)
            includes = sortIncludes(includes, filename)

            content = checkContent(content)

            file.close()

            # file = open( filename, 'w')
            # newContent = writeHeader( commentBlockTags, commentBlockRest, includes, walberlaBaseDir, isHeader )
            # newContent += "\n\n"
            # newContent += "".join( content )
            # file.write( newContent )

        except ParsingException as exception:
            print(filename + " " + str(exception))
