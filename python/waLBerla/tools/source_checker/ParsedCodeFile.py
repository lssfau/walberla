import re
import sys
import itertools


class ParsingException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class ParsedCodeFile:
    separatorRe = re.compile(r'//={80,}')  # A separator line has at least 80 "=" characters
    propertiesRe = re.compile(r'//! \\(file|brief|author|ingroup) (.*)')

    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename)
        self.isHeaderFile = (filename.endswith(".h") and not filename.endswith(".impl.h"))

        # Parse File
        self._parseHeaderLine()
        self.license = self._parseLicense()
        self.tags = self._parseTags()
        self.doc = self._parseDocumentationAtTop()
        self.includes = self._parseIncludes()
        self.body = self._parseBody()

# Parser Functions

    def _putBackLine(self, line):
        self.file = itertools.chain((line,), self.file)

    def _parseHeaderLine(self):
        firstLine = self.file.readline()
        if not ParsedCodeFile.separatorRe.match(firstLine):
            raise ParsingException(
                "No Header block at beginning of file found. File has to start with //======== in first line")

    def _parseLicense(self):
        parsedLicense = []
        for line in self.file:
            if line.startswith('//') and not line.startswith('//!'):
                parsedLicense.append(line[2:-1].lstrip())
            else:
                self._putBackLine(line)
                break

        while len(parsedLicense) > 0 and parsedLicense[-1] == "":
            del parsedLicense[-1]

        while len(parsedLicense) > 0 and parsedLicense[0] == "":
            del parsedLicense[0]

        return parsedLicense

    def _parseTags(self):
        tags = []
        for line in self.file:
            r = ParsedCodeFile.propertiesRe.match(line)
            if r:
                tags.append((r.group(1).strip(), r.group(2).strip()))
            else:
                self._putBackLine(line)
                break
        return tags

    def _parseDocumentationAtTop(self):
        result = []
        for line in self.file:
            if line.startswith('//') and not line.startswith('//===='):
                if line.startswith('//!'):
                    strippedLine = line[3:].rstrip()
                else:
                    strippedLine = line[2:].rstrip()
                if len(strippedLine) > 0:
                    result.append(strippedLine)
            else:
                if not ParsedCodeFile.separatorRe.match(line):
                    raise ParsingException("Header block at beginning of file does not end with separator line")
                break

        strippedResult = []
        if len(result) > 0:
            whiteSpacesInFront = sys.maxsize
            for line in result:
                if line.strip() == "":
                    continue
                whiteSpacesInFront = min(whiteSpacesInFront, len(line) - len(line.lstrip()))

            for line in result:
                strippedResult.append(line[whiteSpacesInFront:])
        else:
            strippedResult = result

        return strippedResult

    def _parseIncludes(self):
        includes = []
        pragmaOnceFound = False

        for line in self.file:
            cleanLine = line.strip()
            if cleanLine.startswith("#include"):
                includes.append(cleanLine[8:].strip())
            elif cleanLine == "#pragma once":
                if not pragmaOnceFound:
                    pragmaOnceFound = True
                else:
                    raise ParsingException("Double 'pragma once' found")
            elif cleanLine == "":
                if len(includes) > 0 and includes[-1] != "":
                    includes.append("")
                continue
            else:
                self._putBackLine(line)
                break

        if len(includes) > 0 and includes[-1] == "":
            del includes[-1]

        if (self.isHeaderFile and not pragmaOnceFound):
            raise ParsingException("No 'pragma once' in Header file")

        return includes

    def _parseBody(self):
        alreadyWarned = False
        body = []
        for line in self.file:
            body.append(line)
            cleanLine = line.strip()
            if cleanLine.startswith("#include"):
                includedFile = cleanLine[8:].strip()
                includedFile = includedFile[1:-1]
                if not includedFile.endswith(".impl.h") and not alreadyWarned:
                    alreadyWarned = True
                    # print(f"{self.filename} Warning: undetected include {cleanLine} because of line {body[0]}")
        return body


class OldStyleParsedCodeFile:
    separatorRe = re.compile(r'//={80,}')  # A separator line has at least 80 "=" characters
    propertiesRe = re.compile(r'\s*\\(file|brief|author|ingroup) (.*)')
    replacerRe = re.compile(r'\/\*!?|\*\/|\*|//')

    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename)
        self.isHeaderFile = (filename.endswith(".h") and not filename.endswith(".impl.h"))

        # Parse File
        self._parseHeaderLine()
        self._parseTopCommentBlock()
        self.includes = self._parseIncludes()
        self.body = self._parseBody()
        self.license = ""

# Parser Functions

    def _putBackLine(self, line):
        self.file = itertools.chain((line,), self.file)

    def _parseHeaderLine(self):
        firstLine = self.file.readline()
        if not ParsedCodeFile.separatorRe.match(firstLine):
            raise ParsingException(
                "No Header block at beginning of file found. File has to start with //======== in first line")

    def _parseTopCommentBlock(self):
        tags = []
        doc = []
        for line in self.file:
            if OldStyleParsedCodeFile.separatorRe.match(line):
                break

            # replacer = re.compile(r'\/\*!?|\*\/|\*|//')
            line = OldStyleParsedCodeFile.replacerRe.sub("", line)

            r = OldStyleParsedCodeFile.propertiesRe.search(line)
            if r:
                tags.append((r.group(1).strip(), r.group(2).strip()))
            else:
                if line.strip() != "":
                    doc.append(line[:-1])

        self.tags = tags
        self.doc = doc

    def _parseIncludes(self):
        includes = []
        pragmaOnceFound = False

        for line in self.file:
            cleanLine = line.strip()
            if cleanLine.startswith("#include"):
                includes.append(cleanLine[8:].strip())
            elif cleanLine == "#pragma once":
                if not pragmaOnceFound:
                    pragmaOnceFound = True
                else:
                    raise ParsingException("Double 'pragma once' found")
            elif cleanLine == "":
                if len(includes) > 0 and includes[-1] != "":
                    includes.append("")
                continue
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
                    self._putBackLine(line)
                    break

        if len(includes) > 0 and includes[-1] == "":
            del includes[-1]

        if (self.isHeaderFile and not pragmaOnceFound):
            raise ParsingException("No 'pragma once' in Header file")

        return includes

    def _parseBody(self):
        alreadyWarned = False
        body = []
        for line in self.file:
            body.append(line)
            cleanLine = line.strip()
            if cleanLine.startswith("#include"):
                includedFile = cleanLine[8:].strip()
                includedFile = includedFile[1:-1]
                if not includedFile.endswith(".impl.h") and not alreadyWarned:
                    alreadyWarned = True
                    # print(f"{self.filename} Warning: undetected include {cleanLine} because of line {body[0]}")
        return body
