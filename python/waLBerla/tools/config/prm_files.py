from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np


def _comment_replacer(match):
    start, mid, end = match.group(1, 2, 3)
    if mid is None:
        # single line comment
        return ''
    elif start is not None or end is not None:
        # multi line comment at start or end of a line
        return ''
    elif '\n' in mid:
        # multi line comment with line break
        return '\n'
    else:
        # multi line comment without line break
        return ' '


def _remove_comments(text):
    import re
    comment_re = re.compile(
        r'(^)?[^\S\n]*/(?:\*(.*?)\*/[^\S\n]*|/[^\n]*)($)?',
        re.DOTALL | re.MULTILINE
    )
    return comment_re.sub(_comment_replacer, text)


def __parse_nested_list(value, dtype=float):
    from pyparsing import Word, Group, Forward, OneOrMore, Optional, alphanums, Suppress
    number = Word(alphanums + ".e-").setParseAction(lambda s, l, t: dtype(t[0]))
    arr = Forward()
    element = number | arr
    arr << Group(Suppress('[') + (OneOrMore(element + Optional(Suppress(",")))) + Suppress(']'))
    return arr.parseString(value, parseAll=True).asList()[0]


def __convert_value(value):
    def try_type(thetype, value):
        try:
            return thetype(value)
        except ValueError:
            return None

    val = try_type(int, value)
    if val: return val
    val = try_type(float, value)
    if val: return val

    if value == "True": return True
    if value == "False": return False

    value = value.strip()
    if value[0] == '[' and value[-1] == "]":
        return np.array(__parse_nested_list(value))
    if value[0] == '<' and value[-1] == ">":
        return tuple([__convert_value(s) for s in value[1:-1].split(",")])
    return value


def fromPrm(fileAsString):
    """Parses a prm file to a nested python dict """
    from pyparsing import Word, Group, Forward, Dict, ZeroOrMore, Optional, alphanums, Suppress

    # Grammar
    identifier = Word(alphanums + "_-")
    value = Word(alphanums + ".e-,<>_[] \n").setParseAction(lambda s, l, t: __convert_value(t[0]))
    key_value_pair = Group(identifier + Optional(value, default="") + Suppress(";"))
    block = Forward()
    block_content = Dict(ZeroOrMore(key_value_pair | block))
    block << Group(identifier + Suppress("{") + block_content + Suppress("}"))

    return block_content.parseString(_remove_comments(fileAsString), parseAll=True).asDict()


def __format(value, level, key=""):
    if type(value) is float:
        return "%.10g" % value
    else:
        return str(value).replace("\n", "\n" + "\t" * level + " " * len(key) + " ")


def toPrm(configDict, level=0):
    """Returns a prm string from a nested python dict - Inverse of parse"""
    result = ""
    for key, value in sorted(configDict.items()):
        if type(value) is dict:
            if len(value) == 0:
                result += "\t" * level + str(key) + " {}\n"
            else:
                result += "\t" * level + str(key) + "\n"
                result += "\t" * level + "{\n"
                result += toPrm(value, level + 1)
                result += "\t" * level + "}\n"
        elif (type(value) is list) and len(value) > 0 and type(value[0]) is dict:
            for e in value:
                result += "\t" * level + str(key) + "\n"
                result += "\t" * level + "{\n"
                result += toPrm(e, level + 1)
                result += "\t" * level + "}\n"
        elif type(value) is tuple:
            result += "\t" * level + str(key) + " "
            result += "< " + ",".join([__format(v, level) for v in value]) + " >;\n"
        elif type(value) is float:
            result += "\t" * level + str(key) + " " + __format(value, level) + ";\n"
        else:
            result += "\t" * level + str(key) + " " + __format(value, level, key) + ";\n"

    return result


if __name__ == "__main__":
    test = """
    blockId1 { 
        firstKey 5; // first key
        secondKey 2; // second key
        myFloat 5.2e-4;
        myVector3 <1.0, 2, 3>;
        myVector2 <hi, there>;
        myArray [ [ 1 , 2]
                  [ 3  4] ];
        keyWithoutValue;
        blockId11 {
            nestedKey 25;
            }
    }  
    blockId2 {
        firstKey 25; 
        secondKey 22; 
    }
    outerKey  hallo;
    """

    res = fromPrm(test)
    assert res['blockId1']['firstKey'] == 5
    assert res['blockId1']['myVector3'] == (1.0, 2, 3)
    assert res['blockId1']['keyWithoutValue'] == ''
    assert (res['blockId1']['myArray'] == np.array([[1, 2], [3, 4]])).all()
