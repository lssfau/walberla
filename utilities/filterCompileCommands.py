#!/usr/bin/env python3

import json
import sys


def compileCommandSelector(x):
    return not (("extern" in x["file"]) or ("tests" in x["file"]))


def removePrecompiler(x):
    pos = x.find("clang++")
    if pos != -1:
        return x[pos:]
    else:
        return x


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: ./filterCompileCommands.py compile_commands.json")
        exit(-1)

    filename = sys.argv[1]
    print("loading compile commands file: {}".format(filename))

    fin = open(filename, "r")
    cc = json.load(fin)
    fin.close()

    print("compile commands read: {}".format(len(cc)))

    cc_filtered = list(filter(compileCommandSelector, cc))
    for x in cc_filtered:
        x["command"] = removePrecompiler(x["command"])

    print("compile commands filtered: {}".format(len(cc_filtered)))

    fout = open(filename, "w")
    json.dump(cc_filtered, fout)
    fout.close()
