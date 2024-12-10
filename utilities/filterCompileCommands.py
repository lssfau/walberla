#!/usr/bin/env python3

import argparse
import pathlib
import json
import sys


class QualifiedSequence(argparse.Action):
    """Append qualified values from different arguments into the same destination."""
    def __call__(self, parser, namespace, values, option_string=None):
        accumulator = getattr(namespace, self.dest, None) or []
        assert option_string is not None
        mode = "include" if option_string in ("-i", "--include") else "exclude"
        accumulator.append((mode, values))
        setattr(namespace, self.dest, accumulator)


parser = argparse.ArgumentParser(description="Filter out source files from CMake database.")
parser.add_argument("-f", "--file", action="store", type=str, required=True,
                    help="Database file to edit")
parser.add_argument("-i", "--include", action=QualifiedSequence, dest="filters",
                    nargs="+", help="Include paths containing these folder names")
parser.add_argument("-e", "--exclude", action=QualifiedSequence, dest="filters",
                    nargs="+", help="Exclude paths containing these folder names")


def compileCommandSelector(x, filters=None):
    if filters is None:
        filters = [("exclude", ("extern", "tests"))]
    path = "/".join(pathlib.Path(x["file"]).parts)[1:]
    keep = True
    for mode, components in filters:
        for component in components:
            subpath = "/".join(("", ) + pathlib.Path(component).parts + ("", ))
            if subpath in path or component == "*":
                keep = (mode == "include")
                break
    return keep


def removePrecompiler(x):
    pos = x.find("clang++")
    if pos != -1:
        return x[pos:]
    else:
        return x


if __name__ == "__main__":
    args = parser.parse_args()

    print(f"loading compile commands file: {args.file}")

    with open(args.file, "r") as f:
        cc = json.load(f)

    print(f"compile commands read: {len(cc)}")

    cc_filtered = list(filter(lambda x: compileCommandSelector(x, args.filters), cc))
    for x in cc_filtered:
        x["command"] = removePrecompiler(x["command"])

    print(f"compile commands filtered: {len(cc_filtered)}")

    with open(args.file, "w") as f:
        json.dump(cc_filtered, f, indent=2)
