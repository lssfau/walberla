from collections import Counter
from argparse import ArgumentParser
import pathlib
import re
import pprint

# WARNING_PATTERN = re.compile(r"\[([^[]+)(?:,-warnings-as-errors)?\]\n")
WARNING_PATTERN = re.compile(r"\[[a-z-]+,-warnings-as-errors\]\n")
TRAILING = len(",-warnings-as-errors]\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("clang_tidy_output", type=str)
    args = parser.parse_args()

    output_fp = pathlib.Path(args.clang_tidy_output)
    clang_tidy_log = output_fp.read_text()
    matches = WARNING_PATTERN.findall(clang_tidy_log)
    matches = [m[1:-TRAILING] for m in matches]
    counter = Counter(matches)

    pprint.pp(counter)
