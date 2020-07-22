#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import os.path

DIR1 = sys.argv[1]
DIR2 = sys.argv[2]


def is_in_folder(filename, folder):
    # normalize both parameters
    fn = os.path.normpath(filename)
    fd = os.path.normpath(folder)

    if fn == fd:
        return True

    # get common prefix
    commonprefix = os.path.commonprefix([fn, fd])
    if commonprefix == fd:
        # in case they have common prefix, check more:
        sufix_part = fn.replace(fd, '')
        sufix_part = sufix_part.lstrip('/')
        new_file_name = os.path.join(fd, sufix_part)
        if new_file_name == fn:
            return True
        pass
    # for all other, it's False
    return False


if is_in_folder(os.getcwd(), DIR1):
    rel = os.path.relpath(os.getcwd(), DIR1)
    print(DIR2 + "/" + rel)
elif is_in_folder(os.getcwd(), DIR2):
    rel = os.path.relpath(os.getcwd(), DIR2)
    print(DIR1 + "/" + rel)
else:
    print(os.getcwd())
