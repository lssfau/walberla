#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import os.path

# Helper script for bash completion

# Call this from .bashrc like this

#    export COMPLETE_SCRIPT="thisScript.py"

#  #cdb goes directly to build directory
#  cdb() {
#    cd /home/me/mybuildDir/$1
#  }
#  # here comes the completion function
# _cdb() {
#   cur="${COMP_WORDS[COMP_CWORD]}"
#   COMPREPLY=(`$COMPLETE_SCRIPT /home/me/mybuildDir $cur`)
# }
# complete -o filenames  -o nospace -F _cdb cdb


# Example:
# ./folderComplete.py /home/bauer/devel/walberla tests/g


# Input:   directory name, and beginning of word
base_dir = sys.argv[1]

if len(sys.argv) == 2:
    to_complete = ""
else:
    to_complete = sys.argv[2]

complete_path = base_dir + "/" + to_complete

# In example valid_path = "/home/bauer/devel/walberla/tests "
#            tail = "g"
valid_path, tail = complete_path.rsplit("/", 1)

# We look for all subdirs in valid_path, only if they start with tail
subdirs = [o for o in os.listdir(valid_path) if os.path.isdir(valid_path + "/" + o) and o.startswith(tail)]

# Make all paths relative to base_dir
result = [os.path.relpath(valid_path + "/" + p, base_dir) + "/" for p in subdirs]

# Print the result
print(" ".join(result))
