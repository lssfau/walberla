#!/usr/bin/env python3

import os

error = False
for root, dirnames, filenames in os.walk(".."):
    for filename in filenames:
        if filename.endswith(".h") and not filename.endswith(".impl.h"):
            if "extern" not in root:
                file = os.path.join(root, filename)
                if "#pragma once" not in open(file).read():
                    print(file)
                    error = True

if error:
    exit(-1)
