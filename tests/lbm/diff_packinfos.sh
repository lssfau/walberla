#!/bin/bash

REGEX='^((#include)|(void)|(uint_t))'
cd default_codegen
diff -u -B <(grep -vP "$REGEX" FromKernelPackInfoPull.cpp)  <(grep -vP "$REGEX" AccessorBasedPackInfoEven.cpp) || exit 1
diff -u -B <(grep -vP "$REGEX" FromKernelPackInfoPush.cpp)  <(grep -vP "$REGEX" AccessorBasedPackInfoOdd.cpp) || exit 1
