#!/bin/bash

REGEX='^((#include)|(void)|(uint_t))'
cd default_codegen
diff -u -B <(tail -n +20 FromKernelPackInfoPull.cpp | grep -vP "$REGEX")  <(tail -n +20 AccessorBasedPackInfoEven.cpp | grep -vP "$REGEX") || exit 1
diff -u -B <(tail -n +20 FromKernelPackInfoPush.cpp | grep -vP "$REGEX")  <(tail -n +20 AccessorBasedPackInfoOdd.cpp | grep -vP "$REGEX") || exit 1
