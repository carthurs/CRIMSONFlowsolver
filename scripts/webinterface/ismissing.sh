#!/bin/bash
#script ismissing.sh.  prints 1 if the file is missing, 0 if it exists.
test -e $1
echo $?
