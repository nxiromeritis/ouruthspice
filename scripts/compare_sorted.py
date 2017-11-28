#!/usr/bin/python

# this scripts takes the spice solution outputs and compares the results
# WARNING: the file contents must be lower and sorted (use lc_sort.sh)
# NOTE: line with possible node named g or gnd is ignored

import sys
from itertools import izip
from math import fabs


if len(sys.argv) != 3:
    print "Error: Invalid number of arguments."
    print "Use: " + sys.argv[0] + " <file1> <file2"
    sys.exit()

eps = 10e-10
pass_flag = 1
lines_parsed = 0
lines_failed = 0

with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
    for fline1, fline2 in izip(f1, f2):
        line1 = fline1.split()
        line2 = fline2.split()

        # check if files are identically sorted
        if line1[0] != line2[0]:
            # print "Error. Node names are not identical"
            # print "line1: " + fline1
            # print "line2: " + fline2
            # sys.exit()
            lines_parsed += 1
            lines_failed += 1
            continue

        # files are supposed to be identically sorted
        # bypass whitespace lines
        if len(line1) == 0 or len(line2) == 0:
            continue


        # files are supposed to be identically sorted
        if len(line1) == 1 or len(line2) == 1:
            print "Warning: Line (" + fline1 + ") contains only 1 word"
            continue


        # files are supposed to be identically sorted
        if line1[0] == 'g' or line2[0] == 'g':
            # print "Note: Grounding node found. Bypassing line"
            continue

        lines_parsed += 1

        # convert the values into doubles

        d1 = float(line1[1]);
        d2 = float(line2[1]);

        # debug
        # print "f1: " + line1[0] + " = " + str(d1)
        # print "f2: " + line2[0] + " = " + str(d1)
        # print

        if d1 != 0:
            if (fabs(d1 - d2)*100.0) > d1:
                # print "Error: Node percentage difference larger than 1%"
                # print "f1: " + line1[0] + " = " + str(d1)
                # print "f2: " + line2[0] + " = " + str(d2)
                pass_flag = 0
                lines_failed += 1
        else:
            if fabs(d1 - d2) > eps:
                # print "Error: Node value difference larger than EPS"
                # print "f1: " + line1[0] + " = " + str(d1)
                # print "f2: " + line2[0] + " = " + str(d2)
                pass_flag = 0
                lines_failed += 1


if pass_flag == 1:
    msg = "PASS (failed " + str(lines_failed) + "/" + str(lines_parsed) + ")"
else:
    msg = "FAIL (failed " + str(lines_failed) + "/" + str(lines_parsed) + ")"

print msg + ": " + sys.argv[1] + ", " + sys.argv[2]
