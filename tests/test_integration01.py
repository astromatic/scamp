#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import math,os,subprocess,sys
import numpy as np
from astropy.io.votable import parse

conf = "test_integration01.conf"
xml  = "test_integration01.xml"

# Exit status 0 (= OK) by default
status = 0

def scamp_meanrms(dx):
  return math.sqrt((dx[0]*dx[0] + dx[1]*dx[1]) / 2.0)

def scamp_test(str, passed):
  global status
  if passed:
    print("%-24s ... passed" %str)
  else:
    status += 1
    print("%-24s ... failed" %str)
  return passed

catlist = ["extra/744331p.cat", "extra/744332p.cat", "extra/744333p.cat", "extra/744334p.cat", \
           "extra/744335p.cat", "extra/744336p.cat", "extra/744337p.cat"]

# Remove any preexisting XML file
try:
  os.remove(xml)
except OSError:
  pass

print("Running Test #1 ...")
ret = subprocess.call(["../src/scamp", "-c", conf] + catlist)
scamp_test("SCAMP execution", ret == 0)


# Parse VOTable
votable = parse(xml)

# Fields table
fields = votable.get_table_by_id("Fields")

# Get AS and XY contrasts
contrast_as_min = np.amin(fields.array["AS_Contrast"])
contrast_xy_min = np.amin(fields.array["XY_Contrast"])
scamp_test("Minimum AS contrast: " + str(contrast_as_min) + " expected: > 22.00", contrast_as_min > 22.00)
scamp_test("Minimum XY contrast: " + str(contrast_xy_min) + " expected: > 22.90", contrast_xy_min > 22.90)

# FGroups table
fgroups = votable.get_table_by_id("FGroups")

# Test number of matched detections
nmatched_int = fgroups.array["AstromNDets_Internal"][0]
scamp_test("# of matched detections: " + str(nmatched_int) + " expected: > 180100", nmatched_int > 180100)

# Test number of matched reference sources
nmatched_ref = fgroups.array["AstromNDets_Reference"][0]
scamp_test("# of matched references: " + str(nmatched_ref) + " expected: > 1640", nmatched_ref > 1640)

# Test mean internal RMS errors
sigma_int = scamp_meanrms(fgroups.array["AstromSigma_Internal_HighSN"][0])
scamp_test("Internal dispersion: " + str(sigma_int) + " expected: < 0.035", sigma_int < 0.035)

# Test mean reference RMS errors
sigma_ref = scamp_meanrms(fgroups.array["AstromSigma_Reference_HighSN"][0])
scamp_test("Reference dispersion: " + str(sigma_ref) + " expected: < 0.097", sigma_ref < 0.097)

# Exit with status and message
if status == 0:
  print("All OK.")
else:
  print("%d errors." %status)
sys.exit(status)

