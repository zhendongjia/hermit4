#!/usr/bin/env python
import os
import math

rab = [float(line.split()[0]) for line in open('rab').readlines()]
ri = [float(line.split()[0]) for line in open('ri').readlines()]

mab_tot = 0.0
for i in rab:
    mab_tot = mab_tot + i**3

mi_tot = 0.0
for i in ri:
    mi_tot = mi_tot + i**3

md = 1.0 - mab_tot/mi_tot

print md

