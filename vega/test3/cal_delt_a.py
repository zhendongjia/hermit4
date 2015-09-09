#!/usr/bin/env python
import random
import math
import os
import shutil
import sys

mj = 0.000954786
ms = 2.0

f = open('aj_ej_b','w')
for i in range(1000):
    while True:
        aj = random.uniform(10,120)
        ej = random.uniform(0.0, 0.8)
        delt_a = 1.8 * aj * ej**0.2 * (mj/ms)**0.2
        a_in = aj - delt_a 
        a_out = aj + delt_a
        if abs(a_out - 110) < 0.1:
            print >> f, aj, ej
            break
        if abs(a_in - 14.) < 0.1:
            print >>f, aj, ej
            break
f.close()
