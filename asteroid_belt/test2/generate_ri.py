#!/usr/bin/env python
import random
import math
import os
import shutil
import sys

def y(x,k):
    return x**(-k)/(-k)
 
def g(x,k):
    return (x*(-k))**(-1.0/k)

def bin(x, binsize):
  return 10 ** (int(math.log10(x)/binsize + 0.5) * binsize)

xmin = 10**1.3 + 1e-9
xmax = 10**2.3 - 1e-9
k = 0

f = open('ran','w')
for i in range(1000):

    r = random.uniform(0, 1)


    if k > 0: 
	temp = r*y(xmin,k) + (1-r)*y(xmax,k)
        radius = g(temp,k)
    else:
	radius = 10**random.uniform(1.3,2.3)
    print >> f, radius, bin(radius, 0.1)

f.close()

