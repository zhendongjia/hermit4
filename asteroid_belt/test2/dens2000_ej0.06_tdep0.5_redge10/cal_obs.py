#!/usr/bin/env python
import os
import math

file_name = 'asteroids_data'

d = [float(line.split()[0]) for line in open(file_name).readlines()]
n = [float(line.split()[1]) for line in open(file_name).readlines()]

f = open('r_nobs','w')
for i in range(len(n)):
    print >> f, d[i]/2, n[i]
f.close()
 
n_cul = 0
f = open('r_nobs_cul','w')
for i in range(len(n)):
    n_cul = n_cul + n[i]
    print >> f, d[i]/2, n_cul
f.close()

n_20_200 = []
r_20_200 = []
n_20_200_tot = 0
f = open('r_fobs_20_200','w')
for i in range(len(n)):
    if d[i]/2 > 19.0 and d[i]/2 < 200.0:
       n_20_200.append(n[i])
       r_20_200.append(d[i]/2)
       n_20_200_tot = n_20_200_tot + n[i] 
for i in range(len(n_20_200)):
    print >> f, r_20_200[i], float(n_20_200[i])/n_20_200_tot
f.close()



