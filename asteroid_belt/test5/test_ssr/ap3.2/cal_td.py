x#!/usr/bin/env python
import os
import math


a = 3.2
rp = 100
if rp > 100: densp = 5.5
else:
    densp = 5.5*rp/100
    if densp < 1.0: densp = 1.0
mp = 4./3 * math.pi * densp * rp**3 * 5e-19

td1_upp = a**2.75*rp*densp*9./1e-3
td1_low = a**2.75*rp*densp*9./1e-2
td2 = (1./mp)*a**2*4e-4


tdep = 1e6
t = [i*1e3 for i in range(0,10000)]
t1_upp = [td1_upp*math.exp(i/tdep)/tdep for i in t]
t1_low = [td1_low*math.exp(i/tdep)/tdep for i in t]
t2 = [td2*math.exp(i/tdep)/tdep for i in t]

f = open('t1_t2_t_rp100','w')
for i in range(len(t)):
    print >> f, t1_upp[i], t1_low[i], t2[i], t[i]
f.close()
