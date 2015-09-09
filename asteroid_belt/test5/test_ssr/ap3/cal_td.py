#!/usr/bin/env python
import os
import math


a = 3.0
rp = [10, 50, 100, 200, 500, 1000]
for i in range(len(rp)):
    if rp[i] > 100: densp = 5.5
    else:
        densp = 5.5*rp[i]/100
        if densp < 1.0: densp = 1.0
        
    mp = 4./3 * math.pi * densp * rp[i]**3 * 5e-19

    td1_upp = a**2.75*rp[i]*densp*9./1e-3
    td1_low = a**2.75*rp[i]*densp*9./1e-2
    td2 = (1./mp)*a**2*4e-4

    td1_upp = td1_upp/(2*math.pi)
    td1_low = td1_low/(2*math.pi)
    td2 = td2/(2*math.pi)
    
    tdep = 1e6
    t = [j*1e3 for j in range(0,10000)]
    t1_upp = [td1_upp*math.exp(j/tdep)/tdep for j in t]
    t1_low = [td1_low*math.exp(j/tdep)/tdep for j in t]
    t2 = [td2*math.exp(j/tdep)/tdep for j in t]

    f = open('t1_t2_t_rp%d' % (rp[i]),'w')
    for j in range(len(t)):
        print >> f, t1_upp[j], t1_low[j], t2[j], t[j]
    f.close()


#####################################################
    densp = 3.0
    mp = 4./3 * math.pi * densp * rp[i]**3 * 5e-19

    td1_upp = a**2.75*rp[i]*densp*9./1e-3
    td1_low = a**2.75*rp[i]*densp*9./1e-2
    td2 = (1./mp)*a**2*4e-4

    td1_upp = td1_upp/(2*math.pi)
    td1_low = td1_low/(2*math.pi)
    td2 = td2/(2*math.pi)

    tdep = 1e6
    t = [j*1e3 for j in range(0,10000)]
    t1_upp = [td1_upp*math.exp(j/tdep)/tdep for j in t]
    t1_low = [td1_low*math.exp(j/tdep)/tdep for j in t]
    t2 = [td2*math.exp(j/tdep)/tdep for j in t]


    f = open('t1_t2_t_densp3_rp%d' % (rp[i]),'w')
    for j in range(len(t)):
        print >> f, t1_upp[j], t1_low[j], t2[j], t[j]
    f.close()
    
