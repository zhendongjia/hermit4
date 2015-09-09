#!/usr/bin/env python
import os
import math


t = [0.01,10]

for j in range(len(t)):

    rp = [i for i in range(10,1000)]
    densp = [0]*len(rp)
    for i in range(len(rp)):
        if rp[i] > 100: densp[i] = 5.5
        else:
            densp[i] = 5.5*rp[i]/100
            if densp[i] < 1.0: densp[i] = 1.0

    #densp = [3.0]*len(rp)
    a_low = 2.0
    a_upp = 3.5
    tdep = 1e6
    
    mp = [4./3 * math.pi * densp[i] * rp[i]**3 * 5e-19 for i in range(len(rp))]

    td1_low = [a_low**2.75*rp[i]*densp[i]*9./1e-2*math.exp(t[j])/(2*math.pi)/tdep for i in range(len(rp))]
    td1_upp = [a_upp**2.75*rp[i]*densp[i]*9./2e-3*math.exp(t[j])/(2*math.pi)/tdep for i in range(len(rp))]
    
    td2_low = [(1./mp[i])*a_low**2*4e-4*math.exp(t[j])/(2*math.pi)/tdep for i in range(len(rp))]
    td2_upp = [(1./mp[i])*a_upp**2*4e-4*math.exp(t[j])/(2*math.pi)/tdep for i in range(len(rp))]


    f = open('t1_t2_t%2.2f' % (t[j]),'w')
    for i in range(len(rp)):
        print >> f, td1_upp[i], td1_low[i], td2_upp[i], td2_low[i], rp[i]
    f.close()    
