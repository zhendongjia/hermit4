#!/usr/bin/env python
import os
import math


file_name1 = 'r_nab'
file_name2 = 'r_ni'

rab = [float(line.split()[0]) for line in open(file_name1).readlines()]

nab =  [int(line.split()[1]) for line in open(file_name1).readlines()]

nab_tot = nab[0]

ri = [float(line.split()[0]) for line in open(file_name2).readlines()]

ni =  [int(line.split()[1]) for line in open(file_name2).readlines()]

ni_tot = ni[0]



robs = [float(line.split()[0]) for line in open(file_name1).readlines()]
nobs = [float(line.split()[1]) for line in open(file_name1).readlines()]

robs_20_200 = []
nobs_20_200 = []
for i in range(len(robs)):
    if robs[i]<20 or robs[i]>200: continue
    else: 
        robs_20_200.append(robs[i])
        nobs_20_200.append(nobs[i])

nobs_20_200_cul = []
nobs_20_200_cul.append(nobs_20_200[0])
for i in nobs_20_200:
    temp = i + nobs_20_200_cul[-1]
    nobs_20_200_cul.append(temp)


robs_20_200.sort()
nobs_20_200_cul.sort()
nobs_20_200_tot = nobs_20_200_cul[0]
print nobs_20_200_cul


f = open(file_name1+'_fra','w')
for i in range(len(nab)):
    print >> f, rab[i], float(nab[i])/nab_tot
f.close()

f = open(file_name2+'_fra','w')
for i in range(len(ni)):
    print >> f, ri[i], float(ni[i])/ni_tot
f.close()

f = open('r_nobs_20_200_cul_fra','w')
for i in range(len(nobs_20_200)):
    print >> f, robs_20_200[i], float(nobs_20_200_cul[i])/nobs_20_200_tot
f.close()


r_i = [float(line.split()[0]) for line in open('ri').readlines()]
r_ab = [float(line.split()[0]) for line in open('rab').readlines()]

m_i = [i**3 for i in r_i]
m_ab = [i**3 for i in r_ab]

m_i_tot = 0
m_ab_tot = 0
for i in m_i: m_i_tot = m_i_tot + i
for i in m_ab: m_ab_tot = m_ab_tot + i
fra = 1- m_ab_tot/m_i_tot
print fra

