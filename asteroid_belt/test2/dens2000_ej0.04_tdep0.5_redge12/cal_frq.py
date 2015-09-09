#!/usr/bin/env python
import os
import math

ri = [float(line.split()[0]) for line in open('r_fi').readlines()]
ni = [float(line.split()[1]) for line in open('r_fi').readlines()]
ni_tot = 0
for i in ni: ni_tot = ni_tot + i

rab = [float(line.split()[0]) for line in open('r_fab').readlines()]
nab = [float(line.split()[1]) for line in open('r_fab').readlines()]
nab_tot = 0
for i in nab: nab_tot = nab_tot + i

f_ri = [i/ni_tot for i in ni]
f_rab = [i/nab_tot for i in nab]

print f_ri, f_rab

frq = [f_rab[i]/f_ri[i] for i in range(len(ri))]

f = open('r_f','w')
for i in range(len(ri)):
	print >> f, ri[i], frq[i]
f.close() 
