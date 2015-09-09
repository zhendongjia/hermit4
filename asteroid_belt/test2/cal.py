#!/usr/bin/env python
import os
import math


file_name = 'asteroids_data1'

rab_obs = [float(line.split()[0]) for line in open(file_name).readlines()]
nab_obs = [float(line.split()[1]) for line in open(file_name).readlines()]
nab_obs_tot = nab_obs[-1] 


f = open(file_name+'_fra','w')
for i in range(len(nab_obs)):
    print >> f, rab_obs[i]/2, float(nab_obs[i])/nab_obs_tot
f.close()



