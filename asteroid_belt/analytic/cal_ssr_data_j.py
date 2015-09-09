#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
from scipy.integrate import quad


def integrand(x, s, n, alpha):
  r = math.cos(n*x) / (1 + alpha**2 - 2*alpha*math.cos(x))**s
  r = r/math.pi
  return r

def bsn(s, n, alpha):
  b = quad(integrand, 0 , 2*math.pi, args = (s, n, alpha))
  bf = b[0]
  return bf


def Gjd(dens, kg, zk, aj):
  gjd = - zk * math.pi * dens * aj**(-kg)
  nj = aj**(-1.5) 
  gjd = gjd / (nj * aj)
  return gjd


def Gjd_wg(dens_in, dens_out, kg, aj, ain, aout):
  nj = aj**(-1.5)
  gjd = math.pi / (nj * aj)
  A1 = 1./4
  A2 = 9./64
  A3 = 25./156
  gjd1 = gjd * 3 * A1 * (dens_in * (ain/aj)**(4-kg)/(4-kg) + dens_out * (aj/aout)**(1+kg)/(1+kg))
  gjd2 = gjd * 10 * A2 * (dens_in * (ain/aj)**(6-kg)/(6-kg) + dens_out * (aj/aout)**(3+kg)/(3+kg))
  gjd3 = gjd * 21 * A3 * (dens_in * (ain/aj)**(8-kg)/(8-kg) + dens_out * (aj/aout)**(5+kg)/(5+kg))
  gjd = gjd1 + gjd2 + gjd3
  return gjd


def Gij(ai, aj, mj):
  alpha = min(ai, aj) / max(ai, aj)
  b1 = bsn(1.5, 1, alpha)
  ni = ai**(-1.5)
  gij = mj * alpha * b1 
  gij = gij / (4 * ni * ai**2 * max(ai, aj))
  return gij



m_j = 0.000954786
a_j = 5.202545



dens0 = 1700 * 1.125e-7
k = 1.5
z_k = 1.094
tdep = 1e6
a_in = 4.5
a_out = 11.0

#--------------------------------------------------------------------------------------------
ap = [i*0.01 for i in range(10,400)]
nstep = 390

f = open('ap_expt_t','w')

for i in range(nstep):

  g_pd = Gjd(dens0, k, z_k, ap[i])

  g_jd_wg = Gjd_wg(dens0, dens0, k, a_j, a_in, a_out)

  g_pj = Gij(ap[i], a_j, m_j)
  
  exp_t = g_pj/(g_jd_wg - g_pd)
  t = math.log(abs(1./exp_t))*tdep
   
  print >> f, ap[i], exp_t, t

f.close()
