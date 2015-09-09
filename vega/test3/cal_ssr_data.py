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

def Gjs(ai, aj, mi, mj):
  alpha = min(ai, aj) / max(ai, aj)
  alpha = ai/aj
  b2 = bsn(1.5, 2, alpha)
  ni = ai**(-1.5)  
  nj = aj**(-1.5)
  gjs = - b2 / (4 * aj**3)
  gjs = gjs * (mi * mj / (ni * nj))**0.5
#  print gjs 
  return gjs




m_j = 0.000954786
m_s = 0.000285837
a_j = 5.202545
a_s = 9.554841



dens0 = 1700 * 1.125e-7
k = 1.5
z_k = 1.094
tdep = 1e6
a_in = 4.5
a_out = 11.0


g_js = Gjs(a_j, a_s, m_j, m_s)
g_sj = Gjs(a_s, a_j, m_s, m_j)

#g_sj = Gij(a_s, a_j, m_j)
#g_js = Gij(a_j, a_s, m_s)
#g_js = -9.63435e-4/360.
#g_sj = -6.09908e-3/360.

#--------------------------------------------------------------------------------------------
#ap = [i*0.01 for i in range(10,200)]
#nstep = 190

ap = [0.63, 1.85]
nstep = 2


f = open('ap_expt','w')
f1 = open('ap_expt_wg','w')
f2 = open('ap_expt_sat','w')
f3 = open('ap_expt_wg_sat','w')

for i in range(nstep):

  g_jd = Gjd(dens0, k, z_k, a_j)
  g_sd = Gjd(dens0, k, z_k, a_s)
  g_pd = Gjd(dens0, k, z_k, ap[i])

  g_jd_wg = Gjd_wg(dens0, dens0, k, a_j, a_j, a_j)
  g_sd_wg = Gjd_wg(dens0, dens0, k, a_s, a_s, a_s)

  g_pj = Gij(ap[i], a_j, m_j)
  g_ps = Gij(ap[i], a_s, m_s)

  g_j = g_js
  g_s = g_sj
  g_p = g_ps + g_pj
  
  g_5 = 0.5*(g_j + g_s + ((g_j - g_s)**2 + 4*g_js**2)**0.5)
  g_6 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)


  #if abs((g_5 - g_p)/g_5) < 1e-1: print ap[i] 
  print 'g5=',g_5, 'g6=',g_6,'gp=',g_p, 'gj=', g_j, 'gs=', g_s, ap[i]
  
  exp_t = (g_pj + g_ps - g_js)/(g_jd - g_pd)
  t = math.log(abs(1./exp_t))*tdep

  exp_t_wg = (g_pj + g_ps - g_js)/(g_jd_wg - g_pd)
  t_wg = math.log(abs(1./exp_t_wg))*tdep  

  exp_t_sat = (g_pj + g_ps - g_sj)/(g_sd - g_pd)
  t_sat = math.log(abs(1./exp_t_sat))*tdep
  
  exp_t_wg_sat = (g_pj + g_ps - g_sj)/(g_sd_wg - g_pd)
  t_wg_sat = math.log(abs(1./exp_t_wg_sat))*tdep

  exp_t4 = 6e-3 * ap[i]**4
   
  print >> f, ap[i], exp_t, exp_t4
  print >> f1, ap[i], exp_t_wg
  print >> f2, ap[i], exp_t_sat
  print >> f3, ap[i], exp_t_wg_sat

  #print 'gjd=', g_jd, 'gjd_wg', g_jd_wg, g_jd/g_jd_wg, 'gsd=', g_sd, 'gsd_wg', g_sd_wg, g_sd/g_sd_wg 
  #print 'ap=', ap[i],'gpj+gps-gjs=', g_pj + g_ps - g_js, 'expt=', exp_t, 'gpj=', g_pj, 'gps=', g_ps,'gjs=', g_js, z_k*math.pi*dens0*exp_t*(1./ap[i] - 1./a_j)
  
f.close()
f1.close()
f2.close()
f3.close()
