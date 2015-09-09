#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
from scipy.integrate import quad


ns = 50
A = [0]*ns
const = [0]*ns
for i in range(len(A)):
  A[i] = float(math.factorial(2*i))/2**(2*i)/(math.factorial(i))**2
  A[i] = A[i]**2
  const[i] = i * (2*i + 1)

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


def Gpd_wg(dens, kg, zk, aj, ain):
  nj = aj**(-1.5)
  dens = dens * aj**(-kg)
  gpd = - math.pi * dens / (nj*aj)
  n_tot = 0
  for i in range(len(A)):
    n =  2 * const[i] * A[i] * (aj/ain)**(2*i-1+kg)/(2*i-1+kg)
    n_tot = n_tot + n
  gpd = gpd * (zk + n_tot) 
  return gpd

def Gjd_wg(dens, kg, aj, ain, aout):
  nj = aj**(-1.5)
  dens = dens * aj**(-kg)
  gjd = 2 * math.pi * dens / (nj * aj)
  n_tot = 0
  for i in range(len(A)):
    n = const[i] * A[i] * ((aj/aout)**(2*i-1+kg)/(2*i-1+kg) + (ain/aj)**(2*i+2-kg)/(2*i+2-kg))
    n_tot = n_tot + n
  gjd = gjd * n_tot
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
  b2 = bsn(1.5, 2, alpha)
  ni = ai**(-1.5)*(1.0 + mi)  
  nj = aj**(-1.5)*(1.0 + mj)
  gjs =  b2 / (4 * aj**3)
  gjs = gjs * (mi * mj / (ni * nj))**0.5
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

g_js = Gij(a_j, a_s, m_s)
g_sj = Gij(a_s, a_j, m_j)

open_v5 = 0
open_gj = 0
open_v5_we = 0
open_v6_we = 1
open_v5_wg = 0
open_gj_wg = 0
open_v5_wg_we = 0
open_v6_wg_we = 0
open_v5_png = 0
open_v5_png_we = 0
open_v6_png_we = 1
open_gj_png = 0

##############################################################################
if open_v5 == 1:
  
  f = open('v5','w')

  ap = [i*0.001 for i in range(500,4500)]
  for i in range(len(ap)):
    t = [10**(j*0.001-4) for j in range(0,5000)]
    for j in range(len(t)):
      
      g_jd = Gjd(dens0, k, z_k, a_j)
      g_sd = Gjd(dens0, k, z_k, a_s)
      g_pd = Gjd(dens0, k, z_k, ap[i])

      g_pj = Gij(ap[i], a_j, m_j)
      g_ps = Gij(ap[i], a_s, m_s)

      g_j = g_js + g_jd*t[j]
      g_s = g_sj + g_sd*t[j]

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t[j]

      print >> f, ap[i], t[j], abs((g_p-g_5))
    
      print ap[i], t[j]

  f.close()

##########33################################################################
if open_v5_we == 1:

  f = open('v5_we','w')

  for i in range(1000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd(dens0, k, z_k, a_j)
      g_sd = Gjd(dens0, k, z_k, a_s)
      g_pd = Gjd(dens0, k, z_k, ap)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_5)/g_5) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()

#############################################################################  
if open_v6_we == 1:

  f = open('v6_we','w')

  for i in range(1000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd(dens0, k, z_k, a_j)
      g_sd = Gjd(dens0, k, z_k, a_s)
      g_pd = Gjd(dens0, k, z_k, ap)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_6 = 0.5*(g_j + g_s + ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_6)/g_6) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()   
#############################################################################  
if open_gj == 1:

  ap = [i*0.01 for i in range(1, 450)]
  nstep = 449
  
  f = open('gj','w')
  for i in range(nstep):

    g_jd = Gjd(dens0, k, z_k, a_j)
    g_sd = Gjd(dens0, k, z_k, a_s)
    g_pd = Gjd(dens0, k, z_k, ap[i])

    g_pj = Gij(ap[i], a_j, m_j)
    g_ps = Gij(ap[i], a_s, m_s)

    exp_t = (g_pj + g_ps - g_js)/(g_jd - g_pd)

    print >> f, ap[i], exp_t, math.log(abs(1.0/exp_t))
    
  f.close()

###############################################################################
if open_v5_wg == 1:
  
  f = open('v5_wg','w')

  ap = [i*0.001 for i in range(500,4500)]
  for i in range(len(ap)):
    t = [10**(j*0.001-4) for j in range(0,5000)]
    for j in range(len(t)):

      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out )
      g_pd = Gpd_wg(dens0, k, z_k, ap[i], a_in)

      g_pj = Gij(ap[i], a_j, m_j)
      g_ps = Gij(ap[i], a_s, m_s)

      g_j = g_js + g_jd*t[j]
      g_s = g_sj + g_sd*t[j]

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t[j]

      print >> f, ap[i], t[j], abs(g_p - g_5)
      print ap[i], t[j]

  f.close()

########################################################################
if open_v5_wg_we == 1:

  f = open('v5_wg_we','w')

  for i in range(2000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
      g_pd = Gpd_wg(dens0, k, z_k, ap, a_in)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_5)/g_5) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()

if open_v6_wg_we == 1:

  f = open('v6_wg_we','w')

  for i in range(1000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
      g_pd = Gpd_wg(dens0, k, z_k, ap, a_in)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_6 = 0.5*(g_j + g_s + ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_6)/g_6) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()  
  
########################################################################
if open_gj_wg == 1:

  ap = [i*0.01 for i in range(100, 450)]
  nstep = 350
  
  f = open('gj_wg','w')
  for i in range(nstep):

    g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
    g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
    g_pd = Gpd_wg(dens0, k, z_k, ap[i], a_in)

    g_pj = Gij(ap[i], a_j, m_j)
    g_ps = Gij(ap[i], a_s, m_s)

    exp_t = (g_pj + g_ps - g_js)/(g_jd - g_pd)

    print >> f, ap[i], exp_t, math.log(abs(1.0/exp_t))
    
  f.close()

#############################################################################
if open_v5_png == 1:

  f = open('v5_png','w')

  ap = [i*0.001 for i in range(500,4500)]
  for i in range(len(ap)):
    t = [10**(j*0.001-4) for j in range(0,5000)]
    for j in range(len(t)):

      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
      g_pd = Gjd(dens0, k, z_k, ap[i])

      g_pj = Gij(ap[i], a_j, m_j)
      g_ps = Gij(ap[i], a_s, m_s)

      g_j = g_js + g_jd*t[j]
      g_s = g_sj + g_sd*t[j]

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t[j]

      print >> f, ap[i], t[j], abs(g_p - g_5)
      print ap[i], t[j]

  f.close()

############################################################################
if open_v5_png_we == 1:

  f = open('v5_png_we','w')

  for i in range(1000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
      g_pd = Gjd(dens0, k, z_k, ap)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_5 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_5)/g_5) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()

#############################################################################
if open_v6_png_we == 1:

  f = open('v6_png_we','w')

  for i in range(1000):
    while True:
      ap = random.uniform(0.5, 4.5)
      t = random.uniform(0,5)
      t = 10**(t-4)
      g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
      g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
      g_pd = Gjd(dens0, k, z_k, ap)

      g_pj = Gij(ap, a_j, m_j)
      g_ps = Gij(ap, a_s, m_s)

      g_j = g_js + g_jd*t
      g_s = g_sj + g_sd*t

      g_6 = 0.5*(g_j + g_s + ((g_j - g_s)**2 + 4*g_js**2)**0.5)

      g_p = g_ps + g_pj + g_pd*t

      if abs((g_p - g_6)/g_6) < 1e-3:
        print >> f, ap, t, math.log(abs(1.0/t))
        print ap, t
        break
      
  f.close()   
  
#############################################################################
if open_gj_png == 1:

  ap = [i*0.01 for i in range(100, 450)]
  nstep = 350
  
  f = open('gj_png','w')
  for i in range(nstep):

    g_jd = Gjd_wg(dens0, k, a_j, a_in, a_out)
    g_sd = Gjd_wg(dens0, k, a_s, a_in, a_out)
    g_pd = Gjd(dens0, k, z_k, ap[i])

    g_pj = Gij(ap[i], a_j, m_j)
    g_ps = Gij(ap[i], a_s, m_s)

    exp_t = (g_pj + g_ps - g_js)/(g_jd - g_pd)

    print >> f, ap[i], exp_t, math.log(abs(1.0/exp_t))

  f.close()
