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




  
def Fpd_wg(dens, kg, zk, aj, aout):
  dens = dens * aj**(-kg)
  dens_base = dens * (-kg)
  fr = - 4 * math.pi * dens
  fr_base = -4 * math.pi * dens_base
  n_tot = 0
  n_tot_base = 0
  for i in range(len(A)):
    n = i * A[i] * (aj/aout)**(2*i-1+kg)/(2*i-1+kg)
    n_base = n * (2*i-1+kg)
    n_tot = n_tot + n
    n_tot_base = n_tot_base + n_base
  fr = fr * (zk + n_tot)
  fr_base = fr_base * (zk + n_tot) + fr * n_tot_base
#  fr_dot = fr_base * (aj_dot/aj)
#  fpd = fr * exp(-t/td)
#  fpd_dot = fr_dot * exp(-t/td) + fr * exp(-t/td) * (-1.0/td)/(2*math.pi)
#  fpd_dot = fr_base * (aj_dot/aj) * exp(-t/td) + fr * exp(-t/td) * (-1.0/td)/(2*math.pi)
#  print aj, fr, fr_base
  return fr, fr_base




def Fjd_wg(dens, kg, aj, ain, aout):
  dens = dens * aj**(-kg) 
  dens_base = dens * (-kg)
  fr = 2 * math.pi * dens
  fr_base = 2 * math.pi * dens_base
  n_tot = 0
  n_tot_base = 0
  for i in range(len(A)):
    n = A[i] * (2*i*(aj/aout)**(2*i-1+kg)/(2*i-1+kg) - (2*i+1)*(ain/aj)**(2*i+2-kg)/(2*i+2-kg))
    n_base = A[i] * (2*i*(aj/aout)**(2*i-1+kg) + (2*i+1)*(ain/aj)**(2*i+2-kg))
    n_tot = n_tot + n
    n_tot_base = n_tot_base + n_base
  fr = fr * n_tot
  fr_base = fr_base * n_tot + fr * n_tot_base
#  fr_dot = fr_base * (aj_dot/aj)
#  fjd = fr * exp(-t/td)
#  fjd_dot = fr_dot * exp(-t/td) + fr * exp(-t/td) * (-1.0/td)/(2*math.pi)
#  fjd_dot = fr_base * (aj_dot/aj) * exp(-t/td) + fr * exp(-t/td) * (-1.0/td)/(2*math.pi)
#  print aj, fr, fr_base
  return fr, fr_base



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

open_wg = 1
##############################################################################
if open_wg == 1:
  
  f = open('ap_f_fdot','w')
  ap = [i*0.0001 for i in range(1000,45000)]
  for i in range(len(ap)):
      fp,fp_dot = Fpd_wg(dens0, k, ap[i], a_in, a_out)
      print >> f, ap[i], fp, fp_dot
  f.close()

  f = open('aj_f_fdot','w')
  apj = [a_j + (i-5000)*0.0001 for i in range(10000)]
  for i in range(len(apj)):
      fj,fj_dot = Fjd_wg(dens0, k, apj[i], a_in, a_out)
      print >> f, apj[i], fj, fj_dot
  f.close()

  f = open('as_f_fdot','w')
  aps = [a_s + (i-5000)*0.0001 for i in range(10000)]
  for i in range(len(aps)):
      fj,fj_dot = Fjd_wg(dens0, k, aps[i], a_in, a_out)
      print >> f, aps[i], fj, fj_dot
  f.close()
