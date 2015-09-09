#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
import Gnuplot as gp
gp.GnuplotOpts.default_term = 'post enh solid color eps'
import Gnuplot.funcutils
from scipy import integrate


def integrand(sita, s, n, alpha):
  r = math.cos(n*sita) / (1 + alpha**2 - 2*alpha*math.cos(sita))**s
  r = 2*r / math.pi
  return r

def bsn(s, n, alpha):
  b = integrate.quad(integrand, 0 , math.pi, args = (s, n, alpha))
  return b


def Gjd(dens0, k, zk, aj):
  gjd = - zk * math.pi * dens0 * aj**(-k)
  nj = aj**(-1.5)
  gjd = gjd / (nj * aj)
  return gjd


def Gjd_wg(dens0, k, zk, aj, a_in, a_out):
  nj = aj**(-1.5)
  gjd = math.pi * dens0 / (nj * aj)
  A1 = 1./4
  A2 = 9./64
  gjd1 = gjd * 3 * A1 * ((a_in/aj)**(4-k)/(4-k) + (aj/a_out)**(1+k)/(1+k))
  gjd2 = gjd * 10 * A2 * ((a_in/aj)**(6-k)/(6-k) + (aj/a_out)**(3+k)/(3+k))
  gjd = gjd1 + gjd2
  return gjd

def Gij(ai, aj, mj):
  alpha = min(ai, aj) / max(ai, aj)
  b1 = bsn(3/2, 1, alpha)
  ni = ai**(-1.5)
  gij = mj * alpha * b1[0]
  gij = gij / (4 * ni * ai**2 * max(ai, aj))
  return gij

def Gjs(a_j, a_s, m_j, m_s):
  b2 = bsn(3/2, 2, a_j/a_s)
  n_j = a_j**(-1.5)
  n_s = a_s**(-1.5)
  gjs = - b2[0] * (a_j / a_s) / (4 * a_s**3)
  gjs = gjs * (m_j * m_s / (n_j * n_s))
  return gjs


  

ap = [i*0.01 for i in range(50,400)]
nstep = 350 

m_j = 0.0009546
m_s = 0.0002857
a_j = 5.2
a_s = 9.5

b2 = bsn(3/2, 2, a_j/a_s)

g_js = Gjs(a_j, a_s, m_j, m_s)


# minimum mass nebular model with no gap
dens0 = 1700 * 1.125e-7
k = 1.5
zk = 1.094
tdep = 1e6

g_jd = Gjd(dens0, k, zk, a_j)

g_pd = [None] * nstep
g_pj = [None] * nstep
g_ps = [None] * nstep
exp_t = [None] * nstep
t = [None] * nstep
for i in range(nstep): 
  g_pd[i] = Gjd(dens0, k, zk, ap[i])
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  t[i] = math.log(1/exp_t[i])*tdep

data = [(ap[i], 6e-3*ap[i]**4) for i in range(nstep)]
data_t = [(math.log(1/(6e-3*ap[i]**4))*tdep, ap[i]) for i in range(nstep)]

data1 = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t1 = [(t[i], ap[i]) for i in range(nstep)]
data_gpd1 = [(t[i], Gjd(dens0,k, zk, 3.0)*exp_t[i]) for i in range(nstep)]
data_gpj1 = [(t[i], Gij(3.0, a_j, m_j)) for i in range(nstep)]
data_gjd1 = [(t[i], g_jd*exp_t[i]) for i in range(nstep)]
data_gjs1 = [(t[i], g_js) for i in range(nstep)]
data_gp1 = [(t[i], Gjd(dens0,k, zk, 3.0)*exp_t[i] + Gij(3.0, a_j, m_j)) for i in range(nstep)]
data_gj1 = data_gjd1



# nebular model with gap
dens0 = 1700 * 1.125e-7
k = 1.5
zk = 1.094
tdep = 1e6
a_in = 0.0
a_out = 11.0

g_jd = Gjd_wg(dens0, k, zk, a_j, a_in, a_out)

g_pd = [None] * nstep
g_pj = [None] * nstep
g_ps = [None] * nstep
exp_t = [None] * nstep
t = [None] * nstep
for i in range(nstep): 
  g_pd[i] = Gjd(dens0, k, zk, ap[i])
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  t[i] = math.log(1/exp_t[i])*tdep

data2 = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t2 = [(t[i], ap[i]) for i in range(nstep)]
data_gpd2 = [(t[i], Gjd(dens0,k, zk, 3.0)*exp_t[i]) for i in range(nstep)]
data_gpj2 = [(t[i], Gij(3.0, a_j, m_j)) for i in range(nstep)]
data_gjd2 = [(t[i], g_jd*exp_t[i]) for i in range(nstep)]
data_gjs2 = [(t[i], g_js) for i in range(nstep)]
data_gp2 = [(t[i], Gjd(dens0,k, zk, 3.0)*exp_t[i] + Gij(3.0, a_j, m_j)) for i in range(nstep)]
data_gj2 = data_gjd2


# nebular model with gap
dens0 = 4500 * 1.125e-7
k = 1.0
zk = 1.0
tdep = 1e6
a_in = 0.0
a_out = 11.0

g_jd = Gjd_wg(dens0, k, zk, a_j, a_in, a_out)

g_pd = [None] * nstep
g_pj = [None] * nstep
g_ps = [None] * nstep
exp_t = [None] * nstep
t = [None] * nstep
for i in range(nstep): 
  g_pd[i] = Gjd_wg(dens0, k, zk, ap[i], a_in, a_out)
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  t[i] = math.log(1/exp_t[i])*tdep

data3 = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t3 = [(t[i], ap[i]) for i in range(nstep)]



# nebular model with gap
dens0 = 2000 * 1.125e-7
k = 1.5
zk = 1.0
tdep = 1e6
a_in = 0.0
a_out = 11.0

g_jd = Gjd_wg(dens0, k, zk, a_j, a_in, a_out)

g_pd = [None] * nstep
g_pj = [None] * nstep
g_ps = [None] * nstep
exp_t = [None] * nstep
t = [None] * nstep
for i in range(nstep): 
  g_pd[i] = Gjd_wg(dens0, k, zk, ap[i], a_in, a_out)
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  t[i] = math.log(1/exp_t[i])*tdep

data4 = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t4 = [(t[i], ap[i]) for i in range(nstep)]



# nebular model with gap
dens0 = 1700 * 1.125e-7
k = 1.5
zk = 1.09
tdep = 1e6
a_in = 4.5
a_out = 11.0

g_jd = Gjd_wg(dens0, k, zk, a_j, a_in, a_out)

g_pd = [None] * nstep
g_pj = [None] * nstep
g_ps = [None] * nstep
exp_t = [None] * nstep
t = [None] * nstep
for i in range(nstep): 
  g_pd[i] = Gjd(dens0, k, zk, ap[i])
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  t[i] = math.log(1/exp_t[i])*tdep

data5 = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t5 = [(t[i], ap[i]) for i in range(nstep)]


#print Gjd_wg(dens0, k, zk, a_j, 4.5, 11.0)/Gjd_wg(dens0, k, zk, a_j, 5.0, 11.0)


item = gp.Data(data, title = 'a^4', with_='l lw 10')
item1 = gp.Data(data1, title = 'dens0 = 1700, no gap', with_='l lw 10')
item2 = gp.Data(data2, title = 'dens0 = 1700, with gap, R_in = 0.0', with_='l lw 10')
item3 = gp.Data(data3, title = 'dens0 = 4500, with gap', with_='l lw 10')
item4 = gp.Data(data4, title = 'dens0 = 2000, with gap, kg = 1.5', with_='l lw 10')
item5 = gp.Data(data5, title = 'dens0 = 1700, with gap, R_in = 4.5', with_='l lw 10')


item_t = gp.Data(data_t, title = 'a^4', with_='l lw 10')
item_t1 = gp.Data(data_t1, title = 'dens0 = 1700, no gap', with_='l lw 10')
item_t2 = gp.Data(data_t2, title = 'dens0 = 1700, with gap, R_in = 0.0', with_='l lw 10')
item_t3 = gp.Data(data_t3, title = 'dens0 = 4500, with gap', with_='l lw 10')
item_t4 = gp.Data(data_t4, title = 'dens0 = 2000, with gap, kg = 1.5', with_='l lw 10')
item_t5 = gp.Data(data_t5, title = 'dens0 = 1700, with gap, R_in = 4.5', with_='l lw 10')

item_gpd1 = gp.Data(data_gpd1, title = 'gpd', with_='l lw 10')
item_gpd2 = gp.Data(data_gpd2, title = 'gpd', with_='l lw 10')
item_gpj1 = gp.Data(data_gpj1, title = 'gpj', with_='l lw 10')
item_gpj2 = gp.Data(data_gpj2, title = 'gpj', with_='l lw 10')
item_gjd1 = gp.Data(data_gjd1, title = 'gjd', with_='l lw 10')
item_gjd2 = gp.Data(data_gjd2, title = 'gjd', with_='l lw 10')
item_gjs1 = gp.Data(data_gjs1, title = 'gjs', with_='l lw 10')
item_gjs2 = gp.Data(data_gjs2, title = 'gjs', with_='l lw 10')

item_gp1 = gp.Data(data_gp1, title = 'gp', with_='l lw 10')
item_gj1 = gp.Data(data_gj1, title = 'gj', with_='l lw 10')
item_gp2 = gp.Data(data_gp2, title = 'gp', with_='l lw 10')
item_gj2 = gp.Data(data_gj2, title = 'gj', with_='l lw 10')


grav1 = 1.094 * 4 * math.pi * dens0 * 5.2**(-k)
grav2 = (0.5/(1+k)*(3.0/11.0)**(1+k) + 0.5625/(3+k)*(3.0/11.0)**(3+k)) * dens0 * 5.2**(-k)

print grav1/grav2

g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "ap_tssr.eps"')
g('set multiplot layout 2, 2')
g('set key left top')
g('set autoscale')
g('set xlabel "a (AU)"')
g('set ylabel "exp(-t/tdep)"')
g('set log y')
g('set xrange [0.5:]')
g.plot(item1, item2, item3, item5)
g('unset log')
g('set ylabel "a (AU)"')
g('set xlabel "Time (yr)"')
g('set xrange [:5e6]')
g('set xtics 1e6')
g.plot(item_t1, item_t2, item_t3, item_t5)
g('set key bottom right')
g('set ylabel "precession rate"')
g('set title "no gap"')
g.plot(item_gp1, item_gj1)
g('set title "with gap"')
g.plot(item_gp2, item_gj2)
g('unset multiplot')





f = open('t_assr','w')
for i in range(nstep):
  print >> f, data_t5[i][0], data_t5[i][1]
f.close()
