#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
import Gnuplot as gp
gp.GnuplotOpts.default_term = 'post enh color eps'
import Gnuplot.funcutils
from scipy.integrate import quad


def integrand(x, s, n, alpha):
  r = math.cos(n*x) / (1 + alpha**2 - 2*alpha*math.cos(x))**s
  r = r/math.pi
  return r


def bsn(s, n, alpha):
  b = quad(integrand, 0 , 2*math.pi, args = (s, n, alpha))
  bf = b[0]
  return bf


def Gjd(dens0, k, zk, aj):
  gjd = - zk * math.pi * dens0 * aj**(-k)
  nj = aj**(-1.5) 
  gjd = gjd / (nj * aj)
  return gjd


def Gjd_wg(dens_in, dens_out, k, zk, aj, a_in, a_out):
  nj = aj**(-1.5)
  gjd = 2 * math.pi * aj**(-kg)/ (nj * aj)
  A1 = 1./4
  A2 = 9./64
  A3 = 25./256
  gjd1 = gjd * 3 * A1 * (dens_in *(a_in/aj)**(4-k)/(4-k) + dens_out * (aj/a_out)**(1+k)/(1+k))
  gjd2 = gjd * 10 * A2 * (dens_in * (a_in/aj)**(6-k)/(6-k) + dens_out * (aj/a_out)**(3+k)/(3+k))
  gjd3 = gjd * 21 * A3 * (dens_in * (a_in/aj)**(8-k)/(8-k) + dens_out * (aj/a_out)**(5+k)/(5+k))
  gjd = gjd1 + gjd2 + gjd3
  return gjd


def Gjd_sg(dens_in, dens_out, k, zk, aj, ej):
  nj = aj**(-1.5)
  gjd = math.pi / (nj * aj)
  A1 = 1./4
  A2 = 9./64
  A3 = 25./256
  a_in = aj * (1-1.5*ej)
  a_out = aj * (1+1.5*ej)
  gjd1 = gjd * 3 * A1 * (dens_in *(a_in/aj)**(4-k)/(4-k) + dens_out * (aj/a_out)**(1+k)/(1+k))
  gjd2 = gjd * 10 * A2 * (dens_in * (a_in/aj)**(6-k)/(6-k) + dens_out * (aj/a_out)**(3+k)/(3+k))
  gjd3 = gjd * 21 * A3 * (dens_in * (a_in/aj)**(8-k)/(8-k) + dens_out * (aj/a_out)**(5+k)/(5+k))
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
m_s = 0.000285837
a_j = 5.202545
a_s = 9.554841

dens0 = 1700 * 1.125e-7
k = 1.5
zk = 1.094
tdep = 1e6
a_in = 4.5
a_out = 11.0

g_js = Gij(a_j, a_s, m_s)
g_sj = Gij(a_s, a_j, m_j)


#---------------------------------------------------------------------------------
for i in range(100):
  while True:
    ap = random.uniform(5,40)
    ap = 0.1*ap
    t = random.uniform(0,5)
    t = 10**(t-4)

    g_jd = Gjd(dens0, k, z_k, a_j)
    g_sd = Gjd(dens0, k, z_k, a_s)
    g_pd = Gjd(dens0, k, z_k, ap)

    g_pj = Gij(ap, a_j, m_j)
    g_ps = Gij(ap, a_s, m_s)

    g_j = g_js + g_jd*t
    g_s = g_sj + g_sd*t

    g_5 = 0.5*(g_j + g_s + ((g_j - g_s)**2 + 4*g_js**2)**0.5)
    g_6 = 0.5*(g_j + g_s - ((g_j - g_s)**2 + 4*g_js**2)**0.5)

    g_p = g_ps + g_pj + g_pd*t

    if abs((g_p - g_5) / g_5) < 1e-3:
      data_v5
      break
    if abs((g_p - g_6) / g_6) < 1e-3:
      print >> f1, ap, t
      break




for i in range(nstep): 
  g_pd[i] = Gjd(dens0, k, zk, ap[i])
  g_pj[i] = Gij(ap[i], a_j, m_j)
  g_ps[i] = Gij(ap[i], a_s, m_s)
  
  exp_t[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd - g_pd[i])
  exp_t_p4[i] = 6e-3 * ap[i]**4

  t[i] = math.log(abs(1./exp_t[i]))*tdep
  t_p4[i] = math.log(abs(1./exp_t_p4[i]))*tdep


  if ap[i] < a_j*(1-0.2):
    exp_t_2d[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_2d - g_pd[i])
  else:
    exp_t_2d[i] = 0.0

  if ap[i] < a_j*(1-0.3):
    exp_t_3d[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_3d - g_pd[i])
  else:
    exp_t_3d[i] = 0.0

  if ap[i] < a_j*(1-0.4):
    exp_t_4d[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_4d - g_pd[i])
  else:
    exp_t_4d[i] = 0.0
    
  if ap[i] < a_j*(1-0.5):
    exp_t_5d[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_5d - g_pd[i])
  else:
    exp_t_5d[i] = 0.0


  

  exp_t_ns[i] = g_pj[i]/(g_jd - g_pd[i])
  t_ns[i] = math.log(abs(1/exp_t_ns[i]))*tdep

  exp_t_wg[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_wg - g_pd[i])
  t_wg[i] = math.log(abs(1./exp_t_wg[i]))*tdep

  t_2d[i] = math.log(abs(1./exp_t_2d[i]))*tdep

  exp_t_sg[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_sg - g_pd[i])
  t_sg[i] = math.log(abs(1/exp_t_sg[i]))*tdep


  exp_t_sat[i] = (g_pj[i] + g_ps[i] - g_sj)/(g_sd - g_pd[i])

  exp_t_wg_sat[i] = (g_pj[i] + g_ps[i] - g_sj)/(g_sd_wg - g_pd[i])
  t_wg_sat[i] = math.log(abs(1./exp_t_wg[i]))*tdep

  exp_t_sg_sat[i] = (g_pj[i] + g_ps[i] - g_sj)/(g_sd_sg - g_pd[i])
  t_sg_sat[i] = math.log(abs(1/exp_t_sg[i]))*tdep


  exp_t_wg_cons[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_wg_cons - g_pd[i])


  print g_pj[i], g_ps[i], g_js, g_sj, g_pd[i], g_jd, g_sd, g_jd_wg, g_sd_wg

  
  exp_t_ni[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_ni - g_pd[i])
  t_ni[i] = math.log(abs(1/exp_t_ni[i]))*tdep

  g_pd_li[i] = Gjd(dens0*0.2, k, zk, ap[i])
  exp_t_li[i] = (g_pj[i] + g_ps[i] - g_js)/(g_jd_li - g_pd_li[i])
  t_li[i] = math.log(abs(1/exp_t_li[i]))*tdep

#print exp_t_wg[-1], exp_t[-1], ap[-1]  
#print math.log10(1./exp_t_wg[-1]), ap [-1]  
  
  
data_expt = [(ap[i], exp_t[i]) for i in range(nstep)]
data_t = [(t[i], ap[i]) for i in range(nstep)]
data_expt_p4 = [(ap[i], exp_t_p4[i]) for i in range(nstep)]
data_t_p4 = [(t_p4[i], ap[i]) for i in range(nstep)]
data_expt_sat = [(ap[i], exp_t_sat[i]) for i in range(nstep)]

data_expt_2d = [(ap[i], exp_t_2d[i]) for i in range(nstep)]
data_expt_3d = [(ap[i], exp_t_3d[i]) for i in range(nstep)]
data_expt_4d = [(ap[i], exp_t_4d[i]) for i in range(nstep)]
data_expt_5d = [(ap[i], exp_t_5d[i]) for i in range(nstep)]

data_expt_ns = [(ap[i], exp_t_ns[i]) for i in range(nstep)]
data_t_ns = [(t_ns[i], ap[i]) for i in range(nstep)]

data_expt_wg = [(ap[i], exp_t_wg[i]) for i in range(nstep)]
data_t_wg = [(t_wg[i], ap[i]) for i in range(nstep)]

data_expt_sg = [(ap[i], exp_t_sg[i]) for i in range(nstep)]
data_t_sg = [(t_sg[i], ap[i]) for i in range(nstep)]

data_expt_wg_sat = [(ap[i], exp_t_wg_sat[i]) for i in range(nstep)]
data_t_wg_sat = [(t_wg_sat[i], ap[i]) for i in range(nstep)]

data_expt_sg_sat = [(ap[i], exp_t_sg_sat[i]) for i in range(nstep)]
data_t_sg_sat = [(t_sg_sat[i], ap[i]) for i in range(nstep)]

data_expt_wg_cons = [(ap[i], exp_t_wg_cons[i]) for i in range(nstep)]



data_t_2d = [(t_2d[i], ap[i]) for i in range(nstep)]

data_expt_ni = [(ap[i], exp_t_ni[i]) for i in range(nstep)]
data_t_ni = [(t_ni[i], ap[i]) for i in range(nstep)]

data_t_10tdep = [(t_wg[i]*10, ap[i]) for i in range(nstep)]

#data_gjd = [(i*tdep*0.1, g_jd_wg*math.exp(-i*0.1)) for i in range(100)]
data_gjd_ni = [(i*tdep*0.1, g_jd_ni*math.exp(-i*0.1)) for i in range(100)]
data_gjd_ain = [(ain[i], g_jd_ain[i]) for i in range(50)]

data_t_li = [(t_li[i], ap[i]) for i in range(nstep)]

data_gjs = [(t[i], g_js) for i in range(nstep)]
data_gjd = [(t[i], g_jd_wg*math.exp(-t[i]/tdep)) for i in range(nstep)]
data_gpd = [(t[i], g_pd[i]*math.exp(-t[i]/tdep)) for i in range(nstep)]
data_gpj = [(t[i], g_pj[i]) for i in range(nstep)]
data_gps = [(t[i], g_ps[i]) for i in range(nstep)]

#print g_js, g_ps


item_expt_p4 = gp.Data(data_expt_p4, title = '~ a^4', with_='l lw 10')
                

item_expt_2d = gp.Data(data_expt_2d, title = '\delta_{d} = 0.2', with_='l lw 10')
item_expt_3d = gp.Data(data_expt_3d, title = '\delta_{d} = 0.3', with_='l lw 10')
item_expt_4d = gp.Data(data_expt_4d, title = '\delta_{d} = 0.4', with_='l lw 10')
item_expt_5d = gp.Data(data_expt_5d, title = '\delta_{d} = 0.5', with_='l lw 10')
item_t_2d = gp.Data(data_t_2d, title = '\delta_{d} = 0.2', with_='l lw 10')

item_expt_ns = gp.Data(data_expt_ns, title = 'no saturn', with_='l lw 10')




item_expt = gp.Data(data_expt, title = 'no gap', with_='l lt 2 lc 0 lw 5')
item_t = gp.Data(data_t, title = 'no gap', with_='l lt 2 lc 0 lw 5')
item_expt_sat = gp.Data(data_expt_sat, title = '', with_='l lt 2 lc 9 lw 5')
item_t_sat = gp.Data(data_t, title = '', with_='l lt 2 lc 9 lw 5')

item_expt_wg = gp.Data(data_expt_wg, title = 'large gap', with_='l lt 1 lc 0 lw 5')
item_t_wg = gp.Data(data_t_wg, title = 'large gap', with_='l lt 1 lc 0 lw 5')

item_expt_sg = gp.Data(data_expt_sg, title = 'small gap', with_='l lt 0 lc 0 lw 5')
item_t_sg = gp.Data(data_t_sg, title = 'small gap', with_='l lt 0 lc 0 lw 5')

item_expt_wg_sat = gp.Data(data_expt_wg_sat, title = '', with_='l lt 1 lc 9 lw 5')
item_t_wg_sat = gp.Data(data_t_wg, title = '', with_='l lt 1 lc 9 lw 5')

item_expt_sg_sat = gp.Data(data_expt_sg_sat, title = '', with_='l lt 0 lc 9 lw 5')
item_t_sg_sat = gp.Data(data_t_sg_sat, title = '', with_='l lt 0 lc 9 lw 5')

item_expt_wg_cons = gp.Data(data_expt_wg_cons, title = 'cons gap', with_='l lt 3 lc 0 lw 5')



item_t_10tdep = gp.Data(data_t_10tdep, title = 'with gap, tdep = 10Myr', with_='l lw 10')

item_expt_ni = gp.Data(data_expt_ni, title = 'with gap, no inner disk', with_='l lw 10')
item_t_ni = gp.Data(data_t_ni, title = 'with gap, no inner disk', with_='l lw 10')

item_gjd = gp.Data(data_gjd, title = 'with inner disk', with_='line lw 10')
item_gjd_ni = gp.Data(data_gjd_ni, title = 'no inner disk', with_='line lw 10')
item_gjd_ain = gp.Data(data_gjd_ain, title = '', with_='line lw 10')

item_t_li = gp.Data(data_t_li, title = 'reduction 5 time in inner disk', with_='l lw 10')

item_gjs = gp.Data(data_gjs, title = 'saturn', with_='l lw 10')
item_gpd = gp.Data(data_gpd, title = 'disk', with_='l lw 10')
item_gpj = gp.Data(data_gpj, title = 'jupiter', with_='l lw 10')
item_gps = gp.Data(data_gps, title = 'saturn', with_='l lw 10')




g = gp.Gnuplot()
g('set term postscript eps  color')
g('set output "ap_t.eps"')
g('set key left top')
g('set autoscale')
g('set xlabel "a (AU)"')
g('set ylabel "exp(-t/Tdep)"')
g('set yrange [1e-4:5.0]')
g('set log y')
#g('set grid')
g('set yrange [:2]')
g.plot(item_expt, item_expt_sat, item_expt_wg, item_expt_sg,  item_expt_wg_sat, item_expt_sg_sat)
g('unset log')
g('set output')


g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "t_ap.eps"')
g('set key left top')
g('set autoscale')
g('set ylabel "a (AU)"')
g('set xlabel "Time (yr)"')
#g('set log x')
g('set xrange [0:]')
g('set yrange [0.85:3.2]')
g.plot(item_t, item_t_wg, item_t_sg, item_t_sat, item_t_wg_sat, item_t_sg_sat )
g('set output')



g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "ap_tssr_com.eps"')
g('set multiplot layout 2, 2')
g('set key left top')
g('set autoscale')
g('set grid')
g('set xlabel "a (AU)"')
g('set ylabel "exp(-t/tdep)"')
g('set log y')
g('set xrange [0.5:]')
g('set yrange [0.0001:5]')
g.plot(item_expt, item_expt_sat)
g.plot(item_expt, item_expt_wg, item_expt_2d, item_expt_3d, item_expt_4d, item_expt_5d)
g.plot(item_expt_sg, item_expt_wg)
g('set ylabel "a (AU)"')
g('set xlabel "Time (yr)"')
g('unset log')
g('set xtics 2e6')
g('set xrange [:7e6]')
g('set yrange [1:4]')
g.plot(item_t_wg, item_t_sg, item_t)
g('unset multiplot')




g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "g_t.eps"')
g('set multiplot layout 1, 2')
g('set key left top')
g('set autoscale')
g('set xlabel "Time (yr)"')
g('set ylabel "precession rate"')
#g('set log y')
#g('set xrange [2.0:10.0]')
g('set title "Jupiter"')
g.plot(item_gjd, item_gjs)
g('set title "Planetesimal"')
g.plot(item_gpd, item_gpj, item_gps)
g('unset multiplot')



g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "ap_grav.eps"')
#g('set multiplot layout 2, 1')
g('set key left top')
g('set autoscale')
g('set xlabel "a (AU)"')
g('set ylabel "grav"')
#g('set log y')
#g('set xrange [2.0:10.0]')
#g.plot(item_grav_ratio)
g('set output')


f = open('t_assr','w')
for i in range(nstep):
  print >> f, data_t[i][0], data_t[i][1]
f.close()


f = open('t_assr_wg','w')
for i in range(nstep):
  print >> f, data_t_wg[i][0], data_t_wg[i][1]
f.close()

f = open('t_assr_wg_2dep','w')
for i in range(nstep):
  print >> f, data_t_wg[i][0]*2, data_t_wg[i][1]
f.close()


f = open('t_assr_p4', 'w')
for i in range(nstep):
  print >> f, data_t_p4[i][0], data_t_p4[i][1]
f.close()


f = open('t_assr_10tdep','w')
for i in range(nstep):
  print >> f, data_t_10tdep[i][0], data_t_10tdep[i][1]
f.close()

f = open('t_assr_ni','w')
for i in range(nstep):
  print >> f, data_t_ni[i][0], data_t_ni[i][1]
f.close()
