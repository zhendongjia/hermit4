#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
import Gnuplot as gp
gp.GnuplotOpts.default_term = 'post enh solid color eps'
import Gnuplot.funcutils

kg = [1.5, 1.5, 1., 1.0]

redge = [12., 1., 8., 11.]

dens0 = [2000., 1700, 2500., 4500.]

tdep = [5e5, 1e6, 1e6, 1e6]

h = [1., 0.4, 0.4, 0.5]

s = [1., 1.0, 1., 1.0]

rlist = [i*0.5 for i in range(4, 7)]
tlist = [i*1e4 for i in range(1,1000)]

def f_grav(r, t1, t2, i):
  grav1 = 0.5/(1+kg[0])*(r/redge[0])**(1+kg[0]) + 144./(3+kg[0])*(r/redge[0])**(3+kg[0])
  grav1 = grav1*dens0[0]*r**(-kg[0])  
  grav1 = grav1 * math.exp(-t1/tdep[0])
  grav2 = 0.5/(1+kg[i])*(r/redge[i])**(1+kg[i]) + 0.5625/(3+kg[i])*(r/redge[i])**(3+kg[i])
  grav2 = grav2*dens0[i]*r**(-kg[i])
  grav2 = grav2 * math.exp(-t2/tdep[i]) 
  f_grav = grav2/grav1  
  return f_grav

def f_grav_doug(r, t1, t2, i):
  grav1 = 1.094 * 4 * math.pi * dens0[1] * r**(-kg[1]) 
  grav1 = grav1 * math.exp(-t1/tdep[1])
  grav2 = 0.5/(1+kg[i])*(r/redge[i])**(1+kg[i]) + 0.5625/(3+kg[i])*(r/redge[i])**(3+kg[i])
  grav2 = grav2*dens0[i]*r**(-kg[i])
  grav2 = grav2 * math.exp(-t2/tdep[i]) 
  f_grav = grav2/ grav1
  return f_grav




def rate(r, t, i):
  if i == 0: 
    grav = 0.5/(1+kg[i])*(r/redge[i])**(1+kg[i]) + 144./(3+kg[i])*(r/redge[i])**(3+kg[i])
  if i == 1:
    grav = 1.094 * 4 * math.pi * dens0[1] * r**(-kg[1]) 
  if i > 1:
    grav = 0.5/(1+kg[i])*(r/redge[i])**(1+kg[i]) + 0.5625/(3+kg[i])*(r/redge[i])**(3+kg[i])
  grav = grav*dens0[i]*r**(-kg[i])
  t0 = t - 1e4
  rate = grav * (math.exp(-t0/tdep[i]) - math.exp(-t/tdep[i]))
  rate = rate / 1e4
  return rate


def f_td(r, t1, t2, i):
  td1 = r**kg[0]/dens0[0]*0.2/6.28
  td1 = td1/math.exp(-t1/tdep[0])
  td2 = r**kg[i]/dens0[i]*0.2*h[i]**4/s[i]
  td2 = td2/math.exp(-t2/tdep[i])
  return td2/td1

def f_gd(r, t1, t2, i):
  td1 = r**kg[0]/dens0[i]/6.28
  td1 = td1/math.exp(-t1/tdep[0])
  td2 = r**kg[i]/dens0[i]*h[i]/s[i]
  td2 = td2/math.exp(-t2/tdep[i])
  return td2/td1

data1 = [(r, f_grav(r, 0., 0., 2)) for r in rlist]
data2 = [(r, f_grav(r, 0., 0., 3)) for r in rlist]

#data3 = [(r, f_grav(r, 10*tdep[0], 10*tdep[2], 2)) for r in rlist]
#data4 = [(r, f_grav(r, 10*tdep[0], 10*tdep[3], 3)) for r in rlist]

data5 = [(t, rate(3.5, t, 0)) for t in tlist]
data6 = [(t, rate(3.5, t, 2)) for t in tlist]
data7 = [(t, rate(3.5, t, 3)) for t in tlist]
data14 = [(t, rate(3.5, t, 1)) for t in tlist]

data8 = [(r, f_td(r, 0., 0., 2)) for r in rlist]
data9 = [(r, f_td(r, 0., 0., 3)) for r in rlist]

data10 = [(r, f_gd(r, 0., 0., 2)) for r in rlist]
data11 = [(r, f_gd(r, 0., 0., 3)) for r in rlist]

data12 = [(r, f_grav_doug(r, 0.0, 0.0, 2)) for r in rlist]
data13 = [(r, f_grav_doug(r, 0.0, 0.0, 3)) for r in rlist]

data15 = [(r, f_grav_doug(r, 0.0, 0.0, 2)/f_grav(r, 0.0, 0.0, 2)) for r in rlist]


g = gp.Gnuplot()
g('set term postscript eps solid')
g('set output "r_grav.eps"')
g('set multiplot layout 2, 2')
g('set key left top')
g('set autoscale')
g('set log y')

g('set style line 1 lt 1 lw 10 lc rgb "green"')
g('set style line 2 lt 1 lw 10 lc rgb "red"')
g('set style line 3 lt 2 lw 10 lc rgb "green"')
g('set style line 4 lt 2 lw 10 lc rgb "red"')
g('set style line 5 lt 1 lw 10 lc rgb "blue"')



item1 = gp.Data(data1, title = 'no saturn at T = 0 Myr', with_='lines ls 1')
item2 = gp.Data(data2, title = 'with saturn at T = 0 Myr', with_='lines ls 2')
item13 = gp.Data(data13, title = 'compare to Doug gravity', with_='lines ls 5')

#item3 = gp.Data(data3, title = 'no saturn at T = 10 Myr', with_='lines ls 3')
#item4 = gp.Data(data4, title = 'with saturn at T = 10 Myr', with_='lines ls 4')

item5 = gp.Data(data5, title = 'mmnm', with_='lines lc "black" lw 10')
item6 = gp.Data(data6, title = 'no saturn', with_='lines ls 1')
item7 = gp.Data(data7, title = 'with saturn', with_='lines ls 2')
item14 = gp.Data(data14, title = 'Doug', with_='lines ls 5')

item8 = gp.Data(data8, title = 'tidal damping no saturn', with_='line ls 1')
item9 = gp.Data(data9, title = 'tidal damping with saturn', with_='line ls 2')

item10 = gp.Data(data10, title = 'gas drag no saturn', with_='line ls 1')
item11 = gp.Data(data11, title = 'gas drag with saturn', with_='line ls 2')

item15 = gp.Data(data15, title = 'mmn/doug', with_='lines lc rgb "purple" lw 10')

g('set xlabel "a (AU)"')
g('set ylabel "force ratio"')
g.plot(item1, item2, item13, item15)

g('set ylabel "migration rate"')
g('set xlabel "T (yr)"')
g.plot(item6, item5, item7, item14)

g('set ylabel "tidal damping ratio"')
g.plot(item8, item9)

g('set ylabel "gas drag ratio"')
g.plot(item10, item11)

g('unset multiplot')
