#!/usr/bin/env python
import os
import math
import Gnuplot as gp
import Gnuplot.funcutils

file_name = 'asteroids_data'

d = [float(line.split()[0]) for line in open(file_name).readlines()]
n = [float(line.split()[1]) for line in open(file_name).readlines()]


count_i = 0
n_tot = 0
for i in range(len(n)):
    if d[i] > 15.:
        count_i = count_i + 1
        n_tot = n_tot + n[i]


r = [None] * count_i
f_r_nor = [None] * count_i
f_r = [None] * count_i
f_r_cul = [None] * count_i

n_cul = 0
n_max = n[-1]
for i in range(count_i):
    r[i] = d[i]/2
    f_r_nor[i] = n[i]/n_tot
    f_r[i] = n[i]/n_max
    n_cul = n[i] + n_cul
    f_r_cul[i] = n_cul/n_tot

data_rf_nor = [(r[i], f_r_nor[i]) for i in range(len(r))]
data_rf = [(r[i], f_r[i]) for i in range(len(r))]
data_rf_cul = [(r[i], f_r_cul[i]) for i in range(len(r))]

item_rf_nor = gp.Data(data_rf_nor, title = 'normalized', with_='linespoints lt 1 lw 2 pt 7 ps 1.5')
item_rf = gp.Data(data_rf, title = '', with_='linespoints lt 1 lw 2 pt 7 ps 1.5')
item_rf_cul = gp.Data(data_rf_cul, title = 'cul', with_='linespoints lt 1 lw 2 pt 7 ps 1.5')
    
g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "r_fobs.eps"')
g('set xlabel "radius (km)"')
g('set ylabel "frequency"')
g('set log')
g('set yrange[0.0002:]')
g('set xrange[7:]')
g.plot(item_rf_nor)
g('unset output')

f = open('r_fobs','w')
for i in range(len(r)):
    print >> f, r[i], f_r_nor[i]
f.close()

f = open('d_fobs','w')
for i in range(len(r)):
    print >> f, r[i]*2, f_r_nor[i]
f.close()
