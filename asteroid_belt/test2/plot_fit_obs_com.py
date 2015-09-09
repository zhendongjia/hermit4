#!/usr/bin/env python

from collections import Counter
import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh solid color eps'

def weight(x): return x**(-3.5)
def bin_id(x, binsize): return min(int(math.log10(x) / binsize), 22)
def bin(x, binsize): return 10 ** (binsize * (bin_id(x, binsize) + 0.5))

binsize = 0.1

size_min = 20
size_max = 200

def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in counter: print i, 10**(binsize*i), 10**(binsize*(i+1)), counter[i]

data_dir = ['./dens2000_ej0.06_tdep0.5_redge12/','./dens2000_ej0.06_tdep0.5_redge12_add_saturn/','./dens2000_ej0.06_tdep0.5_redge12_densp2/']
data_ab = []
data_i = []

for j in data_dir:

    rab = []
    ri = []

    for line in open(j+'rab'):
        size = float(line.split()[0])
        if size >= size_min and size <= size_max:
            rab.append(size)

    for line in open(j+'ri'):
        size = float(line.split()[0])
        if size >= size_min and size <= size_max:
            ri.append(size)

    print 'rab'
    stat(rab)
    print 'ri'
    stat(ri)

    sum_weight_ab = sum(weight(i) for i in rab)
    sum_weight_i = sum(weight(i) for i in ri)

    counter_ab = Counter()
    counter_i = Counter()

    for i in rab: counter_ab[bin(i, binsize)] += weight(i)
    for i in ri: counter_i[bin(i, binsize)] += weight(i)

    num_ab = Counter()
    num_i = Counter()

    for i in rab: num_ab[bin_id(i, binsize)] += 1
    for i in ri: num_i[bin_id(i, binsize)] += 1

    item = [(min(counter_i) - 1e-6, 1e-6, 1e-6)]
    item.extend([(i, counter_ab[i]/sum_weight_ab,
                   counter_ab[i] / sum_weight_ab /
                   math.sqrt(num_ab[bin_id(i, binsize)]))
                  for i in counter_ab])
    data_ab.append(sorted(item))
    data_i.append(sorted([(i, counter_i[i]/sum_weight_i,
                  counter_i[i] / sum_weight_i / math.sqrt(num_i[bin_id(i, binsize)]))
                 for i in counter_i]))

g = Gnuplot.Gnuplot()

g('set term postscript eps solid')
g('set output "r_f_obs_com.eps"')
g('set xlabel "R(km)"')
g('set log')
g('set ylabel "Frequency"')
g('set key right top')
g('set style fill border')
g('unset boxwidth')
g.set_range('yrange',(1e-3, 1.0))
g.set_range('xrange', (size_min, size_max))


item_default = Gnuplot.PlotItems.Data(
    data_ab[0], title='default model', with_='yerrorbars lc rgb "green"')
item_default_smooth = Gnuplot.PlotItems.Data(
    data_ab[0], using='1:2', smooth='bezier',
    with_='l lw 2 lc rgb "green"', title = '')

item_saturn = Gnuplot.PlotItems.Data(
    data_ab[1], title='with saturn', with_='yerrorbars lc rgb "red"')
item_saturn_smooth = Gnuplot.PlotItems.Data(
    data_ab[1], using='1:2', smooth='bezier',
    with_='l lw 2 lc rgb "red"', title = '')

item_densp = Gnuplot.PlotItems.Data(
    data_ab[2], title='densp = 2', with_='yerrorbars lc rgb "yellow"')
item_densp_smooth = Gnuplot.PlotItems.Data(
    data_ab[2], using='1:2', smooth='bezier',
    with_='l lw 2 lc rgb "yellow"', title = '')

item_obs = Gnuplot.PlotItems.File(
    'r_fobs', using='1:2', title='observation',
    with_='lp lc rgb "black"')
g._add_to_queue([item_default, item_saturn, item_densp, item_obs, item_default_smooth, item_saturn_smooth, item_densp_smooth])

g.refresh()
g.close()
