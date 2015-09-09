#!/usr/bin/env python

from collections import Counter
import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh solid color eps'

def weight(x): return x**(-3.5)
def bin_id(x, binsize): return int(math.log10(x) / binsize)
def bin(x, binsize): return 10 ** (binsize * (bin_id(x, binsize) + 0.5))

binsize = 0.1

size_min = 10.
size_max = 1000.

def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in counter: print i, 10**(binsize*i), 10**(binsize*(i+1)), counter[i]

#rab = [float(i.split()[0]) for i in open('rab')]
#ri = [float(i.split()[0]) for i in open('ri')]

data_ab = []
data_i = []
rab = []
ri = []

for line in open('rab'):
    size = float(line.split()[0])
    if size >= size_min and size <= size_max:
        rab.append(size)

#for line in open('ri_20_200'):
for line in open('ri'):
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

item = []
#item = [(min(counter_i)-9.9e-5, 1e-10, 1e-10)]
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
g('set output "r_f_obs.eps"')
g('set xlabel "size (km)"')
g('set log')
g('set ylabel "Frequency"')
g('set key right top')
g('set style fill border')
g('unset boxwidth')
g.set_range('yrange',(1e-7, 1.0))
g.set_range('xrange', (size_min, 900))

item_model_i = Gnuplot.PlotItems.Data(
    data_i, title='', with_='yerrorlines lw 5 ps 1.5 lc rgb "black"')
item_model_ab = Gnuplot.PlotItems.Data(
    data_ab, title='', with_='yerrorbars pt 7 ps 1.5 lc rgb "black"')
 
item_model_ab_smooth = Gnuplot.PlotItems.Data(
    data_ab, using='1:2 smooth bezier lw 5 lc rgb "black"', title='') 

item_obs = Gnuplot.PlotItems.File(
    'r_fobs', using='1:2', title='observation',
    with_='lp lc rgb "black" ')

item_d_obs = Gnuplot.PlotItems.File(
    'd_fobs', using='1:2', title='observation',
    with_='lp lw 5 lc rgb "grey" ')



g._add_to_queue([item_model_i, item_model_ab, item_model_ab_smooth, item_d_obs])
g('set label "T = 0 Myr" at 20, 0.01')
g('set label "T = 10 Myr" at 200, 0.02')

g.refresh()
g.close()
