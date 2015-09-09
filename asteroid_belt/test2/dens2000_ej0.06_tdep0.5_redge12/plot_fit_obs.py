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

for line in open('ri_20_200'):
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

item = [(min(counter_i)-9.9e-5, 1e-10, 1e-10)]
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
g('set xlabel "R(km)"')
g('set log')
g('set ylabel "Frequency"')
g('set key right top')
g('set style fill border')
g('unset boxwidth')
g.set_range('yrange',(1e-4, 1.0))
g.set_range('xrange', (size_min, size_max))

item_model_i = Gnuplot.PlotItems.Data(
    data_i, title='initial', with_='yerrorlines lc rgb "red"')
item_model_ab = Gnuplot.PlotItems.Data(
    data_ab, title='final', with_='yerrorbars lc rgb "green"')
 
item_model_ab_smooth = Gnuplot.PlotItems.Data(
    data_ab, using='1:2 smooth bezier lw 2 lc rgb "green"', title='') 

item_obs = Gnuplot.PlotItems.File(
    'r_fobs', using='1:2', title='observation',
    with_='lp lc rgb "black" ')

g._add_to_queue([item_model_i, item_model_ab, item_model_ab_smooth, item_obs])

g.refresh()
g.close()
