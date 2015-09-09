#!/usr/bin/env python

from collections import Counter
import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh solid color eps'

def weight(x): return 1.0
def bin_id(x, binsize): return min(int(math.log10(x) / binsize), 22)
def bin(x, binsize): return 10 ** (binsize * (bin_id(x, binsize) + 0.5))

binsize = 0.3
def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in counter: print i, 10**(binsize*i), 10**(binsize*(i+1)), counter[i]

dir = [('dens2000_ej0.02_cd0.44', 'red'),
       ('dens2000_ej0.02_cd0.44_0.0001t', 'blue')]

g2 = Gnuplot.Gnuplot()

g2('set term postscript eps solid')
g2('set output "r_dratio.eps"')
g2('set xlabel "R(km)"')
g2('set ylabel "Depletion ratio"')
g2('set log')
g2.set_range('xrange', (2, 2000))

for loop, color in dir:

    rab = [float(i.split()[0]) for i in open(loop + '/rab')]
    ri = [float(i.split()[0]) for i in open('%s/ri' % loop)]
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

    def f(x):
        base = 1 - counter_ab[x] / counter_i[x]
        return x, base

    data = sorted(f(x) for x in counter_ab)
    item = Gnuplot.PlotItems.Data(
        data, title = 'ratio', with_='boxes lc rgb "%s"' % color)
    g2._add_to_queue([item])

g2.refresh()
g2.close()
