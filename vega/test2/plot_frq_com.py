#!/usr/bin/env python

from collections import Counter
import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh solid color eps'

x_min = 10
x_max = 120

def weight(x): return 1.0
def bin_id(x, binsize): return math.floor((x - x_min)/binsize)
def bin(x, binsize): return (x_min + binsize * bin_id(x, binsize))

binsize = 8
def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in sorted(counter): print i, x_min+(binsize*i), x_min+(binsize*(i+1)), counter[i]

dir = ['./aj60_ej0.2', './aj60_ej0.2_ng']
data = []
data_m = []
for loop in dir:
    ap_i = [float(i.split()[0]) for i in open(loop+'/ap_r_i')]
    ap_f = [float(i.split()[0]) for i in open(loop+'/ap_r_f')]
    r_i = [float(i.split()[1]) for i in open(loop+'/ap_r_i')]
    r_f = [float(i.split()[1]) for i in open(loop+'/ap_r_f')]
    print 'ap_i'
    stat(ap_i)
    print 'ap_f'
    stat(ap_f)
    sum_weight_i = sum(weight(i) for i in ap_i)
    sum_weight_f = sum(weight(i) for i in ap_f)
    n_i = len(r_i)
    n_f = len(r_f)
    counter_i = Counter()
    counter_f = Counter()
    counter_m_i = Counter()
    counter_m_f = Counter()
    for i in ap_i: counter_i[bin(i, binsize)] += weight(i)
    for i in ap_f: counter_f[bin(i, binsize)] += weight(i)
    for idx in range(n_i): counter_m_i[bin(ap_i[idx], binsize)] += r_i[idx]**3
    for idx in range(n_f): counter_m_f[bin(ap_f[idx], binsize)] += r_f[idx]**3

    def h(x):
        if counter_i[x] == 0 or counter_f[x] == 0:
            base = 0.0
            error = 0.0
        else:
            base = float(counter_f[x])/counter_i[x]
            if counter_f[x] == 0:
                erro = 0.0
            else:
                error = 1 / math.sqrt(counter_f[x]) + 1/ math.sqrt(counter_i[x])
	x = x + binsize/2.
        return x, base, base * (1 - error), base * (1 + error)

    def h_m(x):
        if counter_m_i[x] == 0 or counter_m_f[x] == 0:
            base = 0.0
            error = 0.0
        else:
            base = float(counter_m_f[x])/counter_m_i[x]
            if counter_m_f[x] == 0:
                erro = 0.0
            else:
                error = 1 / math.sqrt(counter_f[x]) + 1/ math.sqrt(counter_i[x])
        x = x + binsize/2.
	return x, base, base * (1 - error), base * (1 + error)


    xlist = range(x_min, x_max, binsize) 

    data.append(sorted(h(x) for x in xlist))
    data_m.append(sorted(h_m(x) for x in xlist))




##########################################################################################
item0 = Gnuplot.Data(data[0], title = '', with_='boxerrorbars lc rgb "blue"')
item_m0 = Gnuplot.Data(data_m[0], title = '', with_='boxerrorbars lc rgb "blue"')
item_r = Gnuplot.PlotItems.File('./aj60_ej0.2/ap_r_f',title = '', using=' 1:2 axes x1y2')

g = Gnuplot.Gnuplot()
g('set term postscript eps solid')
g('set output "ap_nd_md.eps"')
g('set title "aj = 60, ej = 0.2"')
g('set key center top')
g('set macros')
g('TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.55"')
g('BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"')
g('LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.90"')
g('set multiplot layout 2,1')
g('set ylabel "Residual number fraction"')
g.set_range('xrange', (x_min, x_max))
g.set_range('yrange',(0.0,2.5))
g('@TMARGIN;@LMARGIN')
g('unset xtics')
g('set boxwidth -2')
g.plot(item0)
g('@BMARGIN;@LMARGIN')
g('set xtics')
g('set xlabel "ap (AU)"')
g('set ylabel "Residual mass fraction"')
g('unset title')
g('set y2tics 100 nomirror tc lt 1')
g('set log y2')
g('set y2label "Size (km)" tc lt 1')
g.set_range('y2range',(1e-6, 1e3))
g.plot(item_m0)
g.plot(item_r)
g('unset multiplot')
g.close()

########################################################################
item_m1 = Gnuplot.Data(data_m[1], title = 'no gas ', with_='boxerrorbars lc rgb "blue"')
item_m0 = Gnuplot.Data(data_m[0], title = 'with gas', with_='boxerrorbars lc rgb "blue"')

g = Gnuplot.Gnuplot()
g('set term postscript eps solid')
g('set output "ap_md_gas.eps"')
g('set title "aj = 60, ej = 0.2, mj = 1."')
g('set key center top')
g('set macros')
g('TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.55"')
g('BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"')
g('LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.95"')
g('set multiplot layout 2,1')
g('set ylabel "Residual mass fraction"')
g.set_range('xrange', (x_min, x_max))
g.set_range('yrange',(0.0,1.3))
g('@TMARGIN;@LMARGIN')
g('unset xtics')
g('set boxwidth -2')
g.plot(item_m0)
g('@BMARGIN;@LMARGIN')
g('set xtics')
g('set xlabel "ap (AU)"')
g('unset title')
g.plot(item_m1)
g('unset multiplot')
g.close()
