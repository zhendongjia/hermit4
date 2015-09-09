#!/usr/bin/env python

from collections import Counter
import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh solid color eps'

def weight(x): return 1.0
def bin_id(x, binsize): return min(int(math.log10(x) / binsize), 22)
def bin(x, binsize): return 10 ** (binsize * (bin_id(x, binsize) + 0.5))

binsize = 0.1
def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in counter: print i, 10**(binsize*i), 10**(binsize*(i+1)), counter[i]

data_dir = ['./dens2000_ej0.06_tdep0.5_redge12/', './dens2000_ej0.06_tdep0.5_redge10/','./dens2000_ej0.04_tdep0.5_redge12/', './dens2000_ej0.06_tdep1.0_redge12/','./dens2000_ej0.06_tdep0.5_redge12_add_saturn/','./dens2000_ej0.06_tdep0.5_redge12_densp2/']
data = []

for j in data_dir:

    rab = [float(i.split()[0]) for i in open(j+'rab')]
    ri = [float(i.split()[0]) for i in open(j+'ri')]

    stat(rab)
    print 'rab'
    stat(ri)
    print 'ri'

    sum_weight_ab = sum(weight(i) for i in rab)
    sum_weight_i = sum(weight(i) for i in ri)

    counter_ab = Counter()
    counter_i = Counter()

    for i in rab: counter_ab[bin(i, binsize)] += weight(i)
    for i in ri: counter_i[bin(i, binsize)] += weight(i)

    def f(x):
    	base = counter_ab[x] / counter_i[x]
    	error = 1 / math.sqrt(counter_ab[x]) + 1 / math.sqrt(counter_i[x])
    	return x, base, base * (1 - error), base * (1 + error)

    data.append(sorted(f(x) for x in counter_ab))
########################################################
g = Gnuplot.Gnuplot()
g('set boxwidth -2')
g('set term postscript eps solid')
g('set output "r_df_redge.eps"')
g('set xlabel "R(km)"')
g('set ylabel "Residual fraction"')
g('set log x')
g('set key right top')
g.set_range('xrange', (10, 3000))
g.set_range('yrange', (0.0, 0.05))
g('set style fill pattern border')

item_default = Gnuplot.PlotItems.Data(data[0], title = 'Default model', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')
item_redge10 = Gnuplot.PlotItems.Data(data[1], title = 'Redge = 10AU', with_='boxerrorbars lc rgb "blue"')

g._add_to_queue([item_default, item_redge10])
g.refresh()
g.close()

########################################################
g = Gnuplot.Gnuplot()
g('set boxwidth -2')
g('set term postscript eps solid')
g('set output "r_df_ej.eps"')
g('set multiplot layout 2, 1')
g('set xlabel "R(km)"')
g('set ylabel "Residual fraction"')
g('set log x')
g('set key right top')
g.set_range('xrange', (10, 3000))
g.set_range('yrange', (0.0, 0.07))
g('set grid ytics lc rgb "black" lw 1 lt 0')
g('set grid xtics lc rgb "black" lw 1 lt 0')

item_default = Gnuplot.PlotItems.Data(data[0], title = 'ej = 0.06', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')
item_ej004 = Gnuplot.PlotItems.Data(data[2], title = 'ej = 0.04', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')

g.plot(item_default)
g.plot(item_ej004)
#g._add_to_queue((item_default, item_ej004))
#g.replot()
#g.refresh()
g.close()


##########################################################
g2 = Gnuplot.Gnuplot()
g2('set boxwidth -2')
g2('set term postscript eps solid')
g2('set output "r_df_tdep.eps"')
g2('set multiplot layout 2, 1')
g2('set xlabel "R(km)"')
g2('set ylabel "Residual fraction"')
g2('set log x')
g2('set key right top')
g2.set_range('xrange', (10, 3000))
g2.set_range('yrange', (0.0, 0.03))
g2('set grid ytics lc rgb "black" lw 1 lt 0')
g2('set grid xtics lc rgb "black" lw 1 lt 0')

item_tdep = Gnuplot.PlotItems.Data(data[3], title = 'Tdep = 1 Myr', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')
item_default = Gnuplot.PlotItems.Data(data[0], title = 'Tdep = 0.5 Myr', with_='boxerrorbars lc rgb "black" fillstyle pattern 2 ')

g2.plot(item_default)
g2.plot(item_tdep)
#g2._add_to_queue([item_tdep, item_default])
#g2.refresh()
g2.close()

#########################################################
g2 = Gnuplot.Gnuplot()
g2('set boxwidth -2')
g2('set term postscript eps solid')
g2('set output "r_df_saturn.eps"')
g2('set multiplot layout 2, 1')
g2('set xlabel "R(km)"')
g2('set ylabel "Residual fraction"')
g2('set log x')
g2('set key right top')
g2.set_range('xrange', (10, 3000))
g2.set_range('yrange', (0.0, 0.03))
g2('set grid ytics lc rgb "black" lw 1 lt 0')
g2('set grid xtics lc rgb "black" lw 1 lt 0')

item_saturn = Gnuplot.PlotItems.Data(data[4], title = 'with saturn', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')
item_default = Gnuplot.PlotItems.Data(data[0], title = 'no saturn', with_='boxerrorbars lc rgb "black" fillstyle pattern 2 ')

g2.plot(item_default)
g2.plot(item_saturn)
#g2._add_to_queue([item_tdep, item_default])
#g2.refresh()
g2.close()
#####################################################################
g2 = Gnuplot.Gnuplot()
g2('set boxwidth -2')
g2('set term postscript eps solid')
g2('set output "r_df_densp.eps"')
g2('set multiplot layout 2, 1')
g2('set xlabel "R(km)"')
g2('set ylabel "Residual fraction"')
g2('set log x')
g2('set key right top')
g2.set_range('xrange', (10, 1000))
g2.set_range('yrange', (0.0, 0.03))
g2('set grid ytics lc rgb "black" lw 1 lt 0')
g2('set grid xtics lc rgb "black" lw 1 lt 0')

item_densp = Gnuplot.PlotItems.Data(data[5], title = 'densp = 2', with_='boxerrorbars lc rgb "black" fillstyle pattern 2')
item_default = Gnuplot.PlotItems.Data(data[0], title = 'default model', with_='boxerrorbars lc rgb "black" fillstyle pattern 2 ')

g2.plot(item_default)
g2.plot(item_densp)
#g2._add_to_queue([item_tdep, item_default])
#g2.refresh()
g2.close()




