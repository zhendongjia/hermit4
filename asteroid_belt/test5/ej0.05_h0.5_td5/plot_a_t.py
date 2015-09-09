#!/usr/bin/env python

import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh color eps'



data_dir = range(0,10000)

color_min=1e100
color_max=0

g = Gnuplot.Gnuplot()
g('set term post enh color eps')
g('set output "a_t_r.eps"')
g('set xlabel "T (yr)"')
g('set ylabel "Location (AU)"')
g.set_range('xrange', (1, 1e7))
g.set_range('yrange', (1.5, 3.5))
g.set_range('cbrange', (10, 10**3.0))
g('set pm3d')
g('set logscale zcb')
#g('set logscale x')

items = [ ]
for dirname in data_dir:
    if not os.path.exists('%d'%(dirname)): continue
    print dirname
    os.chdir('%s' % dirname)

    if not os.path.exists('radius'):
        os.chdir('..')
        continue
    if os.stat('radius')[6]==0:
        os.chdir('..')
        continue

    if not os.path.exists('a_t'):
        os.chdir('..')
        continue

    f_color=open('radius')
    radius_open = float(open('radius').readlines()[0].split()[0])
    af_open = float(open('a_t').readlines()[-1].split()[0])
    if radius_open > 1000.:
        os.chdir('..')
        continue
    try:
        color = float(f_color.readlines()[0].split()[0])
    except:
        os.chdir('..')
        continue
    if color < color_min: color_min = color
    if color > color_max: color_max = color
    f_color.close()

    f_data=open('a_t')
    print os.getcwd()
    data = [line.split() for line in f_data.readlines()]
    data = [(float(i[3]),float(i[0])) for i in data if len(i)==4]
    f_data.close()



    #item = Gnuplot.PlotItems.Data(data, with_='lines color palette cb %f' % color)
    item = Gnuplot.PlotItems.Data(data, with_='l lt 1 lc palette frac %f' % ((math.log10(color)-1.0)/2.8) )
    g._add_to_queue([item])
    os.chdir('..')


    item0 = Gnuplot.PlotItems.File('~/Hermit4/asteroid_belt/analytic/gj_png', title = '', with_='l lt 1 lc rgb "green" lw 5' )

g._add_to_queue([item0])

g.refresh()
g.close()
