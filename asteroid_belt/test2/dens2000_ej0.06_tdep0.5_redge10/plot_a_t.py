#!/usr/bin/env python

import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh color eps'

data_dir = range(0,8000)

color_min=1e100
color_max=0

g = Gnuplot.Gnuplot()
g('set term post enh color eps')
g('set output "a_t_r.eps"')
g.set_range('xrange', (1, 5e6))
g.set_range('yrange', (1.0, 4.0))
g.set_range('cbrange', (10, 3000))
g('set pm3d')
g('set logscale zcb')
#g('set logscale x')

for dirname in data_dir:
    os.chdir('%s' % dirname)

    f_color=open('radius')
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
    data = [(float(i[2]),float(i[0])) for i in data if len(i)==3]
    f_data.close()

    #item = Gnuplot.PlotItems.Data(data, with_='lines color palette cb %f' % color)
    item = Gnuplot.PlotItems.Data(data, with_='l lt 1 lc palette frac %f' % ((math.log10(color)-1.0)/2.8) )
    g._add_to_queue([item])
    os.chdir('..')

g.refresh()
g.close()
