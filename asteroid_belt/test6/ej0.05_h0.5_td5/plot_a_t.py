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
g('set xlabel "t (xT_{dep} yr)"')
g('set ylabel "Location (AU)"')
g.set_range('xrange', (0, 10))
g.set_range('yrange', (1.5, 3.7))
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
    item = Gnuplot.PlotItems.Data(data, with_='l lt 1 lw 0.5 lc palette frac %f' % ((math.log10(color)-1.0)/2.8) )
    g._add_to_queue([item])
    os.chdir('..')


f5 = open('v5_wg_we_ord')
data_v5 = [(float(line.split()[2]),float(line.split()[0])) for line in f5.readlines()]
item5 = Gnuplot.PlotItems.Data(data_v5, title = '', with_='l lw 5 lt 1 lc 2 ' )

f6 = open('v6_wg_we_ord')
data_v6 = [(float(line.split()[2]),float(line.split()[0])) for line in f6.readlines()]
item6 = Gnuplot.PlotItems.Data(data_v6, title = '', with_='l lw 5 lt 1 lc 5 ' )
    
g._add_to_queue([item5, item6])

g('set arrow nohead from 9.5, 3.275 to 10, 3.275 lw 5 lc 0 front')
g('set label "2:1" at 9.5,3.33 front')
g('set arrow nohead from 9.5, 2.955 to 10, 2.955 lw 5 lc 0 front')
g('set label "7:3" at 9.5,3.0 front')
g('set arrow nohead from 9.5, 2.822 to 10, 2.822 lw 5 lc 0 front')
g('set label "5:2" at 9.5,2.88 front')
g('set arrow nohead from 9.5, 2.49898 to 10, 2.49898 lw 5 lc 0 front')
g('set label "3:1" at 9.5,2.55 front')
g('set arrow nohead from 9.5, 2.06 to 10, 2.06 lw 5 lc 0 front')
g('set label "4:1" at 9.5,2.1 front')
g('set arrow nohead from 0.0, 2.1 to 10, 2.1 lw 5 lt 2 lc 9 front')
g('set arrow nohead from 0.0, 3.3 to 10, 3.3 lw 5 lt 2 lc 9 front')
g.refresh()
g.close()
