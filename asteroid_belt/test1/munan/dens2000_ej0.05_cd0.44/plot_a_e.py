#!/usr/bin/env python

import math
import os
import Gnuplot
Gnuplot.GnuplotOpts.default_term = 'post enh color eps'

color_min=1e100
color_max=0

g = Gnuplot.Gnuplot()
g('set term post enh color eps')
g('set output "a_e_r.eps"')
g.set_range('xrange', (2, 3.5))
#g.set_range('yrange', (1.0, 4.5))
g.set_range('cbrange', (10, 3000))
g('set pm3d')
g('set logscale zcb')


f_data = open('r_a0_a_e_t_e0_few')
grps = [i.split('\n')
        for i in f_data.read().split('\n\n')
        if len(i) > 1]
f_data.close()
items = []
for grp in grps:
    if len(grp) == 0 or len(grp[0].split()) < 5: continue
    color = float(grp[0].split()[0])
    data = [(float(i.split()[2]), float(i.split()[3])) for i in grp]
    items.append(Gnuplot.PlotItems.Data(
        data,
        with_='lp lt 1 lc palette frac %f' % ((math.log10(color)-1.0)/2.8)))
g._add_to_queue(items)
g.refresh()
g.close()
