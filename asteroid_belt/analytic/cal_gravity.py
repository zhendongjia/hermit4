#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
import Gnuplot as gp


#---------------------------------------------gravity test----------------------------------
k = 1.5

a_g = [5.2, 9.58]
ap = [i*0.01 for i in range(200,400)]
nstep = 200

gap_out = 11.0
gap_in = 4.5
disk_out = 120.
disk_in = 0.01

grav_gap = [((0.5/(1+k)*(i/gap_out)**(1+k)
             + 0.5625/(3+k)*(i/gap_out)**(3+k)
             - 0.75/(4-k)*(gap_in/i)**(4-k)
             - 0.7/(6-k)*(gap_in/i)**(6-k)
             - 1.0/(2-k)*(gap_in/i)**(2-k))*i**(-k)) for i in ap] 

grav_bound = [((- 0.5/(1+k)*(i/disk_out)**(1+k)
               - 0.5625/(3+k)*(i/disk_out)**(3+k)
               + 0.75/(4-k)*(disk_in/i)**(4-k)
               + 0.7/(6-k)*(disk_in/i)**(6-k)
               + 1.0/(2-k)*(disk_in/i)**(2-k)
               - 1./(2-k)
                + 1.25*(1-k)/((4-k)*(1+k))
                + 1.265625*(1-k)/((6-k)*(3+k)))*i**(-k)) for i in ap]

grav_comb = [((0.5/(1+k)*(i/gap_out)**(1+k)
             + 0.5625/(3+k)*(i/gap_out)**(3+k)
             - 0.75/(4-k)*(gap_in/i)**(4-k)
             - 0.7/(6-k)*(gap_in/i)**(6-k)
             - 1.0/(2-k)*(gap_in/i)**(2-k)
               - 0.5/(1+k)*(i/disk_out)**(1+k)
               - 0.5625/(3+k)*(i/disk_out)**(3+k)
               + 0.75/(4-k)*(disk_in/i)**(4-k)
               + 0.7/(6-k)*(disk_in/i)**(6-k)
               + 1.0/(2-k)*(disk_in/i)**(2-k))*i**(-k)) for i in ap]


data_gap = [(ap[i], grav_gap[i]) for i in range(nstep)]
data_bound = [(ap[i], grav_bound[i]) for i in range(nstep)]
data_comb = [(ap[i], grav_comb[i]) for i in range(nstep)]


item_gap = gp.Data(data_gap, title = 'no disk edge', with_='l lw 10')
item_bound = gp.Data(data_bound, title = 'with disk edge', with_='l lw 10')
item_comb = gp.Data(data_comb, title = 'combine disk edge and disk gap', with_='l lw 10')


g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "ap_grav.eps"')
g('set key left top')
g('set autoscale')
g('set xlabel "a (AU)"')
g('set ylabel "gravity"')
#g('set yrange [1e-4:]')
#g('set log y')
g('set grid')
g.plot(item_gap, item_bound, item_comb)
g('unset log')
g('set output')
