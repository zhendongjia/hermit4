#!/usr/bin/env python
import os
import math
import Gnuplot as gp


tdep = 1e6
vrel = 1e-2
rp = [i for i in range(10, 1000)]

#####################################################
a = 2.0

densp = [0]*len(rp)
for i in range(len(rp)):
    if rp[i] > 100: densp[i] = 5.5
    else:
        densp[i] = 5.5*rp[i]/100.
        if densp[i] < 1.0: densp[i] = 1.0

mp = [float(4./3 * math.pi * densp[i] * rp[i]**3 * 5e-19) for i in range(len(rp))]

td1 = [a**2.75*rp[i]*densp[i]/vrel*9/(2*math.pi*tdep) for i in range(len(rp))]
td2 = [(1./mp[i])*a**2*4e-4/(2*math.pi*tdep) for i in range(len(rp))]

td = [1.0/((1/td1[i])+(1/td2[i])) for i in range(len(rp))]

data = [(rp[i], td[i]) for i in range(len(rp))]




########################################################
densp = [3.0]*len(rp)

mp = [float(4./3 * math.pi * densp[i] * rp[i]**3 * 5e-19) for i in range(len(rp))]

td1 = [a**2.75*rp[i]*densp[i]/vrel*9/(2*math.pi*tdep) for i in range(len(rp))]
td2 = [(1./mp[i])*a**2*4e-4/(2*math.pi*tdep) for i in range(len(rp))]
td = [1./((1/td1[i])+(1/td2[i])) for i in range(len(rp))]

data_densp3 = [(rp[i], td[i]) for i in range(len(rp))]

a = 3.5

mp = [float(4./3 * math.pi * densp[i] * rp[i]**3 * 5e-19) for i in range(len(rp))]

td1 = [a**2.75*rp[i]*densp[i]/vrel*9/(2*math.pi*tdep) for i in range(len(rp))]
td2 = [(1./mp[i])*a**2*4e-4/(2*math.pi*tdep) for i in range(len(rp))]
td = [1./((1/td1[i])+(1/td2[i])) for i in range(len(rp))]
data_densp3_ap35 = [(rp[i], td[i]) for i in range(len(rp))]


########################################################
a = 3.5

densp = [0]*len(rp)
for i in range(len(rp)):
    if rp[i] > 100: densp[i] = 5.5
    else:
        densp[i] = 5.5*rp[i]/100.
        if densp[i] < 1.0: densp[i] = 1.0

mp = [float(4./3 * math.pi * densp[i] * rp[i]**3 * 5e-19) for i in range(len(rp))]

td1 = [a**2.75*rp[i]*densp[i]/vrel*9/(2*math.pi*tdep) for i in range(len(rp))]
td2 = [(1./mp[i])*a**2*4e-4/(2*math.pi*tdep) for i in range(len(rp))]
td = [1./((1/td1[i])+(1/td2[i])) for i in range(len(rp))]

data_ap35 = [(rp[i], td[i]) for i in range(len(rp))]
data_t = [(rp[i], td[i]*math.exp(10)) for i in range(len(rp))]







item1 = gp.Data(data, title = '', with_='l lw 5 lc 0')
item2 = gp.Data(data_densp3, title = '', with_='l lw 10 lc 9')
item3 = gp.Data(data_ap35, title = '', with_='l lw 5 lc 0')
item4 = gp.Data(data_t, title = '', with_='l lw 5 lc 0')
item5 = gp.Data(data_densp3_ap35, title = '', with_='l lw 10 lc 9')


g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "r_tdamp.eps"')
g('set xlabel "size (km)"')
g('set ylabel "Tdamp/Tdep"')
g('set log ')
g('set arrow from 10.0, 1.0 to 1000.0, 1.0 nohead lt 0 lw 10')
#g('set arrow from 20, 0.2 to 20, 1.0')
g('set arrow from 190, 5.0 to 190, 7e4')
g('set label "0 Myr to 10 Myr" at 200, 1000')
g('set label "3.5 AU" at 50, 2')
g('set label "2.0 AU" at 50, 0.05')
g.plot(item1, item2, item3, item4, item5)
g('set output')

g = gp.Gnuplot()
g('set term postscript eps solid color')
g('set output "r_tdamp_test.eps"')
g('set xlabel "size (km)"')
g('set ylabel "Tdamp/Tdep"')
g('set log ')
g('set arrow from 10.0, 1.0 to 1000.0, 1.0 nohead lt 0 lw 10')
#g('set arrow from 20, 0.2 to 20, 1.0')
g('set arrow from 200, 3.0 to 200, 7e4')
g('set label "0 Myr to 10 Myr" at 210, 1000')
g('set label "3.5 AU" at 50, 2')
g('set label "2.0 AU" at 50, 0.05')
g.plot(item1, item3)
g('set output')
