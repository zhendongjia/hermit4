#!/usr/bin/env python
import random
import math
import os
import shutil
import sys
#from goto import goto,label

hermit_location = os.getenv('HOME')+'/Hermit4/hermit4'
runfile_location = os.getenv('HOME')+'/Hermit4/vega/test10'

begin = int(sys.argv[1])
end = int(sys.argv[2])

if not os.path.exists(runfile_location):
    os.mkdir(runfile_location)

for loop in range(begin, end + 1):
    dirname='%s/%d'%(runfile_location, loop)
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)

    print dirname

    k_start=1
    t_comp=100000.0

    n=2
    n_rand=0
    eta=1e-3
    delta_t=1e4
    t_crit=1e8
    dtmin=2.0e-5
    rmin=0.002
    etau=0.05
    gmin=1.0e-6

    kz = [1, 0, 1, 0, 0, 1, 0, 0, 1, 2]

    inc=0
    m_j=0.0009546
    m_sun_g=1.989E33
#radius in km unit
#    radius=[10**random.uniform(-3,3)]
    radius = [1.]
	
    dens_p = 1.0
    mp=[dens_p*4*math.pi*(i*1e5)**3/(3*m_sun_g) for i in radius] + [2*m_j]

#    while True:
#    	ej = random.uniform(0.1,0.9)
#    	aj = random.uniform(50,100)
#    	r_edge = aj*(1+ej)
#    	r_in = aj*(1-ej)
#	if r_in > 14 and r_edge < 110: break; 

#    ap = 14.

    aj = 90.
    ej = 0.2 
    r_in = aj*(1-ej)
    r_edge = aj*(1+ej)
    print r_in, r_edge
#    while True:
#	ap = random.uniform(10,120)
#	if ap <= r_in or ap >= r_edge: break;

    ap = r_in - 5.
    print 'radius of planetesimal = ', radius, 'location of planetesimal = ', ap
    print 'ecc of planet = ', ej, 'location of planet = ', aj 	

    a=[ap, aj]
    ecc = [0.0, ej]
    rp=[a[i]*(1-ecc[i]) for i in range(n)]
    theta=[random.uniform(0,2*math.pi) for i in range(n)]
    x=[rp[i]*math.cos(theta[i]) for i in range(n)]
    y=[rp[i]*math.sin(theta[i]) for i in range(n)]
    z=[0]*n

    g_p=2
    g_d=1
    g_r=0

    t_dep = 1e7

    dens0 = 2000
    m_crit = 1e-20
    r_esc = 200
    kg = 1.5
    cd = 0.165
    r_mhs = 9.0
    three_d = 0
    n_min = n-1

#    p=[math.sqrt(a[i]**3)/(1+mp[i]) for i in range(n)]

    m_unit = 2.0
    a_unit = 1.0
    t_unit = math.sqrt(a_unit**3/m_unit)
    v_unit = math.sqrt(m_unit/a_unit)

    vp=[math.sqrt((1+ecc[i])*(m_unit+mp[i])/rp[i]) for i in range(n)]

    vx=[-vp[i]*math.sin(theta[i]) for i in range(n)]
    vy=[vp[i]*math.cos(theta[i]) for i in range(n)]
    vz=[0]*n


    f=open('input','w')
    print >>f, k_start, t_comp
    print >>f, n, n_rand, eta, delta_t/t_unit, t_crit/t_unit

    for i in kz: print >>f, i,
    print >>f

    print >>f, dtmin, rmin, etau, gmin
    for i in range(n):
          print >>f,mp[i]/m_unit,x[i]/a_unit,y[i]/a_unit,z[i]/a_unit,vx[i]/v_unit,vy[i]/v_unit,vz[i]/v_unit
    print >>f, g_p, g_d, g_r, t_dep/t_unit, r_edge/a_unit, r_in/a_unit, dens0/(m_unit/a_unit**2), dens_p/(m_unit/a_unit**3), m_crit/m_unit, r_esc/a_unit, kg, ej, cd, r_mhs, three_d, n_min
    print >>f, m_unit, a_unit, t_unit, v_unit
    f.close()

    os.system("%s<input>output" % hermit_location)


    os.chdir('%s'%runfile_location)
