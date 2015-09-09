#!/usr/bin/env python
import random
import math
import os
import shutil
import sys

hermit_location = os.getenv('HOME')+'/Hermit4/hermit4_s0.1_td2'
runfile_location = os.getenv('HOME')+'/Hermit4/asteroid_belt/test5/ej0.05_ws_s0.1_td2'

begin = int(sys.argv[1])
end = int(sys.argv[2])

#print 'Runing from %d to %d\n' % (begin, end)

if not os.path.exists(runfile_location): os.system('mkdir -p %s' % runfile_location)

for loop in range(begin, end+1):
    loc = '%s/%d' % (runfile_location, loop)
    if not os.path.exists(loc): os.mkdir(loc)
    os.chdir(loc)
    print loc

    k_start=1
    t_comp=100000.0

    n=3
    n_rand=0
    eta=0.003
    delta_t=5e4
    t_crit=5e6
    dtmin=1.0e-8
    rmin=0.002
    etau=0.05
    gmin=1.0e-6

    kz = [1, 0, 1, 0, 0, 1, 0, 0, 1, 2]

    g_p=1
    g_d=1
    g_r=0
    dens_p_temp=5.5
    t_dep=1e6
    r_edge=10.5
    r_in = 4.5
    dens0=8000
    m_crit=1e-5
    r_esc=20.0
    kg = 1.0
    cd = 0.165
    r_mhs = 1.4

    inc=0
    m_j=0.0009546
    m_s=0.000285716656
    m_sun_g=1.989E33
#radius in km unit

    def y(x,k):
        return x**(-k)/(-k)
    def g(x,k):
        return (x*(-k))**(-1.0/k)


    alpha = 0
    radius = [0]*(n-1)
    for i in range(n-1):
	 if alpha > 0:
		xmax = 200
		xmin = 20
		temp1 = random.uniform(0, 1)
         	temp2 = temp1*y(xmin,alpha) + (1-temp1)*y(xmax,alpha)
	 	radius[i] = g(temp2,alpha)
	 else:
		radius[i] = 10**(random.uniform(1.0,3.0))
    if radius[0] >= 100: dens_p = dens_p_temp
    else: 
	dens_p = dens_p_temp*radius[0]/100.0
   	if dens_p < 1.0: dens_p = 1.0

  #  radius = [10.]
    mp=[dens_p*4*math.pi*(radius[i]*1e5)**3/(3*m_sun_g) for i in range(n-2)] + [m_j] +[m_s]
    ej = 0.05
    es = 0.00
    ecc=[0.000]*(n-2) + [ej] + [es] 
    aj = 5.2
    a=[random.uniform(3.5**(-0.5), 2**(-0.5))**(-2)]+[aj]+[9.5]
 #   a = [3.5, 5.2, 9.5]
    rp=[a[i]*(1-ecc[i]) for i in range(n)]
    theta=[random.uniform(0.0, 2*math.pi)]*n
    x=[rp[i]*math.cos(theta[i]) for i in range(n)]
    y=[rp[i]*math.sin(theta[i]) for i in range(n)]
    z=[0]*n


    print 'radius = ', radius[0], 'location = ', a[0] 

    vp=[math.sqrt((1+ecc[i])*(1+mp[i])/rp[i]) for i in range(n)]

    vx=[-vp[i]*math.sin(theta[i]) for i in range(n)]
    vy=[vp[i]*math.cos(theta[i]) for i in range(n)]
    vz=[0]*n

    p=[math.sqrt(a[i]**3)/(1+mp[i]) for i in range(n)]

    three_d = 0
    nleast = n-1
    twogg = 0

    t_s = 5e5
    m_s_int = 0.0001
    r_s = 10.5
    r_edge_int = r_edge

    f=open('input','w')
    print >>f, k_start, t_comp
    print >>f, n, n_rand, eta, delta_t, t_crit

    for i in kz: print >>f, i,
    print >>f

    print >>f, dtmin, rmin, etau, gmin
    for i in range(n):
          print >>f,mp[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]
    print >>f, g_p, g_d, g_r, t_dep, r_edge, r_in, dens0, dens_p, m_crit, r_esc, kg, ej, cd, r_mhs, three_d, nleast, twogg
    if twogg > 0:
        print >>f, m_s, m_s_int, t_s, r_s, r_edge_int 
    f.close()



    os.system("%s<input>output" % hermit_location)


    os.chdir('%s'%runfile_location)
