#!/usr/bin/env python
import random
import math
import os
import shutil
import sys

hermit_location = os.getenv('HOME')+'/Hermit4/hermit4'
runfile_location = os.getenv('HOME')+'/Hermit4/asteroid_belt'

begin = int(sys.argv[1])
end = int(sys.argv[2])

print 'Runing from %d to %d\n' % (begin, end)

for loop in range(begin-1, end):

    if not os.path.exists('%d'%(loop+1)):
        os.mkdir('%d'%(loop+1))
    os.chdir('%s/%d'%(runfile_location,(loop+1)))

    print '%s/%d'%(runfile_location,(loop+1))

    k_start=1
    t_comp=100000.0

    n=2
    n_rand=0
    eta=0.003
    delta_t=1e5
    t_crit=1e7
    dtmin=1.0e-8
    rmin=0.002
    etau=0.05
    gmin=1.0e-6

    kz = [1, 0, 1, 0, 0, 1, 0, 0, 1, 2]

    g_p=1
    g_d=1
    g_r=0
    dens_p=5.5
    t_dep=1e6
    r_edge=10
    r_in = 4.5
    dens0=2000
    m_crit=1e-20
    r_esc=10.0
    kg = 1.5
    cd = 0.165

    inc=0
    m_j=0.0009546
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
		radius[i] = 10**(random.uniform(1.0,3.3))
         
    mp=[dens_p*4*math.pi*(radius[i]*1e5)**3/(3*m_sun_g) for i in range(n-1)] + [m_j]
    ej = 0.02
    ecc=[0.000]*(n-1) + [ej]
    a=[random.uniform(3.5**(-0.5), 2**(-0.5))**(-2)]+[5.2]
    rp=[a[i]*(1-ecc[i]) for i in range(n)]
    theta=[random.uniform(0.0, 2*math.pi)]*n
    x=[rp[i]*math.cos(theta[i]) for i in range(n)]
    y=[rp[i]*math.sin(theta[i]) for i in range(n)]
    z=[0]*n

    vp=[math.sqrt((1+ecc[i])*(1+mp[i])/rp[i]) for i in range(n)]

    vx=[-vp[i]*math.sin(theta[i]) for i in range(n)]
    vy=[vp[i]*math.cos(theta[i]) for i in range(n)]
    vz=[0]*n

    p=[math.sqrt(a[i]**3)/(1+mp[i]) for i in range(n)]


    f=open('input','w')
    print >>f, k_start, t_comp
    print >>f, n, n_rand, eta, delta_t, t_crit

    for i in kz: print >>f, i,
    print >>f

    print >>f, dtmin, rmin, etau, gmin
    for i in range(n):
          print >>f,mp[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]
    print >>f, g_p, g_d, g_r, t_dep, r_edge, r_in, dens0, dens_p, m_crit, r_esc, kg, ej, cd
    f.close()



    os.system("%s<input>output" % hermit_location)


    os.chdir('%s'%runfile_location)
