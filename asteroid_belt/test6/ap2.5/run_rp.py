#!/usr/bin/env python
import random
import math
import os
import shutil
import sys

hermit_location = os.getenv('HOME')+'/Hermit4/hermit4'
runfile_location = os.getenv('HOME')+'/Hermit4/asteroid_belt/test6/ap2.5'

begin = int(sys.argv[1])
end = int(sys.argv[2])


##########################################################################################
ns = 50
A = [0]*ns
for i in range(len(A)):
  A[i] = float(math.factorial(2*i))/2**(2*i)/(math.factorial(i))**2
  A[i] = A[i]**2


def Fpd_wg(aj):
  if aj <= ain or aj >= aout:
    dens = aj**(-kg)
    dens_base = dens * (-kg)
    f = - 4 * math.pi * dens
    f_base = -4 * math.pi * dens_base
    n_tot = 0
    n_tot_base = 0
    for i in range(len(A)):
      if aj <= ain:
        n = i * A[i] * (aj/ain)**(2*i-1+kg)/(2*i-1+kg)
        n_base = n * (2*i-1+kg)
      if aj >= aout:
        n = -(2*i+1)/2. * A[i] * (aout/aj)**(2*i+2-kg)/(2*i+2-kg)
        n_base = - n * (2*i+2-kg)
      n_tot = n_tot + n
      n_tot_base = n_tot_base + n_base
    fr = f * (zk + n_tot)
    fr_base = f_base * (zk + n_tot) + f * n_tot_base
  if aj > ain and aj < aout:
    fr,fr_base = Fjd_wg(aj) 
  return fr, fr_base


def Fjd_wg(aj):
  dens = aj**(-kg) 
  dens_base = dens * (-kg)
  f = 2 * math.pi * dens
  f_base = 2 * math.pi * dens_base
  n_tot = 0
  n_tot_base = 0
  for i in range(len(A)):
    n = A[i] * (2*i*(aj/aout)**(2*i-1+kg)/(2*i-1+kg) - (2*i+1)*(ain/aj)**(2*i+2-kg)/(2*i+2-kg))
    n_base = A[i] * (2*i*(aj/aout)**(2*i-1+kg) + (2*i+1)*(ain/aj)**(2*i+2-kg))
    n_tot = n_tot + n
    n_tot_base = n_tot_base + n_base
  fr = f * n_tot
  fr_base = f_base * n_tot + f * n_tot_base
  return fr, fr_base



m_j = 0.000954786
m_s = 0.000285837
a_j = 5.202545
a_s = 9.554841

dens0 = 1700 
kg = 1.5
zk = 1.094
tdep = 1e6
ain = 4.5
aout = 11.0

f = open('ap_f_fdot','w')
ap = [i*0.0001 for i in range(1000,45000)]
for i in range(len(ap)):
    fp,fp_dot = Fpd_wg(ap[i])
    print >> f, ap[i], fp, fp_dot
f.close()

f = open('aj_f_fdot','w')
apj = [ain + (a_j-ain)*2*i/10000. for i in range(10000)]
for i in range(len(apj)):
    fj,fj_dot = Fjd_wg(apj[i])
    print >> f, apj[i], fj, fj_dot
f.close()

f = open('as_f_fdot','w')
aps = [2*a_s-aout + (aout-a_s)*2*i/10000. for i in range(10000)]
for i in range(len(aps)):
    fj,fj_dot = Fjd_wg(aps[i])
    print >> f, aps[i], fj, fj_dot
f.close()

#####################################################################################################
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
    delta_t=1e4
    t_crit=1e7
    dtmin=1.0e-8
    rmin=0.002
    etau=0.05
    gmin=1.0e-6

    kz = [1, 0, 1, 0, 0, 1, 0, 0, 1, 2]

    g_p=2
    g_d=1
    g_r=0
    dens_p_temp=5.5
    
    m_crit = 1e-5
    r_esc = 20.
    cd = 0.165
    r_mhs = 0.1

    inc=0
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
         radius[0] = float(sys.argv[3])      
    if radius[0] >= 100: dens_p = dens_p_temp
    else: 
	dens_p = dens_p_temp*radius[0]/100.0
   	if dens_p < 1.0: dens_p = 1.0


    mp=[dens_p*4*math.pi*(radius[i]*1e5)**3/(3*m_sun_g) for i in range(n-2)] + [m_j] + [m_s]
    ej = 0.05
    es = 0.00
    ecc=[0.000]*(n-2) + [ej] + [es]
    a = [2.5, a_j, a_s]
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
    r_out_int = aout

    f=open('input','w')
    print >>f, k_start, t_comp
    print >>f, n, n_rand, eta, delta_t, t_crit

    for i in kz: print >>f, i,
    print >>f

    print >>f, dtmin, rmin, etau, gmin
    for i in range(n):
          print >>f,mp[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]
    print >>f, g_p, g_d, g_r, tdep, aout, ain, dens0, dens_p, m_crit, r_esc, kg, ej, cd, r_mhs, three_d, nleast, twogg
    if twogg > 0:
        print >>f, m_s, m_s_int, t_s, r_s, r_out_int 
    f.close()

#    print 'finish write output'


    os.system("%s<input>output" % hermit_location)


    os.chdir('%s'%runfile_location)
