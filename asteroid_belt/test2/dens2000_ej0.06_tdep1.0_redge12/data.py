#!/usr/bin/env python
import os
import math

bin_width = 0.2

data_dir = range(20000,30000)


open_tidal = 1
open_1 = 1
open_2 = 0
open_3 = 1
open_escape = 1


class Orbit:
        def __init__(self, line, time):
                segments = line.split()
                if segments[0] == 'ORBIT':
                        self.name = int(segments[9])
                        self.time = time
                        self.ecc = float(segments[10])
                        self.r = float(segments[11])
                        self.a = float(segments[12])
                        self.w = float(segments[14])
                        self.x = float(segments[15])
                        self.y = float(segments[16])
                        self.vx = float(segments[17])
                        self.vy = float(segments[18])
                        self.t_tidal1 = toFloat(segments[19])
                        self.t_tidal2 = toFloat(segments[20])
                        self.radius = toFloat(segments[21])
                        self.mp = toFloat(segments[22])
                        self.fd = toFloat(segments[23])
                elif segments[0] == 'ESCAPE':
                        self.name = int(segments[8])
                        self.ecc = float(segments[10])
                        self.a = abs(float(segments[12]))
                        self.r = float(segments[11])
                        self.time = time
        


def toFloat(s):
    try: return float(s)
    except: return None


def get_delt_w(i,j):
        if j==0: i[j].delt_w = i[j].w - i[j].w
        if j>0:
           if i[j].w < 0: i[j].w = i[j].w + 360
           if i[j-1].w < 0: i[j-1].w = i[j-1].w + 360
           i[j].delt_w = i[j].w - i[j-1].w
           if i[j].delt_w >= 180: i[j].delt_w = i[j].delt_w - 360
           if i[j].delt_w <= -180: i[j].delt_w = i[j].delt_w + 360
           return i[j].delt_w

def get_angular_momentum(i):
    l_tot = 0
    for j in i:
        if j is None: continue
        l_tot = l_tot + math.sqrt(j.a*abs((1-j.ecc**2))/(1+j.mp))*j.mp
    return l_tot


def get_w_per_time(i,j):
    if i == 0: orbit[i][j].fw = 0
    if i > 0:
       delt_t = orbit[i][j].time - orbit[i-1][j].time
       if orbit[i][j].w < 0: orbit[i][j].w=orbit[i][j].w+360
       if orbit[i-1][j].w < 0: orbit[i][j].w=orbit[i-1][j].w+360
       delt_w = orbit[i][j].w - orbit[i-1][j].w
       if delt_w >= 180: delt_w = delt_w - 360
       if delt_w <= -180: delt_w = delt_w + 360
       orbit[i][j].fw = delt_w/delt_t
       return orbit[i][j].fw


def get_resonance_a(mj,aj,dens0,fgas):
        ai = 1/(aj*(1 + mj))
        if fgas != 0.0000:
                ai = ai + mj/(4*aj**(2.0)*1.094*3.14*dens0*fgas)
                ai = 1/ai
        else: ai = None
        return ai

def get_infor(orbit, name):
  for i in orbit:
    if i[name-1] is not None:
      return i[name-1].mp, i[name-1].radius, i[name-1].t_tidal1, i[name-1].t_tidal2, i[name-1].w, i[name-1].fd

cum_gt={}
cum_gt_intl={}
if open_1 > 0:
    f0 = open('r_nab','w')
    f19 = open('r_ni','w')

if open_1 > 0:
    f18 = open('rab','w')
    f20 = open('ri','w')
    f22 = open('r_af','w')
    f23 = open('r_ai','w')
    f24 = open('af2.5_t','w')


if open_tidal == 1:
    f6 = open('r_t1_t2_td_i','w')
    f13 = open('r_t1_t2_td_f','w')



for loop in data_dir:

    if not os.path.exists('%d'%(loop)): continue
#    if os.listdir('%s' %loop) == '': continue 

    os.chdir('%s'% loop)

    print loop
 
    if not os.path.isfile('input'):
        os.chdir('..')
        continue

    if os.stat('input')[6]==0 or os.stat('output')[6]==0:
	os.chdir('..')
	continue
    
    orbit = []
    time = 0.0
    delt = float(open('input').readlines()[1].split()[3])
    n = int(open('input').readlines()[1].split()[0])


    for line in open('output').readlines():
        if line.startswith(' YRS = '):
            time = float(line.split()[2])
            orbit.append([None]*n)
        elif line.startswith(' ORBIT    I NAM ECC R A S'):
            current_orbit = Orbit(line, time)
            orbit[-1][int(current_orbit.name)-1] = current_orbit
        elif line.startswith(' ESCAPE'):
            time += 100.0
            current_orbit = Orbit(line, time)
            current_orbit.mp, current_orbit.radius, current_orbit.t_tidal1, current_orbit.t_tidal2, current_orbit.w, current_orbit.fd = get_infor(orbit, current_orbit.name)
            orbit.append([None]*n)
            orbit[-1][int(current_orbit.name)-1] = current_orbit

    f1 = open('e_t','w')
    f2 = open('a_t','w')
    f3 = open('n_t','w')
    f4 = open('q_a_Q_t','w')
    f5 = open('l_t','w')
    f16 = open('ar_t','w')
    f21 = open('radius','w')
    
    for time_orbit in orbit:
        if len(time_orbit) == 0: continue
        num = n - 1
        for i in time_orbit:
            if i is None:
                num = num - 1
                print >> f1, None,
                print >> f2, None,
                print >> f4, None, None, None,
                print >> f5, None,
                continue;
            else:
                i.p = 2*3.14*math.sqrt(i.a**3)/(1+i.mp)
                print >> f1, i.ecc,
                print >> f2, i.a,
                print >> f4, i.a*(1-i.ecc), i.a, i.a*(1+i.ecc),
                print >> f5, math.sqrt(i.a*abs((1-i.ecc**2))),
                if i.mp > 1e-4:
                        mj = i.mp
                        aj = i.a
                        dens0 = 2000*1.125e-7
                        fgas = i.fd
                        a_re = get_resonance_a(mj,aj,dens0,fgas)
                        print >> f16, a_re, i.time
        l_tot = get_angular_momentum(time_orbit)
        print >> f3, num,
        for j in range(n):
            if time_orbit[j] != None:
                print >> f1, time_orbit[j].time
                print >> f2, time_orbit[j].time
                print >> f3, time_orbit[j].time
                print >> f4, time_orbit[j].time
                print >> f5, time_orbit[j].time
                break
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f16.close()

    if open_tidal == 1:
            if len(orbit)<=1:
                    continue
            for i in orbit[1]:
                    if i is None: print >> f6, None, None, None, None
                    else:
                            if i.mp < 3e-5 and i.t_tidal1 > 0:
                                    force = 1/i.t_tidal1 + 1/i.t_tidal2
                                    print >> f6, i.radius, i.t_tidal1, i.t_tidal2, 1/force
                                    print >> f21, i.radius
            for i in orbit[-1]:
                    if i is None: print >> f13, None, None, None, None
                    elif i.t_tidal1 > 0:
                            if i.mp < 3e-5:
                                    print >> f13, i.radius, i.t_tidal1, i.t_tidal2, 1/(1/i.t_tidal1+ 1/i.t_tidal2)
    f21.close()

        
    f7 = open('fw_t','w')
    for i in range(len(orbit)):
          for j in range(len(orbit[i])):
              if orbit[i][j] is None:
                 print >> f7, None,
                 continue;
              orbit[i][j].fw = get_w_per_time(i,j)
              print >> f7, orbit[i][j].fw,
          for k in range(len(orbit[i])):
              if orbit[i][k] != None:
                 print >> f7, orbit[i][k].time
                 break
    f7.close()

    f8 = open('fd_t','w')
    for i in orbit:
    	for j in i:
	    if j is None: continue
	    else:
		print >> f8, j.fd, j.time
		break
    f8.close()



    if open_3 == 1:
            f17 = open('af_r','w')
            for i in orbit[-1]:
                if i is None:
                        print >> f17, None, None
                        continue
                if i.mp>3e-5: continue
                print >> f17, i.a, i.radius
            f17.close()




    if open_2 == 1:
            a_in = 1.5
            a_out = 4.0
            e_crit = 0.0
            f9 = open('ec0_r_a0_a_e_t','w')
            f10 = open('ec0_r_a0_dt','w')
            f11 = open('ec0_a0_e0','w')
            for j in range(n-1):
                    print >> f11, orbit[0][j].a, orbit[0][j].ecc
                    state_start = 0
                    state_stop = 0
                    for i in orbit:
                            if i[j] is None and state_start == 0: break
                            if i[j] is None and state_start == 1:
                                    t_stop = t_temp
                                    break
                            if i[j].a >= a_in and i[j].a <= a_out:
                                    if i[j].ecc >= e_crit:
                                            print >> f9, i[j].radius, orbit[0][j].a, i[j].a, i[j].ecc, i[j].time
                                    if i[j].ecc >= e_crit and state_start == 0:
                                            t_start = i[j].time
                                            state_start = 1
                                    if i[j].ecc < e_crit and state_stop == 0 and state_start == 1:
                                            t_stop = i[j].time
                                            state_stop = 1
                                            break
                            t_temp = i[j].time
                    print >> f9
                    if state_start == 1 and state_stop == 0:
                            if orbit[-1][j]is None:
                                    delt_t = t_temp - t_start
                                    print >> f10, orbit[0][j].radius, orbit[0][j].a, delt_t
                                    break
                            delt_t =  orbit[-1][j].time - t_start
                            print >> f10, orbit[0][j].radius, orbit[0][j].a, delt_t
                    if state_start == 1 and state_stop == 1:
                            delt_t = t_stop - t_start
                            print >> f10, orbit[0][j].radius, orbit[0][j].a, delt_t

            f9.close()
            f10.close()
            f11.close()



    if open_1 > 0:
             for j in orbit[-1]:
                    if j is None: continue
                    if j.a >= 2.0 and j.a <= 5.2 and j.mp < 1e-4:
                        b = int(math.log10(j.radius)/bin_width+0.5)
                        cum_gt[b] = cum_gt.get(b,0)+1
             for j in orbit[0]:
                    if j is None: continue
                    if j.a >= 2.0 and j.a <= 5.2 and j.mp < 1e-4:
                        c = int(math.log10(j.radius)/bin_width+0.5)
                        cum_gt_intl[c] = cum_gt_intl.get(c,0)+1


             for j in orbit[-1]:
                    if j is None: continue
                    if j.mp < 1e-4 and j.a <= 5.2:
			if j.a >= 2.0:                       
 				print >> f18, j.radius
                        	print loop, j.a, j.radius
                    	print >> f22, j.radius, j.a
             for j in orbit[0]:
                    if j is None: continue
                    if j.mp < 1e-4 and j.a <= 5.2:
			if j.a >= 2.0:
                        	print >> f20, j.radius
			print >> f23, j.radius, j.a


             af25 = 0
             for j in orbit[-1]:
                    if j is None: continue
                    else: 
                            if abs(j.a - 2.5) < 1e-2 and j.mp < 1e-4: af25 = 1
             if af25 == 1:
                    for j in orbit:
                            for i in j:
                                    if i is None or i.mp > 1e-4: continue
                                    print >> f24, i.a, i.time
                    print >> f24, ''



    os.chdir('..')



if open_1 > 0:
        keys = sorted(cum_gt.keys(), reverse=True)
        for i in range(1, len(keys)):
                cum_gt[keys[i]] += cum_gt[keys[i-1]]
        keys.reverse()
        for k in keys:
                print >> f0, 10**(bin_width*(k+0.5)), cum_gt[k]
        f0.close()

        keys = sorted(cum_gt_intl.keys(), reverse=True)
        for i in range(1, len(keys)):
                cum_gt_intl[keys[i]] += cum_gt_intl[keys[i-1]]
        keys.reverse()
        for k in keys:
                print >> f19, 10**(bin_width*(k+0.5)), cum_gt_intl[k]
        f19.close()


        f18.close()
        f20.close()
        f6.close()
        f13.close()
	f22.close()
	f23.close()
        f24.close()
