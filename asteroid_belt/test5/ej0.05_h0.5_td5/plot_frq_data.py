#!/usr/bin/env python

from collections import Counter
import math
import os

def weight(x): return 1.0
def bin_id(x, binsize): return int(math.log10(x) / binsize)
def bin(x, binsize): return 10 ** (binsize * (bin_id(x, binsize) + 0.5))

binsize = 0.1
def stat(r):
    counter = Counter()
    for i in r: counter[bin_id(i, binsize)] += 1
    for i in sorted(counter): print i, 10**(binsize*i), 10**(binsize*(i+1)), counter[i]

rab = [float(i.split()[0]) for i in open('rab')]
ri = [float(i.split()[0]) for i in open('ri')]
print 'rab'
stat(rab)
print 'ri'
stat(ri)
sum_weight_ab = sum(weight(i) for i in rab)
sum_weight_i = sum(weight(i) for i in ri)
counter_ab = Counter()
counter_i = Counter()
for i in rab: counter_ab[bin(i, binsize)] += weight(i)
for i in ri: counter_i[bin(i, binsize)] += weight(i)

def f(x):
    base = counter_ab[x] / counter_i[x]
    error = 1 / math.sqrt(counter_ab[x]) + 1 / math.sqrt(counter_i[x])
    #return x, base
    return x, base, base * (1 - error), base * (1 + error)
data = sorted(f(x) for x in counter_ab)

f = open('r_f_err', 'w')
for i in data:
    for j in i:
        print >> f, j,
    print >> f, ''
f.close()
