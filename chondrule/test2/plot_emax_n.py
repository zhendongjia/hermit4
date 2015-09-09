#!/usr/bin/env python

import os

minVal = 0.0
maxVal = 1.0
nBin = 10
binSize = (maxVal - minVal) / nBin

def bin(x): return int((x - minVal) / binSize)
def binVal(i): return (i + 0.5) * binSize + minVal

subdir = 10, 50, 100, 500, 1000, 2000
data = [ ]
for d in subdir:
    cnt = [ 0 ] * nBin
    for line in open('./rp%d/emax' % d):
        val = float(line)
        if val >= minVal and val <= maxVal:
            cnt[bin(val)] += 1
    data.append(cnt)
fout = open('emax_n0.dat', 'w')
print >> fout, 'number', ' '.join('"r = %d km"' % d for d in subdir)
for i in range(nBin):
    print >> fout, binVal(i), ' '.join(str(cnt[i]) for cnt in data)
fout.close()
os.system('./test.gplt')

