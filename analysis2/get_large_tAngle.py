# 提取大角度转动前后的两次指向，以及曝光时间，由此分析出最浪费时间的序列出现在何处。

import os, sys
from pylab import *


if len(sys.argv) < 2:
    print('usage: %s survey.dat'%(sys.argv[0]))
    sys.exit(0)

dat = loadtxt(sys.argv[1])

size_all = len(dat)

t = dat[:,0]
ra = dat[:,2]
dec= dat[:,1]
expTime= dat[:,16]
tAngle = dat[:,17]

fp1 = open('trans_50.stat','w')
fp2 = open('trans_75.stat','w')
fp3 = open('trans_120.stat','w')
fp4 = open('trans_150.stat','w')

size1=0
size2=0
size3=0
size4=0

for i in range(1,len(t)):
    if tAngle[i] >= 160 and tAngle[i] < 165:
        fp1.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n'%(ra[i-1],dec[i-1],ra[i],dec[i],expTime[i]))
        size1 +=1

    if tAngle[i] >= 165 and tAngle[i] < 170:
        fp2.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n'%(ra[i-1],dec[i-1],ra[i],dec[i],expTime[i]))
        size2 +=1

    if tAngle[i] >= 170 and tAngle[i] < 175:
        fp3.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n'%(ra[i-1],dec[i-1],ra[i],dec[i],expTime[i]))
        size3 +=1

    if tAngle[i] >= 175:
        fp4.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n'%(ra[i-1],dec[i-1],ra[i],dec[i],expTime[i]))
        size4 +=1

fp1.close()
fp2.close()
fp3.close()
fp4.close()

print('size1 / size_all = %g'%(size1/size_all))
print('size2 / size_all = %g'%(size2/size_all))
print('size3 / size_all = %g'%(size3/size_all))
print('size4 / size_all = %g'%(size4/size_all))
