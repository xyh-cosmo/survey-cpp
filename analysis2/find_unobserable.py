import os,sys
from pylab import *

if len(sys.argv) < 2:
    sys.exit(0)


d=loadtxt(sys.argv[1])

t=d[:,0]
dt=t[1:]-t[0:len(t)-1]


n,b=histogram(dt,bins=100)

b=0.5*(b[0:len(b)-1]+b[1:])


plot(b,log10(n),'+')


show()