import os, sys

from pylab import *
from numpy import *
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

if len(sys.argv) < 2:
		print('usage:\n\t'+sys.argv[0]+'  ecliptic_coordinates.txt')
		sys.exit(0)

# read Ecliptic coordinates (l,b)
print('reading (l,b) coordinates from: %s'%(sys.argv[1]))
l_b = loadtxt(sys.argv[1])

lb_size = l_b.shape[0]


radius = 4

l = []
b = []
lc = []
bc = []


# generate more (l,b) coordiantes that sorround the given l-b coordinates
size_r	= 20
size_phi = 20
for i in range(lb_size):
	l.append(l_b[i][0])
	b.append(l_b[i][1])
	lc.append(l_b[i][0])
	bc.append(l_b[i][1])
	for m in range(size_r):
		r = (m+1)*radius/size_r
		for n in range(size_phi):
			phi = n*2*pi/(size_phi-1)
			l_tmp = l_b[i][0] + r*cos(phi)
			b_tmp = l_b[i][1] + r*sin(phi)
			l.append(l_tmp)
			b.append(b_tmp)

# convert l,b to np.array
l = array(l)
b = array(b)
lc = array(lc)
bc = array(bc)

figure()
m = Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90,90,30), labels=[1,0,1,1], fontsize=14)
m.drawmeridians(np.arange(-180,180,60))

cMap = ListedColormap(['b'])

#m.scatter(l,b,1,marker='x',color='b',latlon=True,alpha=0.5)

# 2-steps method:
#   1) use m to convert (l,b) to (x,y)
#   2) use m.scatter to plot (x,y), but set latlon=False
x,y = m(l,b)
m.scatter(x,y,1,marker='x',color='b',latlon=False,alpha=0.5)

# add legends to the ultra-deep survey positions
xc,yc = m(lc,bc)
dy = yc.max()*0.05
ax = gca()
for i in range(len(lc)):
    ax.text(xc[i],yc[i]+dy,'s',fontsize=12)

show()

