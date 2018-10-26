import os, sys

from pylab import *
from numpy import *
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

area1 = loadtxt('area_normal.txt')
area2 = loadtxt('area_new.txt')

ra1,dec1 = area1[:,0], area1[:,1]
ra2,dec2 = area2[:,0], area2[:,1]

figure(figsize=(6,4));
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))
m.scatter(ra1,dec1,1,marker='.',color='b',latlon=True,alpha=0.5)
m.scatter(ra2,dec2,1,marker='.',color='g',latlon=True,alpha=0.5)

#show()

fig = gcf()
fig.subplots_adjust(left=0.075,right=0.9,top=1,bottom=0.0)
fig.savefig('sky_map.png')

show()
