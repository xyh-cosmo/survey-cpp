import os, sys

from pylab import *
from numpy import *
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

idx_deep_area = 12
idx_deep_area = 13
idx_time = 0

if len(sys.argv) < 2:
	print('usage:\n'+sys.argv[0]+' survey_result.txt')
	sys.exit(0)

fileName = str(sys.argv[1])

data = loadtxt(fileName);

#idx = data[:,0] <= 2463419.07492934
#data = data[idx,:]

udeep_ids = data[:,14]>=0
deep_ids  = data[:,14]<0

figure();
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))

x1,y1 = m(data[deep_ids,2],data[deep_ids,1]);
x2,y2 = m(data[udeep_ids,2],data[udeep_ids,1]);

cMap = ListedColormap(['b'])

m.scatter(x1,y1,1,marker='x',color='b',latlon=False,alpha=0.25);
m.scatter(x2,y2,1,marker='+',color='r',latlon=False,alpha=0.1);

show();
