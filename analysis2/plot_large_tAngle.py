from pylab import *
from numpy import *
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

d1 = loadtxt('trans_50.stat')
d2 = loadtxt('trans_75.stat')
d3 = loadtxt('trans_120.stat')
d4 = loadtxt('trans_150.stat')


figure()

subplot(2,2,1)
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))
m.scatter(d1[:,0],d1[:,1],1,marker='+',color='r',latlon=True,alpha=0.5)
m.scatter(d1[:,2],d1[:,3],1,marker='+',color='b',latlon=True,alpha=0.5)
ax=gca()
ax.set_title(str(len(d1)))


subplot(2,2,2)
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))
m.scatter(d2[:,0],d2[:,1],1,marker='+',color='r',latlon=True,alpha=0.5)
m.scatter(d2[:,2],d2[:,3],1,marker='+',color='b',latlon=True,alpha=0.5)
ax=gca()
ax.set_title(str(len(d2)))

subplot(2,2,3)
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))
m.scatter(d3[:,0],d3[:,1],1,marker='+',color='r',latlon=True,alpha=0.5)
m.scatter(d3[:,2],d3[:,3],1,marker='+',color='b',latlon=True,alpha=0.5)
ax=gca()
ax.set_title(str(len(d3)))

subplot(2,2,4)
m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
m.drawmeridians(np.arange(-180.,180.,60.))
idx=d4[:,3]*d4[:,1]>0
m.scatter(d4[:,0],d4[:,1],1,marker='+',color='r',latlon=True,alpha=0.5)
m.scatter(d4[:,2],d4[:,3],1,marker='+',color='b',latlon=True,alpha=0.5)
ax=gca()
ax.set_title(str(len(d4)))

#figure()

#ax = subplot(1,1,1,projection='3d')
#ax.scatter(d1[:,0],d1[:,1],d1[:,4],c='r',marker='.',alpha=0.5)
#ax.scatter(d1[:,2],d1[:,3],d1[:,4],c='b',marker='.',alpha=0.5)
#ax.set_zlim(0,275)

figure()

subplot(2,2,1)
plot(d1[:,2]-d1[:,0],d1[:,3]-d1[:,1],'.',label='ra-diff vs dec-diff')
legend()

subplot(2,2,2)
plot(d2[:,2]-d2[:,0],d2[:,3]-d2[:,1],'.',label='ra-diff vs dec-diff')
legend()

subplot(2,2,3)
plot(d3[:,2]-d3[:,0],d3[:,3]-d3[:,1],'.',label='ra-diff vs dec-diff')
legend()

subplot(2,2,4)
idx=abs(d4[:,3]-d4[:,1])<15
#plot(d4[idx,2]-d4[idx,0],d4[idx,3]-d4[idx,1],'.',label='ra-diff vs dec-diff')
plot(d4[idx,0],d4[idx,1],'.',label='ra-diff vs dec-diff')
plot(d4[idx,2],d4[idx,3],'x',label='ra-diff vs dec-diff')
legend()

show()
