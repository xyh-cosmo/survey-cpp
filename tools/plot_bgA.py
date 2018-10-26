from pylab import *
from matplotlib import cm

d = loadtxt('beta_bgal_area.txt')

X,Y = meshgrid(linspace(10,20,41),linspace(10,20,41))
Z   = d[:,2].reshape(41,41)
Ls  = linspace(Z.min(),Z.max(),20)

fig = figure(figsize=(6,6))
ax = fig.add_subplot(1,1,1)

cntr = ax.contourf(X,Y,Z,levels=Ls,cmap=cm.coolwarm)

fig.colorbar(cntr,shrink=0.5,aspect=5)

show()
