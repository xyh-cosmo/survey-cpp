from pylab import *

# this is the original version used by Zhang Xin
def cmg_old(tAngle):
    tm_k=0
    cmg_use = 0
    if tAngle <= 5.0:
        cmg_use = 0
    elif tAngle <= 10.0:
        tm_k = (1/29.0-0)/(10.0-5.0)
        cmg_use = tm_k*(tAngle-5)
    elif tAngle <= 20.0:
        tm_k = (1/19.0-1/29.0)/(20.0-10.0)
        cmg_use = tm_k*(tAngle-10)+1/29.0
    elif tAngle <= 35.0:
        tm_k = (1/13.0-1/19.0)/(35.0-20.0)
        cmg_use = tm_k*(tAngle-20)+1/19.0
    elif tAngle <= 45.0:
        tm_k = (1/10.0-1/13.0)/(45.0-35.0)
        cmg_use = tm_k*(tAngle-35)+1/13.0
    elif tAngle <= 75.0:
        tm_k = (1/6.0-1/10.0)/(75.0-45.0)
        cmg_use = tm_k*(tAngle-45)+1/10.0
    elif tAngle <= 90.0:
        tm_k = (1/5.0-1/6.0)/(90.0-75.0)
        cmg_use = tm_k*(tAngle-75)+1/6.0
    elif tAngle <= 135.0:
        tm_k = (1/3.0-1/5.0)/(135.0-90.0)
        cmg_use = tm_k*(tAngle-90)+1/5.0
    elif tAngle <= 180.0:
        tm_k = (1/2.0-1/3.0)/(180.0-135.0)
        cmg_use = tm_k*(tAngle-135)+1/3.0
    
    return cmg_use

# newer version used by Youhua Xu
def cmg_new(tAngle):
    cmg_use = 0
    if tAngle < 5.0:
        cmg_use = 0
    elif tAngle >=  5.0 and tAngle < 10.0:
        cmg_use = 1./29.
    elif tAngle >= 10.0 and tAngle < 20.0:
        cmg_use = 1./19.
    elif tAngle >= 20.0 and tAngle < 35.0:
        cmg_use = 1./13.
    elif tAngle >= 35.0 and tAngle < 45.0:
        cmg_use = 1./10.
    elif tAngle >= 45.0 and tAngle < 75.0:
        cmg_use = 1./6.
    elif tAngle >= 75.0 and tAngle < 90.0:
        cmg_use = 1./5.
    elif tAngle >= 90.0 and tAngle < 135.0:
        cmg_use = 1./3.
    elif tAngle >= 135.0 and tAngle <= 180.0:
        cmg_use = 1./2.

    return cmg_use    

# new version-2 under construction
def cmg_new2_exp(tAngle,use_a,use_b,Angle_a,Angle_b,idx=5.0):
    factor = idx*(tAngle-Angle_a)/(Angle_b-Angle_a)
    cmg_use = use_a
    cmg_use += (use_b-use_a)*(1 - exp(-factor) + exp(-idx))
    return cmg_use

def cmg_new2(tAngle,idx=5):
    cmg_use = 0
    if tAngle <= 5.0:
        cmg_use = 0
    elif tAngle <= 10.0:
        cmg_use = cmg_new2_exp(tAngle,0,1/29.,5,10,idx=idx)
    elif tAngle <= 20.0:
        cmg_use = cmg_new2_exp(tAngle,1./29,1./19,10,20,idx=idx)
    elif tAngle <= 35.0:
        cmg_use = cmg_new2_exp(tAngle,1./19,1./13,20,35,idx=idx)
    elif tAngle <= 45.0:
        cmg_use = cmg_new2_exp(tAngle,1./13,1./10,35,45,idx=idx)
    elif tAngle <= 75.0:
        cmg_use = cmg_new2_exp(tAngle,1./10,1./6,45,75,idx=idx)
    elif tAngle <= 90.0:
        cmg_use = cmg_new2_exp(tAngle,1./6,1./5,75,90,idx=idx)
    elif tAngle <= 135.0:
        cmg_use = cmg_new2_exp(tAngle,1./5,1./3,90,135,idx=idx)
    elif tAngle <= 180.0:
        cmg_use = cmg_new2_exp(tAngle,1./3,1./2,135,180,idx=idx)

    return cmg_use

angle = linspace(0,180,500)

cmg1 = []
cmg2 = []
cmg3 = []
cmg4 = []
cmg5 = []

fp=open('cmg_use0.txt','w')

for i in range(len(angle)):
    cmg1.append(cmg_old(angle[i]))
    fp.write('%5d  %15.12f\n'%(angle[i],cmg_old(angle[i])))
    cmg2.append(cmg_new(angle[i]))
    cmg3.append(cmg_new2(angle[i],5))
    cmg4.append(cmg_new2(angle[i],7.5))
    cmg5.append(cmg_new2(angle[i],10))

fp.close()

#plot(angle, cmg1, '-.r', label='CMG-old')
#plot(angle, cmg2, ':b', label='CMG-new')
#plot(angle, cmg3, '-g', label='CMG-new2,idx=5.0')
#plot(angle, cmg4, '--m', label='CMG-new2,idx=7.5')
#plot(angle, cmg5, '--c', label='CMG-new2,idx=10.')
#xlabel(r'$\theta$',fontsize=15)
#ylabel(r'CMG cost',fontsize=15)
#legend()

#show()
