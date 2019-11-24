# Ref: https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html

import matplotlib.pyplot as plt
import numpy as np

Problem1='capacitor'
Problem2='mirror_charge'
Problem3='charge_plate'

File1=Problem1+".dat"
File2=Problem2+".dat"
File3=Problem3+".dat"

X1,Y1,Charge1,P1 = np.loadtxt(File1, unpack=True)
X2,Y2,Charge2,P2 = np.loadtxt(File2, unpack=True)
X3,Y3,Charge3,P3 = np.loadtxt(File3, unpack=True)

X1 = X1.reshape(128,128)
Y1 = Y1.reshape(128,128)
P1 = P1.reshape(128,128)

X2 = X2.reshape(128,128)
Y2 = Y2.reshape(128,128)
P2 = P2.reshape(128,128)

X3 = X3.reshape(128,128)
Y3 = Y3.reshape(128,128)
P3 = P3.reshape(128,128)


x1 = X1[0::5,0::5]
y1 = Y1[0::5,0::5]
x2 = X2[0::5,0::5]
y2 = Y2[0::5,0::5]
x3 = X3[0::5,0::5]
y3 = Y3[0::5,0::5]



# Plots direction of the electrical vector field
dX1, dY1 = np.gradient( P1 )
dX2, dY2 = np.gradient( P2 )
dX3, dY3 = np.gradient( P3 )

fig,ax=plt.subplots(3,2, figsize=(12,18))



dx1 = dX1[0::5,0::5]
dy1 = dY1[0::5,0::5]
dx2 = dX2[0::5,0::5]
dy2 = dY2[0::5,0::5]
dx3 = dX3[0::5,0::5]
dy3 = dY3[0::5,0::5]


#ax[0][0].quiver(x1, y1, dx1, dy1, color='r',  zorder=1, width=0.007, headwidth=3., headlength=4., units='x')
ax[0][0].quiver(x1, y1, dx1, dy1, color='r')
ax[0][1].contourf(X1, Y1, P1, 100, cmap='RdGy')

ax[0][0].set_title(Problem1+": E lines")
ax[0][1].set_title(Problem1+": potential")

ax[0][0].set_xlim(0.0,1.0)
ax[0][0].set_ylim(0.0,1.0)
ax[0][1].set_xlim(0.0,1.0)
ax[0][1].set_ylim(0.0,1.0)

#ax[1][0].quiver(x2, y2, dx2, dy2, color='r',  zorder=1, width=0.007, headwidth=3., headlength=4., units='x')
ax[1][0].quiver(x2, y2, dx2, dy2, color='r')
ax[1][1].contourf(X2, Y2, P2, 100, cmap='RdGy')

ax[1][0].set_title(Problem2+": E lines")
ax[1][1].set_title(Problem2+": potential")

ax[1][0].set_xlim(0.0,1.0)
ax[1][0].set_ylim(0.0,1.0)
ax[1][1].set_xlim(0.0,1.0)
ax[1][1].set_ylim(0.0,1.0)


#ax[2][0].quiver(x3, y3, dx3, dy3, color='r', units='xy', zorder=1, width=0.007, headwidth=3., headlength=4.)
ax[2][0].quiver(x3, y3, dx3, dy3, color='r')
ax[2][1].contourf(X3, Y3, P3, 100, cmap='RdGy')

ax[2][0].set_title(Problem3+": E lines")
ax[2][1].set_title(Problem3+": potential")

ax[2][0].set_xlim(0.0,1.0)
ax[2][0].set_ylim(0.0,1.0)
ax[2][1].set_xlim(0.0,1.0)
ax[2][1].set_ylim(0.0,1.0)



#plt.show()
plt.savefig('fff.png')
