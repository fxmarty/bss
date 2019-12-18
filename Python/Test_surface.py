import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
## Surface
grid_z = np.zeros((6,7))
grid_z[3,3] = 1
grid_z[3,4] = 0.5
grid_z[3,5] = 1
##
    
grid_freq, grid_tps = np.meshgrid([0,1,2,3,4,5,6],[0,0.1,0.2,0.3,0.4,0.5])
    
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')

ax.plot_surface(grid_freq,grid_tps, grid_z)

plt.xlabel('fréquence (Hz)')
plt.ylabel('temps (s)')
ax.set_zlabel('res')

plt.show()

##
X = np.arange(1, 4, 0.2)
Y = np.copy(X)
X, Y = np.meshgrid(X, Y)

Z1 = np.copy(X)*4 - 8
Z2 = 2/X

fig = plt.figure()
ax = fig.gca(projection='3d')
surf1 = ax.plot_surface(X, Y, Z1, color='g')
surf2 = ax.plot_surface(X, Y, Z2, color='r')

plt.show()

##

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
import numpy as np
import mpmath
mpmath.dps = 5

# Use instead of arg for a continuous phase
def arg2(x):
    return mpmath.sin(mpmath.arg(x))

#f = lambda z: abs(mpmath.loggamma(z))
#f = lambda z: arg2(mpmath.exp(z))
#f = lambda z: abs(mpmath.besselj(3,z))
f = lambda z: arg2(mpmath.cos(z))

fig = pylab.figure()
ax = Axes3D(fig)
X = np.arange(-5, 5, 0.125)
Y = np.arange(-5, 5, 0.125)
X, Y = np.meshgrid(X, Y)
xn, yn = X.shape
W = X*0
for xk in range(xn):
    for yk in range(yn):
        try:
            z = complex(X[xk,yk],Y[xk,yk])
            w = float(f(z))
            if w != w:
                raise ValueError
            W[xk,yk] = w
        except (ValueError, TypeError, ZeroDivisionError):
            # can handle special values here
            pass
    print(xk, xn)

# can comment out one of these
ax.plot_surface(X, Y, W, rstride=1, cstride=1, cmap=cm.jet)
#ax.plot_wireframe(X, Y, W, rstride=5, cstride=5)

pylab.show()

##
freqmin = 100
freqmax = 1000

grid_z = np.zeros((20,2000))
    
grid_freq, grid_tps = np.meshgrid(freq1[0:2000],temps_seg1[0:20])
    
for j in range(len(grid_tps[0])): #freq
    for i in range(4): #tps
        a = np.log10(1/var[j][i][2])
        for k in range(4):
            grid_z[4 * i + 2 + k][j] = a
            
            
        #grid_z[0][j] = var[j][0][2]
        #grid_z[1][j] = var[j][0][2]
        #p = len(temps_seg1) - 4 * (nb_fenetre - 1) - 2 - 3

for i in range(len(grid_z)):
    for j in range(len(grid_z[0])):
        if grid_freq[i][j] < freqmin or grid_freq[i][j] > freqmax:
            grid_z[i][j] = NaN

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')

#grid_freq = np.log10(grid_freq)
#xticks = [100, 2e2, 2e3,2e4]
#ax.set_xticks(np.log10(xticks))
#ax.set_xticklabels(xticks)

import matplotlib.ticker as mticker
def log_tick_formatter(val, pos=None):
    return "{:.0e}".format(10**val)
ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))


ax.plot_surface(grid_freq,grid_tps, grid_z, cmap=cm.jet,vmin=np.nanmin(grid_z),rstride=1, cstride=1,vmin=np.nanmin(grid_z),vmax=np.nanmax(grid_z))

#,rstride=1, cstride=1
#surf = ax.plot_surface(grid_freq,grid_tps, grid_z,cstride=1,vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z))


#surf = ax.plot_surface(grid_freq,grid_tps, grid_z,shade = False,antialiased=False,cmap=cm.coolwarm,vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z)-4)
#cmap=cm.coolwarm, shade = False)
#linewidth=0,rstride=1, cstride=1

#shade=False,antialiased=False,rstride=1, cstride=1

plt.xlabel('fréquence (Hz)')
plt.ylabel('temps (s)')
ax.set_zlabel('1/variance (échelle log)')

#ax.set_xlim(np.log10(freqmin),np.log10(freqmax))
ax.set_xlim(freqmin,freqmax)

plt.show()
