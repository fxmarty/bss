freq_echant, x1 = read("/home/felix/Simulation salle/Test_salle_0.1m/m-01.wav")
freq_echant, x2 = read("/home/felix/Simulation salle/Test_salle_0.1m/m-02.wav")
#freq_echant, x1 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance1m/m-01.wav")
#freq_echant, x2 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance1m/m-02.wav")
nb_echantillons = len(x1)

temps = np.linspace(0,nb_echantillons/freq_echant,nb_echantillons)

freq1,temps_seg1,Zxx_x1 = stft(x1,freq_echant,nperseg = 4000)
freq2,temps_seg2,Zxx_x2 = stft(x2,freq_echant,nperseg = 4000)
stft_x1 = np.abs(Zxx_x1)
stft_x2 = np.abs(Zxx_x2)
## Sources
freq_echant, s1 = read("/home/felix/Flamenco1_U89_cut.wav")
freq_echant, s2 = read("/home/felix/Flamenco2_U89_cut.wav")
nb_echantillons = len(s1)
temps = np.linspace(0,nb_echantillons/freq_echant,nb_echantillons)

freq1,temps_seg1,Zxx_s1 = stft(s1,freq_echant,nperseg = 4000)
freq2,temps_seg2,Zxx_s2 = stft(s2,freq_echant,nperseg = 4000)
stft_s1 = np.abs(Zxx_s1)
stft_s2 = np.abs(Zxx_s2)

##
x1,x2 = melange(s1,s2,0.5,0.8,0.7,0.6)
write("/home/felix/melange1.wav",freq_echant,x1)
write("/home/felix/melange2.wav",freq_echant,x2)


##
melanges_stft = np.array([np.abs(Zxx_x1),np.abs(Zxx_x2)])
stft_melanges_comp = np.array([Zxx_x1,Zxx_x2])
melanges = np.array([x1,x2])

rapp = rapport(melanges_stft[1],melanges_stft[0],freq_echant)
var,nb_fenetre = variance(rapp,20)

var_trie = trier(var)
visi,indmax,rappseul = detecter_sources_visibles(rapp,var_trie,5e-3,0.2,20,2)
#rappseul[rappseul[:,2].argsort()]

mat_mixage = creer_mat_mixage(2,visi,rapp,20)
mat_sourcesextraites = retrouver_sources(mat_mixage,melanges)
##
for i in range(0,15):
    #print(rappseul[i][2])
    print(freq1[rappseul[i][0]])

#freq1[rappseul[i][0]],,var_trie[i]

##
mat_mixage = creer_mat_mixage(2,visi,rapp,20)


stft_res = np.zeros((2,len(freq1),len(temps_seg1)),dtype=complex)
for fk in range(len(freq1)):
    mat_mixage_inv = np.linalg.inv(np.array([[1,1],[mat_mixage[1][0]*np.exp(-1j*2*np.pi*(0)*freq1[fk]/freq_echant),mat_mixage[1][1]*np.exp(-1j*2*np.pi*(37)*freq1[fk]/freq_echant)]]))
    L = np.dot(mat_mixage_inv,stft_melanges_comp[:,fk])
    stft_res[0][fk] = L[0]
    stft_res[1][fk] = L[1]


tps1,res1 = istft(stft_res[0],freq_echant,nperseg = 4000)
tps2,res2 = istft(stft_res[1],freq_echant,nperseg = 4000)

write("/home/felix/sourceseparee1_1m.wav",freq_echant,res1)
write("/home/felix/sourceseparee2_1m.wav",freq_echant,res2)


##
wav_sources(mat_sourcesextraites,freq_echant)
##
def wav_sources(mat_sourcesextraites,freq_echant):
    N = len(mat_sourcesextraites)
    for i in range(0,N):
        maax = np.max(np.abs(mat_sourcesextraites[i]))
        mat_sourcesextraites[i] = mat_sourcesextraites[i]/(maax + 0.3)
        write("/home/felix/TIPE/signal_separe%d.wav" %i,freq_echant,
        mat_sourcesextraites[i])

## Surface
freqmin = 0
freqmax = grid_freq[0][-1]
temps_bis = temps_seg1[:170]
freq_bis = freq1[:100]
grid_z = np.zeros((len(temps_bis),len(freq_bis)))
    
grid_freq, grid_tps = np.meshgrid(freq_bis,temps_bis)
    
for j in range(len(grid_tps[0])): #freq
    for i in range(nb_fenetre): #tps
        a = np.log10(1/var[j][i][2])
        for k in range(10):
            grid_z[10 * i + k][j] = a
            
"""         
for i in range(len(grid_z)):
    for j in range(len(grid_z[0])):
        if grid_freq[i][j] < freqmin or grid_freq[i][j] > freqmax:
            grid_z[i][j] = NaN
"""
fig = plt.figure()
ax = fig.gca(projection='3d')

import matplotlib.ticker as mticker
def log_tick_formatter(val, pos=None):
    return "{:.0e}".format(10**val)

ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))


ax.yaxis.set_ticks([0,2,4,6])
##

light = LightSource(azdeg=-20, altdeg=60)
#green_surface = light.shade_rgb(rgb * green, grid_z)

illuminated_surface = light.shade(grid_z, cmap=cm.jet)
ax.plot_surface(grid_freq,grid_tps, grid_z,rstride=1,cstride=1,antialiased=False, facecolors=illuminated_surface)
#, cmap=cm.jet
#
plt.xlabel('\nfréquence (Hz)', linespacing=0.8)
plt.ylabel('\ntemps (s)', linespacing=0.8)

ax.tick_params(axis='z', which='major', pad=11)
ax.set_zlabel('\n1/variance (échelle log)',linespacing=5)
#,vmin=np.nanmin(grid_z)-100,vmax=np.nanmax(grid_z)

#ax.set_xlim(np.log10(freqmin),np.log10(freqmax))
ax.set_xlim(freqmin,freqmax)
ax.set_zlim(-2,np.max(grid_z))
plt.show()

##
from matplotlib.colors import LightSource
from matplotlib import cm

##
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(.5*R)

Z[X+Y>4.] = np.nan  # the diagonal slice

surf = ax.plot_surface(X, Y, Z, cmap=cm.binary,vmin=np.nanmin(Z), vmax=np.nanmax(Z))

ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
##
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

x = linspace(-5, 5, 200)
y = x
X,Y = meshgrid(x, y)
Z = bivariate_normal(X, Y)

fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap = cm.binary)
plt.show()

##
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xp,yp,zp = [],[],[]

for i in range(len(rapp)):
    for j in range(len(rapp[0])):
        xp.append(freq1[i])
        yp.append(temps_seg1[j])
        zp.append(1/var[i][j][2])

ax.scatter(tps,fre,zp)
plt.axis([0,10,0,1000])
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

plt.show()