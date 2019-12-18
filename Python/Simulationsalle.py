
import time
import numpy as np
import pyroomacoustics as pra
import matplotlib.pyplot as plt
from scipy.io import wavfile
##
pra.constants.set('c', 345)

fs, s1 = wavfile.read('/home/felix/Simulation salle/Test_salle_1m/ref_500Hz_pure.wav')
fs, s2 = wavfile.read('/home/felix/Simulation salle/Test_salle_1m/ref_880Hz_pure.wav')

#fs, s1 = wavfile.read('/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.25m/melange-01.wav')
#fs, s2 = wavfile.read('/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.25m/melange-02.wav')

#fs, s1 = wavfile.read('/home/felix/Simulation salle/Test reberv = 0,2 et 5 images/ref_1200Hz_pure.wav')
#fs, s2 = wavfile.read('/home/felix/Simulation salle/s1-02.wav')

# room dimension
room_dim = [5, 5]

# Create the shoebox
shoebox = pra.ShoeBox(
    room_dim,
    absorption=0.9999,
    fs=fs,
    max_order=0,
    )

# source and mic locations



shoebox.add_source([2, 3.01], signal=s1)
shoebox.add_source([2.02, 2.99], signal=s2)


R = np.array([[2,2.02],[3,3]])
bf = pra.MicrophoneArray(R,shoebox.fs)
shoebox.add_microphone_array(bf)


# run ism
shoebox.simulate()

shoebox.plot()
plt.grid()
plt.show()

audio_reverb = shoebox.mic_array.to_wav('/home/felix/test_reverb.wav', norm=True, bitdepth=np.int16)
##
#fs, m1 = wavfile.read('/home/felix/0,5m-01.wav')
#fs, m2 = wavfile.read('/home/felix/0,5m-02.wav')

fs,m1 = wavfile.read('/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-01.wav')
fs, m2 = wavfile.read('/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-02.wav')

def retarder(m1,n):
    p = len(m1)
    res = np.zeros(p+n)
    for i in range(p):
        res[n+i] = m1[i]
    return res

def trouver_retard(m1,m2):
    i = 0
    cond = True
    while cond:
        if m1[i] == 0:
            i = i +1
        else:
            cond = False
    
    j = 0
    cond = True
    while cond:
        if m2[j] == 0:
            j = j +1
        else:
            cond = False
    return j - i,i,j

##
f, t, Zxx1 = stft(m1, fs, nperseg=4000) # Utilise scipy
f, t, Zxx2 = stft(m2, fs, nperseg=4000) # Utilise scipy
#f, t, Zxx2 = stft(m2,nperseg=4000) # Utilise scipy

#plt.pcolormesh(t, f, np.abs(Zxx), vmin=0, vmax=amp)
plt.pcolormesh(t, f, np.abs(Zxx1))
plt.title('STFT Magnitude')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
#plt.yscale('log')
plt.axis([0,t[len(t)-1],0,5000])
plt.show()

##
def modifier_stft(Zxx,L,M):
    nb_freq,nb_tps = Zxx.shape
    res = np.zeros(np.shape(Zxx),dtype=complex)
    for k in L:
        f = k[0]
        t = int(M/2) * k[1]
        for j in range(M):
            res[f][t+j] = Zxx[f][t+j]
    return res

##
L = visi[0]
#L = visi[1]
Zxx_m1 = modifier_stft(Zxx_x1,L,20)
tps1,res1 = istft(Zxx_m1,freq_echant,nperseg = 4000)
#tps1,res1 = istft(Zxx1,fs,nperseg = 4000)

write("/home/felix/1freq_m1.wav",freq_echant,res1)

Zxx_m2 = modifier_stft(Zxx_x2,L,20)
tps2,res2 = istft(Zxx_m2,freq_echant,nperseg = 4000)
#tps2,res2 = istft(Zxx2,fs,nperseg = 4000)
write("/home/felix/1freq_m2.wav",freq_echant,res2)


## ANCIEN
Zxx_m1 = modifier_stft(Zxx1,24,3,158,20)
tps1,res1 = istft(Zxx_m1,fs,nperseg = 4000)
#tps1,res1 = istft(Zxx1,fs,nperseg = 4000)

write("/home/felix/1freq_m1.wav",fs,res1)

Zxx_m2 = modifier_stft(Zxx2,24,3,158,20)
tps2,res2 = istft(Zxx_m2,fs,nperseg = 4000)
#tps2,res2 = istft(Zxx2,fs,nperseg = 4000)
write("/home/felix/1freq_m2.wav",fs,res2)


##
x = [i for i in range(len(res1))]
plt.plot(x,res1,linestyle = 'None',marker = 'o',markersize = 1)
plt.plot(x,res2,linestyle = 'None',marker = 'o',markersize = 1)
plt.grid()
plt.show()

##
write("/home/felix/1freq_m1_red.wav",fs,res1_red)

write("/home/felix/1freq_m2_red.wav",fs,res2_red)

##
corr = correlate(res1,res2)
abs = [i for i in range(-len(res1)+1,len(res2))]

plt.plot(abs,corr)
plt.grid()
plt.ylabel('Corrélation croisée (x1,x2)')
plt.xlabel('Déphasage s-x1 par rapport à s-x1')
plt.xlim(xmax=20000)  
plt.xlim(xmin=-20000)
plt.show()
##
def nonnul(L):
    i = 0
    j = len(L) - 1
    for k in L:
        if k != 0:
            break
        else:
            i = i + 1
    for k in L[::-1]:
        if k != 0:
            break
        else:
            j = j - 1
    return i,j
##
def reduire1(L):
    i,j = nonnul(L)
    n = len(L)
    newres = np.zeros(j - i + 2*n)
    for k in range(j-i+1):
        newres[n+k] = L[i + k]
    return newres
    
##
def reduire2(L):
    i,j = nonnul(L)
    newres = np.zeros(j - i + 1)
    for k in range(j-i+1):
        newres[k] = L[i + k]
    return newres

def nextpow2(i):
    n = 1
    count = 0
    while n < i:
        n = n*2
        count = count + 1
    return count


## GCC-PHAT
L1 = [1,2,4]
L2 = [-3,2,7]

##
corr = correlate(res2_red,res1_red)
abs = [i for i in range(-len(res1_red)+1,len(res1_red))]

plt.plot(abs,corr)
plt.grid()
#plt.xlim(xmax=2000)  
#plt.xlim(xmin=-2000)
plt.show()


##
res1_red = reduire2(res1)
res2_red = reduire2(res2)
res1_fft = np.fft.fft(res1_red)
res2_fft = np.fft.fft(res2_red)
prod = np.conjugate(res1_fft) * res2_fft
sign = 1/(np.abs(prod)) * prod
GCC = np.fft.ifft(prod)
GCC_PHAT = np.fft.ifft(sign)

##
abs = [i for i in range(-int(len(GCC)/2),int(len(GCC)/2))]

plt.plot(abs,np.fft.fftshift(np.real(GCC)))
#plt.plot(abs,np.real(GCC))

plt.grid()
#plt.xlim(xmax=2000)  
#plt.xlim(xmin=-2000)
plt.show()
