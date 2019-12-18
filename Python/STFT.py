import matplotlib
import numpy as np
from scipy import signal
##
freq, tps, stft_x = signal.stft(x, fs = 44100,nperseg = 2000)
plt.pcolormesh(tps, freq, np.abs(stft_x), vmin=0, vmax=1)
plt.title('STFT Magnitude')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

##
fs = 44100
N = fs*10
amp = 2 * np.sqrt(2)
noise_power = 0.01 * fs / 2
time = np.arange(N) / float(fs)
mod = 500*np.cos(2*np.pi*0.25*time)
carrier = amp * np.sin(2*np.pi*3e3*time + mod)
noise = np.random.normal(scale=np.sqrt(noise_power), size=time.shape)
noise *= np.exp(-time/5)
x = carrier + noise
##
plt.plot(time,x)
plt.show()
plt.axis()


##
f, t, Zxx = signal.stft(x, fs, nperseg=3000) # Utilise scipy
#plt.pcolormesh(t, f, np.abs(Zxx), vmin=0, vmax=amp)
plt.pcolormesh(t, f, np.abs(Zxx))
plt.title('STFT Magnitude')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
#plt.yscale('log')
plt.axis([0,10,1,6000])
plt.show()

## Stft mélanges
#plt.pcolormesh(temps_seg1, freq1, stft_x1) #para vmin,vmax
plt.pcolormesh(temps_seg2, freq2, stft_x2)
plt.title('STFT Magnitude')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
#plt.yscale('log')
plt.axis([0,temps_seg1[-1],100,2000])
plt.show()

## STFT sources
plt.pcolormesh(temps_seg1, freq1, stft_s1) #para vmin,vmax
#plt.pcolormesh(temps_seg2, freq2, stft_s2)
plt.title('')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
#plt.yscale('log')
plt.axis([0,temps_seg1[-1],0,5000])
plt.show()

##
matplotlib.rcParams.update({'font.size': 13})
## Sources
plt.figure(1)
plt.subplot(211)
plt.pcolormesh(temps_seg1, freq1, np.abs(Zxx_m1))
#plt.pcolormesh(temps_seg1, freq1, stft_x1)
#plt.legend()
plt.ylabel('Fréquence (Hz)')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
plt.axis([0,temps_seg1[-1],0,1300])

plt.subplot(212)
plt.pcolormesh(temps_seg2, freq2, np.abs(Zxx_m1))
plt.axis([0,temps_seg1[-1],0,1300])


#plt.title('Signaux source')
plt.xlabel('Temps (s)')
plt.ylabel('Fréquence (Hz)')
#plt.legend()
plt.show()

#color='#31e016'
#,color='#166fff'
