import numpy as np
import matplotlib.pyplot as plt

## Aim: create a single signal composed of a sum of weighted sinus

# f0 = fundamental frequency of the signal
# n = harmonic order (starting from 1)
# sampl_f = sampling frequency
# start_t = starting time
# end_t = ending time
# ampli = amplitude of the harmonic

# createHarmonic: Creates one harmonic of frequency n*f0 with default amplitude 1
def createHarmonic(f0,n,sampl_f,start_t,end_t,ampli = 1):
    duration = end_t - start_t
    nb_samples = int(sampl_f*duration)
    time = np.linspace(start_t,end_t,nb_samples)
    samples = ampli * np.sin(2 * np.pi * n * f0 * time)
    return samples
    
# num_harmonics = List of the harmonics orders to include in the signal to create

# createSignal: Creates a signal of fundamental frequency f0, with its harmonics 
# num_harmonics associated to the amplitudes ampli
def createSignal(f0,sampl_f,start_t,end_t,num_harmonics,ampli):
    duration = end_t - start_t
    nb_samples = int(sampl_f*duration)
    res = np.zeros(nb_samples)
    time = np.linspace(start_t,end_t,nb_samples)
    for i,k in enumerate(num_harmonics):
        harmo = createHarmonic(f0,k,sampl_f,start_t,end_t,ampli[i])
        res = np.add(res,harmo)
    maxi = np.max(np.abs(res))
    res = res / (maxi + 0.1)
    return time,res

# visualiseSignal: Show with matplotlib the signal created using createSignal
def visualiseSignal(f0,sampl_f,start_t,end_t,harmoniques,ampli):
    time,res = createSignal(f0,sampl_f,start_t,end_t,harmoniques,ampli)
    plt.plot(time,res)
    plt.grid()
    plt.show()