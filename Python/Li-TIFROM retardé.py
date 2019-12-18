import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.signal import stft
from scipy.signal import square
from scipy.signal import istft
from scipy.signal import correlate
from scipy.io.wavfile import read
from scipy.io.wavfile import write
from contextlib import contextmanager
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import NaN
import cmath as cm
from matplotlib import ticker
from matplotlib.colors import LightSource
from scipy import stats
import matplotlib


## LI-TIFROM retardé (main)
def detecter_sources_visibles_freq(rapp,correlation_max,M,temps_stft,freq_stft,stft_melanges):
    
    var_freq = variance_freq(rapp,M)
    tabvar = np.reshape(var_freq,(len(var_freq)*len(var_freq[0]),3))
    tabvar = tabvar[np.argsort(tabvar[:,2])]#[::-1]
    
    # tabvar de longueur nb_fenetre * nb_tps et contenant les [t,fi,var] trié
    # par variance croissante.
    
    return tabvar
    
    """
    var_freq = correlation(rapp,M,freq_stft)
    tabvar = np.reshape(var_freq,(len(var_freq)*len(var_freq[0]),3))
    tabvar = tabvar[np.argsort(tabvar[:,2])][::-1]
    
    # tabvar de longueur nb_fenetre * nb_tps et contenant les [t,fi,var] trié
    # par variance croissante.
    
    return tabvar
    """
    
    
    ##
    res = []
    indicemax_tabvar = 0 # Indice de la zone limite considérée (var assez faible).
    n = len(tabvar)
    while indicemax_tabvar < (n-1):
        if tabvar[indicemax_tabvar + 1][2] < correlation_max:
            indicemax_tabvar = indicemax_tabvar + 1
        else:
            break
    
    rapports_seuls = np.empty((indicemax_tabvar,2,M), dtype=np.object)
    
    for i in range(0,indicemax_tabvar):
        rapports_seuls[i][0][0] = tabvar[i][0]
        rapports_seuls[i][0][1] = tabvar[i][1]
        
        for j in range(M):
            rapports_seuls[i][1][j] = rapp[tabvar[i][1],tabvar[i][0]]

    return rapports_seuls
  
"""
rapports_seuls[i][0][0] : num tps
rapports_seuls[i][0][1] : num plus petite freq de la fenêtre
rapports_seuls[i][1] : tableau des angles du rapport
"""

def estimer_colonnes_phase(rapports_seuls,temps_stft,freq_stft,N,stft_melanges,seuil):
    indicemax_tabvar = len(rapports_seuls)
    colonnes_estimees_abs = np.zeros((indicemax_tabvar,N-1))
    colonnes_estimees_arg = np.zeros((indicemax_tabvar,N-1))
    
    for i in range(indicemax_tabvar):
        tps = rapports_seuls[i][0][0]
        freq = rapports_seuls[i][0][1]
        pulsations = 2 * np.pi * freq_stft[rapports_seuls[i][0][1]:rapports_seuls[i][0][1] + M] # x
        
        for j in range(1,N):
            rapport_etudie = np.zeros(M)
            
            for k in range(M):
                rapport_etudie[k] = stft_melanges[j][freq + k][tps] / stft_melanges[0][freq + k][tps]
                
            phases_unwrapped = np.unwrap(np.angle(rapport_etudie)) # y
            # Régression linéaire ax + b = y
            A = np.vstack([pulsations, np.ones(len(pulsations))]).T
            a,b = np.linalg.lstsq(A, rapports_seuls[j][2])[0]
            
            alpha = np.mean(np.absolute(rapport_etudie))
            
            colonnes_estimees_abs[i][j] = alpha
            colonnes_estimees_arg[i][j] = -a
    
    
    mat_LI_TIFROM_transp = creer_mat_mixage(N,matrice_rapports,melanges_stft,M)
    
    colonnes_estimees_abs_bon = np.zeros((N,N-1))
    for k in range(N):
        colonnes_estimees_abs_bon[k] = mat_LI_TIFROM_transp[k][1:]
    
    commun = [[] for k in range(N)]
    
    for i in range(indicemax_tabvar):
        colonne_a_mapper = colonnes_estimees_abs[i]
        dist = 0
        num_commun = -1
        trop_loin = False
        for j in range(N):
            
            for k,m in enmuerate(colonnes_estimees_abs_bon[j]):
                if np.abs(colonne_a_mapper[k] - m) > seuil:
                    trop_loin = True
                    break
                    
            if trop_loin == False:
                d_new = numpy.linalg.norm(colonne_a_mapper - colonnes_estimees_abs_bon[j])
                if d_new < dist:
                    num_commun = j
        
        if num_commun != -1:
            commun[num_commun].append(indicemax_tabvar)
    
    if [] in num_commun:
        return False
    
    matrice_retards = np.zeros((N,N))
    for j in range(N):
        for i in range(1,N):
            liste_retards = []
            for k in commun[j]:
                liste_retards.append(colonnes_estimees_arg[k][i])
            
            matrice_retards[i][j] = np.mean(liste_retards)

    stft_sources_retrouves = np.zeros((freq_stft,temps_stft))

    for fk in freq_stft:
        matrice_mixage = np.zeros((N,N))
        
        for j in range(N):
            matrice_mixage[0][j] = 1
    
        for i in range(1,N):
            for j in range(N):
                matrice_mixage[i][j] = colonne_estimees_abs_bon[j][i - 1] * cm.exp(complex(0,fk*matrice_retards[i][j]))
        
        stft_sources_retrouves[fk] = np.dot(np.linalg.inv(matrice_mixage),stft_melanges[:,fk])

    sources_retrouves = istft(stft_sources_retrouves,freq_echant,nperseg = 2)
    #A ECRIRE

    return sources_retrouves




## LI-TIFROM retardé, auxiliaire
def variance_freq(rapp,M): #renvoie une matrice de taille nb_fenetre * nb_tps
                           #donnant la variance pour chaque fenêtre (à temps CSTE).
    nb_freq,nb_tps = rapp.shape
    nb_fenetre = 0
    k = 0
    
    while k < (nb_freq - M): # donne le nb de fenêtre fréquentielles à analyser
        nb_fenetre = nb_fenetre + 1
        k = k + M/2
  
    res = np.zeros((nb_tps,nb_fenetre,3),dtype=object) 
                                        # pos 0 : num temps (=/= temps réel)
                                        # pos 1 : num fenêtre en fréquence
                                        # pos 2 : variance
    
    for t in range(0,nb_tps):
        for fi in range(0,nb_fenetre):
            res[t][fi][2] = calc_variance_freq(rapp,M,fi,t)
            res[t][fi][1] = fi
            res[t][fi][0] = t
    
    return res


def calc_variance_freq(rapp,M,fi,t): # calcul la variance de rapp pour la
                                     # fi-eme fenêtre fréquentielle
                                     # de longueur M, au temps t
    moy = 0
    for j in range(0,M):
        moy = moy + np.abs(rapp[int(fi*M/2 + j)][t])
    moy = 1/M * moy

    var = 0
    
    for j in range(0,M):
        var = var + np.abs(np.abs(rapp[int(fi*M/2 + j)][t]) - moy)**2
    var = 1/M * var
    return var

def calc_moyenne_freq(rapp,M,fi,t):
    moy = 0
    for j in range(0,M):
        moy = moy + np.abs(rapp[int(fi*M/2 + j)][t])
    moy = 1/M * moy
    return moy

def rapport_freq(stft_x1,stft_x2,freq_echant):
    nb_freq,nb_tps = stft_x1.shape
    res = np.zeros((nb_freq,nb_tps),dtype=complex)
    for f in range(0,nb_freq):
        for n in range(0,nb_tps):
            if stft_x2[f][n] != 0:
                res[f][n] = stft_x1[f][n]/stft_x2[f][n]
            else:
                res[f][n] = NaN
    return res

def correlation(rapp,M,freq_stft):
    nb_freq,nb_tps = rapp.shape
    nb_fenetre = 0
    k = 0
    
    while k < (nb_freq - M): # donne le nb de fenêtre fréquentielles à analyser
        nb_fenetre = nb_fenetre + 1
        k = k + M/2
    
    res = np.zeros((nb_tps,nb_fenetre,3),dtype=object)
    
    for t in range(nb_tps):
        for fj in range(nb_fenetre):
            res[t][fj][1] = fj
            res[t][fj][0] = t
            pulsations = 2 * np.pi * freq_stft[fj:fj + M] # x
            rapport_etudie = np.zeros(M,dtype = complex)
            
            for k in range(M):
                rapport_etudie[k] = stft_melanges[1][int(M/2)*fj + k][t] / stft_melanges[0][int(M/2)*fj + k][t]
            phases_unwrapped = np.unwrap(np.angle(rapport_etudie)) # y
            # Régression linéaire ax + b = y
            slope, intercept, r_value, p_value, std_err = stats.linregress(pulsations,phases_unwrapped)
            res[t][fj][2] = np.abs(r_value)
    
    return res
            