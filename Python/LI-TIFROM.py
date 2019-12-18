import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.signal import stft
from scipy.signal import square
from scipy.io.wavfile import read
from scipy.io.wavfile import write
from contextlib import contextmanager
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import NaN
import math

from matplotlib import ticker
from matplotlib.colors import LightSource
##
@contextmanager
def show_complete_array():
    oldoptions = np.get_printoptions()
    np.set_printoptions(threshold=np.inf)
    yield
    np.set_printoptions(**oldoptions)

def afficher(arai,nb):
    for i in range(nb):
        print(arai[i])
    
## Importation des deux mélanges wav vers deux array x1 et x2
freq_echant, x1 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-01.wav")
freq_echant, x2 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-02.wav")
temps = np.linspace(0,nb_echantillons/freq_echant,nb_echantillons)

##
freq_echant, s1 = read("/home/felix/Flamenco1_U89_cut.wav")
freq_echant, s2 = read("/home/felix/Flamenco2_U89_cut.wav")
nb_echantillons = len(s1)
temps = np.linspace(0,nb_echantillons/freq_echant,nb_echantillons)

x1,x2 = melange(s1,s2,0.6,0.9,0.8,0.5)
write("/home/felix/flamencomelange1.wav",48000,x1)
write("/home/felix/flamencomelange2.wav",48000,x2)

freq1,temps_seg1,Zxx_x1 = stft(x1,freq_echant)
freq2,temps_seg2,Zxx_x2 = stft(x2,freq_echant)
stft_x1 = np.abs(Zxx_x1)
stft_x2 = np.abs(Zxx_x2)
melanges_stft = np.array([stft_x1,stft_x2])
melanges = np.array([x1,x2])

rapp = rapport(stft_x1,stft_x2,freq_echant)
var = variance(rapp,8)

var_trie = trier(var)
visi,indmax,rappseul = detecter_sources_visibles(rapp,var_trie,1e-5,0.1,8)
mat_mixage = creer_mat_mixage(2,visi,melanges_stft,8)
mat_sourcesextraites = retrouver_sources(mat_mixage,melanges)


## Mélange de sources
def melange(s1,s2,coef11,coef12,coef21,coef22):
    x1 = coef11 * s1 + coef12 * s2
    x2 = coef21 * s1 + coef22 * s2
    return x1,x2

## Juste pour les voir
plt.plot(temps,x1,linewidth=0.8)
plt.plot(temps,x2,linewidth=0.8)
plt.show()

## Calculer stft des deux mélanges
freq1,temps_seg1,Zxx_x1 = stft(x1,freq_echant)
freq2,temps_seg2,Zxx_x2 = stft(x2,freq_echant)
stft_x1 = np.abs(Zxx_x1)
stft_x2 = np.abs(Zxx_x2)

## Juste pour voir les spectrogrammes
plt.pcolormesh(temps_seg1, freq1, stft_x1) #para vmin,vmax
#plt.pcolormesh(temps_seg2, freq2, stft_x2)
plt.title('STFT Magnitude')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

##
# Calculer rapport de deux stft de mélanges
def rapport(stft_x1,stft_x2,freq_echant):
    nb_freq,nb_tps = stft_x1.shape
    res = np.zeros((nb_freq,nb_tps))
    for f in range(0,nb_freq):
        for n in range(0,nb_tps):
            if stft_x2[f][n] != 0:
                res[f][n] = stft_x1[f][n]/stft_x2[f][n]
            else:
                res[f][n] = NaN
    return res

#Calculer variance du rapport DE DEUX MELANGES SEULMENT pour chaque fréquence

# fenêtre temporelle de longueur M, supposé pair

def calc_moyenne(rapp,M,f,i):
    moy = 0
    for j in range(0,M):
        moy = moy + rapp[f][int(i*M/2 + j)]
    moy = 1/M * moy
    return moy

def calc_variance(rapp,M,f,i): # calcul la variance de rapp pour la ieme fenêtre
                               # de longueur M, à la fréquence num f
    moy = 0
    for j in range(0,M):
        moy = moy + rapp[f][int(i*M/2 + j)]
    moy = 1/M * moy

    var = 0
    
    for j in range(0,M):
        var = var + np.abs(rapp[f][int(i*M/2 + j)] - moy)**2
    var = 1/M * var
    return var

# i : INDIQUE LE NUMERO DE LA FENETRE

def variance(rapp,M): #renvoie une matrice de taille nb_freq * nb_fenetre
                      #donnant la variance pour chaque fenêtre,fréquence.
    nb_freq,nb_tps = rapp.shape
    nb_fenetre = 0
    k = 0
    
    while k < (nb_tps - M): # donne le nb de fenêtre temporelles à analyser
        nb_fenetre = nb_fenetre + 1
        k = k + M/2
  
    res = np.zeros((nb_freq,nb_fenetre,3),dtype=object) 
                                        # pos 0 : num fréquence (=/= freq réelle)
                                        # pos 1 : num fenêtre
                                        # pos 2 : variance
    
    for i in range(0,nb_fenetre):
        for f in range(0,nb_freq):
            res[f][i][2] = calc_variance(rapp,M,f,i)
            res[f][i][0] = f
            res[f][i][1] = i
    
    return res,nb_fenetre

# Tri 

def trier(A): # A = matrice de taille nb_freq * nb_fen * 3
              # Renvoie une liste de longueur nb_freq * nb_fen triée selon la
              # variance croissante
    atrier = np.reshape(A,(len(A)*len(A[0]),3))
    atrier = atrier[np.argsort(atrier[:,2])]
    return atrier # atrier = tableau des variances triées

# Détection des sources visibles, détection des zones mono-sources

# rapp : rapport de deux mélanges
# tabvar : tableau des variances croissantes associé
def detecter_sources_visibles(rapp,tabvar,variance_max,seuil,M,N):
    res = []
    indicemax_tabvar = 0 # Indice de la zone limite considérée (var assez faible).
    n = len(tabvar)
    while indicemax_tabvar < (n-1):
        if tabvar[indicemax_tabvar + 1][2] < variance_max:
            indicemax_tabvar = indicemax_tabvar + 1
        else:
            break
    
    rapports_seuls = np.copy(tabvar[0:indicemax_tabvar + 1])
    num_source = np.zeros(indicemax_tabvar + 1,dtype = int)
            # Contient pour chaque rapport le numéro de la source associé
            # (complété au fur et à mesure après)
    
    
    # On remplace la variance par les rapports dans rapports_seuls
    for i in range(0,indicemax_tabvar + 1):
        rapports_seuls[i][2] = \
                calc_moyenne(rapp,M,rapports_seuls[i][0],rapports_seuls[i][1])
                
    res.append([[rapports_seuls[0][0],rapports_seuls[0][1]]]) #Initialisation de res
    
    for i in range(1,indicemax_tabvar + 1):
        indice_memesource = -1 # Indice où la source seule est identique.
        for k in range(0,i):
            if np.abs(rapports_seuls[i][2] - rapports_seuls[k][2]) < seuil:
                indice_memesource = k
                break
                
        if indice_memesource != -1:
            if len(res[num_source[indice_memesource]]) < 15:
                res[num_source[indice_memesource]].append([rapports_seuls[i][0],rapports_seuls[i][1]])
                num_source[i] = num_source[indice_memesource]
            else:
                rapports_seuls[i][2] = NaN
        
        else:
            if len(res) < N:
                res.append([[rapports_seuls[i][0],rapports_seuls[i][1]]])
                num_source[i] = len(res) - 1
            else:
                rapports_seuls[i][2] = NaN

    nb_rapports = 0
    for k in res:
        nb_rapports = nb_rapports + len(k)
    
    rapports_seuls_new = np.empty((nb_rapports,3), dtype=np.object)
    count = 0
    for k in rapports_seuls:
        if not math.isnan(k[2]):
            rapports_seuls_new[count] = k
            count = count + 1

    return res,indicemax_tabvar,rapports_seuls_new


"""
        if indice_memesource != -1:
            if len(res[num_source[indice_memesource]]) < 20:
                res[num_source[indice_memesource]].append([rapports_seuls[i][0],rapports_seuls[i][1]])
                num_source[i] = num_source[indice_memesource]
        
        else:
                res.append([[rapports_seuls[i][0],rapports_seuls[i][1]]])
                num_source[i] = len(res) - 1

    return res,indicemax_tabvar,rapports_seuls
"""

    # res : matrice de taille nb_sources_vis * ? * 2
    #       donnant toutes les zones mono-sources par source.
    #       Le "?" change selon la source, mais minimum 1.
    #       - Colonne 0 : num fréquence.
    #       - Colonne 1 : num fenêtre temporelle.


# Calcul moyenne du rapport entre deux mélanges en une zone monosource

#i : num de la fenêtre temporelle
def calc_moy_monosource(stft_x1,stft_x2,M,f,i):
    moy = 0
    for k in range(0,M):
        moy = moy + stft_x1[f][int(M/2 * i + k)]/stft_x2[f][int(M/2 * i + k)]
    moy = 1/M * moy
    return moy



# Cas N = Q = P ***********************

# où N = nb de sources ; Q = nb de sources visibles ; P = nb de mélanges.

#matrice_rapports = le res du résultat précédent
#melanges = matrice contenant tous les STFT des mélanges 0,..,(N-1).
def creer_mat_mixage(N,matrice_rapports,rapp,M):
    if N != len(matrice_rapports):
        return False
    
    else:
        mat_mixage = np.zeros((N,N))
    
        for j in range(0,N):        # 1ere ligne
            mat_mixage[0][j] = 1
    
        for i in range(1,N):
            for j in range(0,N):    # Case i,j
                nb_cij = 0
                for L in matrice_rapports[j]:
                    nb_cij = nb_cij + np.mean(rapp[L[0]][int(M/2)*L[1]:int(M/2)*L[1]+M])
                
                mat_mixage[i][j] = nb_cij/len(matrice_rapports[j])
        return mat_mixage

def retrouver_sources(mat_mixage,melanges):
    return np.dot(np.linalg.inv(mat_mixage),melanges)


##
def wav_sourcesbis(mat_sources_separees,freq_echant):
    N = len(mat_sources_separees)
    for i in range(0,N):
        write("/home/felix/sourceseparee%d.wav" %i,freq_echant,
        mat_sources_separees[i])