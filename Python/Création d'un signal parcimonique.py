# f0 = fréquence fondamentale
# n = numéro d'harmonique
# freq_echant = fréquence d'échantillonnage
# debut = début du signal
# fin = fin du signal
# phase = phase de l'harmonique
# ampli = amplitude de l'harmonique

## Créer signal de sommes de sinusoides choisies "aléatoirement" avec des plats)
import numpy as np
import matplotlib
import IPython.display as ipd
from matplotlib import pyplot as plt
import random
import scipy as sc

##
plt.plot(temps,echant,color='#31e016',label='source2 (s2)')
plt.grid()
plt.show()

##
def creer_harmonique_mult(f0,n,num_debut,num_fin,temps,ampli = 1):
    echantillons = ampli * np.sin(2 * np.pi * n * f0 * temps[num_debut:num_fin] - 2 * np.pi * n * f0 * temps[num_debut])
    return echantillons

# signaux = liste [temps, harmonique1,harmonque2] où les harmoniques sont obtenues
# avec creer_harmonique
def somme_mult(harmonique1,harmonique2):
    res = np.add(harmonique1,harmonique2)
    return res

'''
Effet de bord sur res : sur l'espace de temps entre num_debut et num_fin, crée une somme de sinusoide, de fondamental f0, et d'harmoniques listés dans num_harmoniques.
'''
def creer_sign_mult(f0,freq_echant,num_debut,num_fin,num_harmoniques,temps,res):
    A = res[num_debut:num_fin]
    for k in num_harmoniques:
        harmo = creer_harmonique_mult(f0,k,num_debut,num_fin,temps,1)
        for i in range(len(harmo)):
            res[num_debut + i] = res[num_debut + i] + harmo[i] #+ np.random.normal(0,0.4)



# Renvoie un signal de sommes de sinusoides choisies "aléatoirement" avec nb_plat
# de durée allant de a secondes à b secondes
# où plage_duree_plat = (a,b) !!

# Prérequis : duree est grand devant b



def generer_signal(duree,nb_plat,plage_duree_plat,freq_echant = 44100):
    
    nb_echantillons = int(duree * freq_echant)
    temps = np.linspace(0,duree,nb_echantillons)
    res = np.zeros(nb_echantillons)
    min_duree_plat,max_duree_plat = plage_duree_plat
    
    # On obtient deux array contenant respectivement les numéros d'instants de début de plat
    # et les numéros d'instants de fin de durée de plat
    
    nb_point_max_plat = int(freq_echant * max_duree_plat)
    
    num_instant_debut_plat = np.random.randint(nb_point_max_plat,nb_echantillons - nb_point_max_plat,nb_plat)
    
    num_duree_plat = (freq_echant * np.random.uniform(min_duree_plat,max_duree_plat,nb_plat)).astype(int)
    
    num_instant_fin_plat = np.add(num_instant_debut_plat,num_duree_plat)
    
    # Puis on obtient DANS L'ORDRE les instants de debut et de fin des plats, avec aussi
    # les débuts et fin de l'échantillon
    debut_fin = np.concatenate((num_instant_debut_plat,num_instant_fin_plat,np.array([0,nb_echantillons])))
    debut_fin = np.sort(debut_fin)
    
    # Puis pour chaque période où l'on veut crée une somme de sinusoïde, on la crée grâce
    # aux fonctions auxiliaires précédentes
    harmoniques = []
    for i in range(0,len(debut_fin),2):
        num_deb_harmo = debut_fin[i]
        num_fin_harmo = debut_fin[i + 1] - 1
        duree = temps[num_fin_harmo] - temps[num_deb_harmo]
        num_harmoniques = np.concatenate((np.array([1]),random.sample(range(2,10),2)))
        #num_harmoniques = np.array([1])
        harmoniques.append(num_harmoniques)
        
        creer_sign_mult(1/duree,freq_echant,num_deb_harmo,num_fin_harmo,num_harmoniques,temps,res)
        #creer_sign_mult(20,freq_echant,num_deb_harmo,num_fin_harmo,num_harmoniques,temps,res)
    
    
    plt.plot(temps,res)
    plt.grid()
    plt.show()
    return temps,res,harmoniques,debut_fin

##
def CL(coeff1,signal1,coeff2,signal2): #Effectue le mélange
    res = np.add(coeff1 * signal1,coeff2 * signal2)
    return res

def rapportbis(melange1,melange2):
    long = len(melange1)
    res = np.zeros(long)
    for i in range(long):
        if melange2[i] != 0:
            res[i] = melange1[i]/melange2[i]
    return res

## Load exemples
temps = np.load('/home/felix/TIPE/Expériences/Expérience théorique diapo/time.npy')
source1 = np.load('/home/felix/TIPE/Expériences/Expérience théorique diapo/test1_source1.npy')
source2 = np.load('/home/felix/TIPE/Expériences/Expérience théorique diapo/test1_source2.npy')

#source1 = source1 + np.random.normal(0, 0.001, size=441000)
#source2 = source2 + np.random.normal(0, 0.001, size=441000)


##
x1 = CL(0.9,source1,0.6,source2)
x2 = CL(0.7,source1,0.9,source2)
rapport_m1m2 = rapportbis(x1,x2)

##
freq_echant,source1 = read('/home/felix/echant_cc3.wav')
nb_echantillons = len(source1)-300

nouv_tps = np.linspace(0,nb_echantillons/freq_echant,nb_echantillons)

source1_bis = source1[300:] + np.random.normal(0, 0.005, size=nb_echantillons)
source2_bis = 0.7*source1[:-300] + + np.random.normal(0, 0.005, size=nb_echantillons)


## Sources
plt.figure(1)
plt.subplot(211)
plt.plot(nouv_tps,source1_bis,label='s1',color='royalblue')
plt.legend()
plt.grid()
#plt.ylabel('Signal (V)')
plt.ylim(-0.56,0.35)

"""
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
"""

plt.subplot(212)
plt.plot(nouv_tps,source2_bis,label='s1 atténué et retardé',color='royalblue')

#plt.title('Signaux source')
plt.xlabel('Temps (s)')
#plt.ylabel('Signal (V)')
plt.legend()
plt.grid()
plt.ylim(-0.56,0.35)
plt.show()

#color='#31e016'
#,color='#166fff'

## Mélanges et rapport
#plt.plot(temps,rapport_m1m2,linestyle = 'None', marker='.',markersize=1,markevery=500,color='#f9331d',label='Rapport melange1/melange2')
#plt.plot(temps,x1,color='#45a545',label='melange1 (x1)',markevery=500)
#plt.plot(temps,x2,color='#456aa5',label='melange2 (x2)',markevery=500)
plt.plot(temps,source1,label='source1 (s1)')
plt.plot(temps,source2,label='source2 (s2)')
#plt.plot(temps,rapptps,linestyle = 'None', marker='.',markersize=3,markevery=50)
#plt.plot(x1,x2,linestyle = 'None', marker='.',markersize=5,markevery=200,label='Mode X-Y à tout instant')
#plt.plot(x1_res,x2_res,linestyle = 'None', marker='.',markersize=4,markevery=10,color = 'r',label='Mode X-Y au niveau des zones monosources détectées')
plt.plot(temps,res[0],label='Estimation 1')
plt.plot(temps,res[1],label='Estimation 2')

plt.xlabel('x1',fontsize=16)
plt.ylabel('x2',rotation=0,fontsize=16)
plt.legend(fontsize=14)
#plt.title('Mélange linéaire des sources et leur rapport')
plt.grid()
plt.show()

## Préparation du résultat
Mat_melange = np.array([melange2,melange1])
rapport1 = rapport_m1m2[44100 * 2]
rapport2 = rapport_m1m2[44100 * 3]
Mat_mixage = np.array([[1,1],[rapport1,rapport2]])
res = np.dot(np.linalg.inv(Mat_mixage),Mat_melange)

## Résultat
plt.plot(temps,source1,color='#0ed612',label = 'source1 (s1)')
plt.plot(temps,source2,color='#369638',label = 'source2 (s2)')
plt.plot(temps,res[1], label = 'résultat1 (y1)',color='r')
plt.plot(temps,res[0], label = 'résultat2 (y2)',color='b')
plt.xlabel('temps (s)')
plt.ylabel('Sources et mélanges séparés (V)')
plt.title('Comparasion entre les sources et leur évaluation')

plt.legend()
plt.grid()
plt.show()

##
sr = 44100 # sample rate
T = 2.0    # seconds
t = np.linspace(0, T, int(T*sr), endpoint=False) # time variable
x = 0.5*np.sin(2*np.pi*440*t)                # pure sine wave at 440 Hz

##
framerate = 44100
t = np.linspace(0,5,framerate*5)
data = np.sin(2*np.pi*220*t) + np.sin(2*np.pi*224*t)

