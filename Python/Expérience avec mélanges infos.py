##
x1,x2 = melange_freq(s1,s2,0.6,0.9,0.8,0.5)
##
stft_melanges = np.array([Zxx_x1,Zxx_x2])
melanges = np.array([x1,x2])

rapp_freq = rapport_freq(stft_melanges[1],stft_melanges[0],freq_echant)

rapports_seuls = detecter_sources_visibles_freq(rapp_freq,1e-4,20,temps_seg1,freq1,stft_melanges)

#rappseul[rappseul[:,2].argsort()]

##
M = 20
L = []
for i,k in enumerate(rapports_seuls[0:20]):
    fj = k[1]
    t = k[0]
    #fj = 189
    #t = 244
    pulsations = freq1[int(M/2)*fj:int(M/2)*fj + M] # x
    rapport_etudie = np.zeros(M,dtype = complex)
            
    for p in range(M):
        rapport_etudie[p] = stft_melanges[1][int(M/2)*fj + p][t] / stft_melanges[0][int(M/2)*fj + p][t]
    
    phases_unwrapped = np.unwrap(np.angle(rapport_etudie)) # y
    # Régression linéaire ax + b = y
    slope, intercept, r_value, p_value, std_err = stats.linregress(pulsations,phases_unwrapped)
    
    if np.abs(r_value) > 0.9:
        s = 0
        for j in range(M):
            s = s + np.abs(rapp_freq[k[1] * int(M/2) + j][k[0]])
        s = s/M
        L.append(slope)
        print(i,np.abs(r_value),s,slope)

print(np.mean(L))

##
def wav_sources(mat_sources_separees,freq_echant):
    N = len(mat_sources_separees)
    for i in range(0,N):
        write("/home/felix/sourceseparee%d.wav" %i,freq_echant,
        mat_sources_separees[i])

wav_sources(mat_sourcesextraites,freq_echant)
##
L = rapports_seuls[0:20]
bon = []
M = 20
for i,k in enumerate(L):
    ra = 0
    for j in range(M):
        ra = ra + np.abs(rapp_freq[k[1] * int(M/2) + j][k[0]])
    ra = ra/M
    if ra > 1: #and len(bon) < 20:
        print(i,ra)
        bon.append(i)


## Mélanges non retardés
freq1,temps_seg1,Zxx_x1 = stft(x1,freq_echant,nperseg = 2000)
freq2,temps_seg2,Zxx_x2 = stft(x2,freq_echant,nperseg = 2000)
stft_x1 = np.abs(Zxx_x1)
stft_x2 = np.abs(Zxx_x2)
melanges_stft = np.array([stft_x1,stft_x2])
melanges = np.array([x1,x2])

rapp = rapport(stft_x1,stft_x2,freq_echant)
var,nb_fen = variance(rapp,10)

var_trie = trier(var)
visi,indmax,rappseul = detecter_sources_visibles(rapp,var_trie,5e-3,0.1,10)
mat_mixage = creer_mat_mixage(2,visi,melanges_stft,14)
mat_sourcesextraites = retrouver_sources(mat_mixage,melanges)

#rappseul[rappseul[:,2].argsort()]

## Fake retardé

stft_x1 = np.abs(stft_melanges[0])
stft_x2 = np.abs(stft_melanges[1])
melanges_stft = np.array([stft_x1,stft_x2])
#melanges = np.array([x1,x2])

rapp = rapport(stft_x2,stft_x1,freq_echant)
var,nb_fen = variance(rapp,10)

var_trie = trier(var)
visi,indmax,rappseul = detecter_sources_visibles(rapp,var_trie,1e-3,0.1,10)


#mat_mixage = creer_mat_mixage(2,visi,melanges_stft,14)
mat_mixagebis = np.array([[1,1],[1.33333*np.exp(-complex(0,10)),0.55555]])

stft_res = np.zeros((2,len(freq1),len(temps_seg1)),dtype=complex)
for fk in range(len(freq1)):
    mat_mixage_inv = np.linalg.inv(np.array([[1,1],[1.33333*np.exp(-complex(0,10*2*np.pi*freq1[fk])),0.55555]]))
    L = np.dot(mat_mixage_inv,stft_melanges[:,fk])
    stft_res[0][fk] = L[0]
    stft_res[1][fk] = L[1]


tps1,res1 = istft(stft_res[0],freq_echant,nperseg = 4000)
tps2,res2 = istft(stft_res[1],freq_echant,nperseg = 4000)


mat_sourcesextraites = np.array([res1,res2])
wav_sources(mat_sourcesextraites,freq_echant)

#rappseul[rappseul[:,2].argsort()]

## Fake retardé
freq_echant, s1 = read("/home/felix/TIPE/Expériences/Flamenco/Flamenco1_U89_cut.wav")
freq_echant, s2 = read("/home/felix/TIPE/Expériences/Flamenco/Flamenco2_U89_cut.wav")

freq1,temps_seg1,Zxx_s1 = stft(s1,freq_echant,nperseg = 4000)
freq2,temps_seg2,Zxx_s2 = stft(s2,freq_echant,nperseg = 4000)
Mat_sources = np.array([Zxx_s1,Zxx_s2])


stft_melanges = np.zeros((2,len(freq1),len(temps_seg1)),dtype=complex)
for fk in range(len(freq1)):
    arg = - freq1[fk]
    Mat_mixage = np.array([[1,0.5*np.exp(complex(0,arg*20/freq_echant))],[0.5*np.exp(complex(0,arg*20/freq_echant)),1]])

    L = np.dot(Mat_mixage,Mat_sources[:,fk])
    stft_melanges[0][fk] = L[0]
    stft_melanges[1][fk] = L[1]

rapp_freq = rapport_freq(stft_melanges[1],stft_melanges[0],freq_echant)

rapports_seuls = detecter_sources_visibles_freq(rapp,1e-4,20,temps_seg1,freq1,stft_melanges)



## Fake retardé
freq_echant, s1 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-01.wav")
freq_echant, s2 = read("/home/felix/Simulation salle/Test LI-TIFROM retardé/Distance 0.5 bis/m-02.wav")

#freq_echant, s1 = read("/home/felix/Flamenco1_U89_cut.wav")
#freq_echant, s2 = read("/home/felix/Flamenco2_U89_cut.wav")

x1,x2 = s1,s2

freq1,temps_seg1,Zxx_x1 = stft(x1,freq_echant,nperseg = 4000)
freq2,temps_seg2,Zxx_x2 = stft(x2,freq_echant,nperseg = 4000)
stft_melanges = np.array([Zxx_x1,Zxx_x2])

rapp_freq = rapport_freq(stft_melanges[1],stft_melanges[0],freq_echant)

rapports_seuls = detecter_sources_visibles_freq(rapp,1e-4,20,temps_seg1,freq1,stft_melanges)

##
M = 20
L = []
diff = []
for i,k in enumerate(rapports_seuls[0:20]):
    fj = k[1]
    t = k[0]
    pulsations = freq1[int(M/2)*fj:int(M/2)*fj + M] # x
    rapport_etudie = np.zeros(M,dtype = complex)
            
    for p in range(M):
        rapport_etudie[p] = stft_melanges[1][int(M/2)*fj + p][t] / stft_melanges[0][int(M/2)*fj + p][t]
    
    phases_unwrapped = np.unwrap(np.angle(rapport_etudie)) # y
    # Régression linéaire ax + b = y
    slope, intercept, r_value, p_value, std_err = stats.linregress(pulsations,phases_unwrapped)
    
    if np.abs(r_value) > -1:
        s = 0
        for j in range(M):
            s = s + np.abs(rapp_freq[k[1] * int(M/2) + j][k[0]])
        s = s/M
        L.append(slope)
        print(i,np.abs(r_value),s,slope)

print(np.mean(L))
##
L = rapports_seuls[0:50]
M = 20
for i,k in enumerate(L):
    ra = 0
    for j in range(M):
        ra = ra + np.abs(rapp_freq[k[1] * int(M/2) + j][k[0]])
    ra = ra/M
    #if ra < 1:
    print(i,ra,k[0],k[1])