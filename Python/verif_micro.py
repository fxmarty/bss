freq_echant, x1 = read("/home/felix/TIPE/Test micro/mic1-vol2-07.wav")
freq_echant, x2 = read("/home/felix/TIPE/Test micro/mic2-vol2-07.wav")

x1pos = x1[x1 > 0]
x2pos = x2[x2 > 0]


##
def meanmax(x,freq_sign,freq_echant):
    nb_echant_fenetre = int(freq_echant/freq_sign)
    nb_fenetre = int(len(x) / nb_echant_fenetre)
    res = []
    for i in range(0,nb_fenetre):
        res.append(max(x[i*nb_echant_fenetre:(i+1)*nb_echant_fenetre]))
    return np.mean(res)

def rapp(x1,x2,freq_sign,freq_echant):
    return meanmax(x1,freq_sign,freq_echant)/meanmax(x2,freq_sign,freq_echant)