import math
from operator import itemgetter

def calc_moyenneTEMP(rapp,M,i):
    moy = 0
    for j in range(0,M):
        moy = moy + rapp[int(i*M/2 + j)]
    moy = 1/M * moy
    return moy

def calc_varianceTEMP(rapp,M,i): # calcul la variance de rapp pour la ieme fenêtre
                               # de longueur M, à la fréquence num f
    moy = 0
    for j in range(0,M):
        moy = moy + rapp[int(i*M/2 + j)]
    moy = 1/M * moy

    var = 0
    
    for j in range(0,M):
        var = var + abs(rapp[int(i*M/2 + j)] - moy)**2
    var = 1/M * var
    return var

def varianceTEMP(rapp,M):
    k = 0
    nb_tps = len(rapp)
    nb_fenetre = 0
    
    while k < (nb_tps - M): # donne le nb de fenêtre temporelles à analyser
        nb_fenetre = nb_fenetre + 1
        k = k + M/2
    
    res = np.zeros((nb_fenetre,2)) 
                                        # pos 0 : num fenêtre
                                        # pos 1 : variance
    
    for i in range(0,nb_fenetre):
            res[i][1] = calc_varianceTEMP(rapp,M,i)
            res[i][0] = i
    
    return res,nb_fenetre

def rapportTEMP(melange1,melange2):
    long = len(melange1)
    res = np.zeros(long)
    for i in range(long):
        res[i] = NaN
        if melange2[i] != 0:
            res[i] = melange1[i]/melange2[i]
    return res
##

rapp = rapportTEMP(x1,x2)

res,nb_fenetre = varianceTEMP(rapp,8)

def supprimer_nan(rapp_trie):
    n = len(rapp_trie)
    rapp_trie = rapp_trie.tolist()
    res = []
    for i in range(n):
        if not math.isnan(rapp_trie[i][1]):
            res.append(rapp_trie[i])
        
    return res



def construire_rapptps(variances_fenetres,M,rapport):
    pas = int(M/2)
    n = len(rapport)
    res = np.zeros(n)
    nb_fenetre = len(variances_fenetres)
    for i in range(nb_fenetre):
        for k in range(pas):
                res[i * pas + k] = 1/(variances_fenetres[i][1] * 1e8)
    return res

rapptps = construire_rapptps(res,8,rapp)

##
rapp = supprimer_nan(res)
rapp_trie = sorted(rapp, key=itemgetter(1))
    
## Nuage de points bons
def creernuage(rapp_trie,var_max,x1,x2,M):
    x1_res = []
    x2_res = []
    i = 0
    n = len(rapp_trie)
    while i < n and rapp_trie[i][1] < var_max:
        x1_res.append(x1[int(rapp_trie[i][0] * M/2)])
        x2_res.append(x2[int(rapp_trie[i][0] * M/2)])
        i = i + 1
        
    return x1_res,x2_res

x1_res,x2_res = creernuage(rapp_trie,1e-14,x1,x2,8)

## Retrouver
res = np.dot(np.linalg.inv(np.array([[1,1],[c1,c2]])),[x1,x2])
