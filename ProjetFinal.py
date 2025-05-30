import math
import time
import matplotlib.pyplot as plt
import numpy as np
from memory_profiler import memory_usage



# CREATION MOT
def creer_tableau(fichier):
    x=[]
    y=[]
    f=open(fichier,'r')

    x_size = int(f.readline().strip())     #.strip() retire les '\n'
    y_size = int(f.readline().strip())

    ligne_x=f.readline().strip()
    for c in ligne_x:
        if c!=' ':
            x.append(c)
    ligne_y=f.readline().strip()
    for c in ligne_y:
        if c!=' ':
            y.append(c)
            
    f.close()
    return ((x_size,x),(y_size,y))


# COUTS OPERATIONS 
c_del = 2
c_ins = 2
def c_sub(a,b):
    if(a==b):
        return 0
    elif a == "A" and b == "T" or a == "T" and b == "A":
        return 3
    elif a == "C" and b == "G" or a == "G" and b == "C":
        return 3
    else:
        return 4


# ===================== TACHE A ===================== #
res1=[]
def dist_naif(x,y):
    st = time.time()
    res = dist_naif_rec(x,y,0,0,0,math.inf)
    et = time.time()
    res1.append(et-st)
    return res


def dist_naif_rec(x,y,i,j,c,dist):
    if i==len(x) and j==len(y):
        if c<dist:
            dist=c
    else:
        if i<len(x) and j<len(y):
            dist=dist_naif_rec(x,y,i+1,j+1,c+c_sub(x[i],y[j]),dist)
        if i<len(x):
            dist=dist_naif_rec(x,y,i+1,j,c+c_del,dist)
        
        if j<len(y):
            dist=dist_naif_rec(x,y,i,j+1,c+c_ins,dist)

    return dist

"""
# TESTS TACHE A
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_7.adn")
print(dist_naif(x,y))
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_8.adn")
print(dist_naif(x,y)) 
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_44.adn")
print(dist_naif(x,y)) 
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000012_13.adn")
print(dist_naif(x,y)) 
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000012_32.adn")
print(dist_naif(x,y)) 
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000012_56.adn")
print(dist_naif(x,y)) 
"""

# ===================== TACHE B ===================== #

def DIST_1(x,y): 
    n=len(x)
    m=len(y)
    #Initialisation de la matrice
    T= [[0 for i in range(m+1)] for j in range(n+1)]
    for i in range(n+1):
        T[i][0]=i*c_del
    for j in range(1,m+1):
        T[0][j]=j*c_ins

    for i in range(1,n+1):
        for j in range(1,m+1):
            T[i][j]= min(T[i-1][j]+c_del,T[i][j-1]+c_ins,T[i-1][j-1]+c_sub(x[i-1],y[j-1]))

    return (T,T[n][m])


def SOL_1(x,y,T):
    u=[]
    v=[]
    i=len(x)
    j=len(y)

    while i>0 and j>0:
        if T[i][j]==T[i][j-1]+c_ins:
            u.insert(0,'-')
            v.insert(0,y[j-1])
            j=j-1
            continue
        if T[i][j]==T[i-1][j]+c_del:
            u.insert(0,x[i-1])
            v.insert(0,'-')
            i=i-1
            continue
        if T[i][j]==T[i-1][j-1]+c_sub(x[i-1],y[j-1]):
            u.insert(0,x[i-1])
            v.insert(0,y[j-1])
            i=i-1
            j=j-1
            continue

    while i>0:
        u.insert(0,x[i-1])
        v.insert(0,'-')
        i=i-1
    while j>0:
        u.insert(0,'-')
        v.insert(0,y[j-1])
        j=j-1
    
    return (u,v)


def PROG_DYN(x,y):
    return (DIST_1(x,y)[1],SOL_1(x,y,DIST_1(x,y)[0]))


'''
# TESTS TACHE B
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_7.adn")
print(PROG_DYN(x,y))
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_44.adn")
print(PROG_DYN(x,y))
'''


# ===================== TACHE C ===================== #

def DIST_2(x,y):
    n=len(x)
    m=len(y)
    T= [[0 for i in range(m+1)] for j in range(2)]

    cpt=0
    i=1
    for j in range(m+1):
        T[0][j]=j*c_ins
    T[1][0]= c_del

    while cpt<n:
        for j in range(1,m+1):
            T[i][j]=min(T[(i-1)%2][j]+c_del,T[i][j-1]+c_ins,
                        T[(i-1)%2][j-1]+c_sub(x[cpt],y[j-1]))
        i=(i+1)%2
        T[i][0] = T[(i-1)%2][0] +2
        cpt+=1

    if n%2==0 :
        return (T,T[0][m])
    else : 
        return (T,T[1][m])

'''
# TESTS TACHE C
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_7.adn")
print(DIST_2(x,y)[1])
((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000010_44.adn")
print(DIST_2(x,y)[1])
'''


# ===================== TACHE D ===================== #


def mot_gaps(k):
    return k*'-'


def align_lettre_mot(x,y):
    m=len(y)
    i = 0
    while i < m:
        if y[i] == x :
            return (mot_gaps(i)+x+mot_gaps(m-i-1) , y)  # Cas trivial
        i += 1

    i = 0
    while i < m :
        if c_sub(x, y[i]) == 3 :
            return (mot_gaps(i)+x+mot_gaps(m-i-1) , y)
        i += 1
    if i==m :
        return (mot_gaps(i-1)+str(x) , y)

def coupure(x,y):
    n = len(x)
    m = len(y)
    # Initialisation de deux tableaux T et I de dimensions 2*m
    T= [[0 for i in range(m+1)] for j in range(2)]
    I= [[0 for i in range(m+1)] for j in range(2)]
    coupe = n//2
    cmp = 0 
    i1 = 1
    i2 = 1 
    for j in range(m+1) :
        T[0][j] = j*c_ins
    T[1][0] = c_del #insertion dans y
    for j in range(m+1): 
        I[0][j] = j 
        I[1][j] = j  

    while cmp < n :
        for j in range(1,m+1) :
            T[i1][j] = min(T[(i1-1)%2][j]+c_del , T[i1][j-1]+c_ins,T[(i1-1)%2][j-1]+ c_sub(x[cmp],y[j-1]))
            if cmp >= coupe :
                if T[i1][j]==T[(i1-1)%2][j]+c_del:
                    I[i2][j] = I[(i2-1)%2][j]
                if T[i1][j]==T[i1][j-1]+c_ins :
                    I[i2][j] = I[i2][j-1]
                if T[i1][j] == T[(i1-1)%2][j-1]+c_sub(x[cmp], y[j-1]) :
                    I[i2][j] = I[(i2-1)%2][j-1]

        if cmp >= coupe :
            i2 = (i2+1)%2
        i1 = (i1+1)%2 
        T[i1][0] = T[(i1-1)%(2)][0]+ 2
        cmp += 1 

    if (n-coupe)%2 == 0 :
        return I[0][m]
    return I[1][m]



def SOL_2(x,y,u,v):
    n = len(x)
    m = len(y)
    i = n//2
    if m == 0 and n != 0  :
        u += x
        v += mot_gaps(n)
    if n == 0 and m != 0  :
        u += mot_gaps(m)
        v += y
    if n == 1 and m >= 1 :
        u += align_lettre_mot(x,y)[0]
        v += align_lettre_mot(x,y)[1]
    if n > 1 and m >= 1 :
        j = coupure(x,y)
        SOL_2(x[:i],y[:j],u,v)
        SOL_2(x[i:],y[j:],u,v)
    return (u,v)


# TESTS TACHE D
'''
print(mot_gaps(5))
print(align_lettre_mot('A',"CCCGTAG"))
print(align_lettre_mot('A',"CCTGCCG"))
'''
#print(SOL_2("ATTGTA","ATCTTA",[],[]))

((n,x),(m,y))=creer_tableau("Instances_genome/Inst_0000050_2.adn")
x2=""
for temp in x:
    x2+=(temp)

y2=""
for temp in y:
    y2+=(temp)

#assert(PROG_DYN(x,y)==SOL_2(x,y,[],[]))
#print(PROG_DYN(x,y)[1])
#print('\n')
#print()
print(SOL_2(x2,y2,[],[])[0])
print(SOL_2(x2,y2,[],[])[1])

