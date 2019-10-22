# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:09:18 2019

@author: hindesa
"""
import numpy as np
from numpy import random
from numpy import linalg as LA
import matplotlib.pyplot as plt
from tqdm import tqdm

from math import isnan
from jacobian import computeJac
from systemTypes import typeNode

# M = mangrove biomass (kg), P = peat soil elevation (mm)
# S = soil salinity (ppm?

# Linspace range
n = 100000

# ---------------
# Timescales
# ---------------
#years
alphaM0 = 3/12
alphaM1 = 5

alphaP0 = 6/12
alphaP1 = 10

alphaS0 = 6/12
alphaS1 = 10

mangLifespan= 4 # average leaf lifespan #https://academic.oup.com/treephys/article/30/9/1148/1641261
unitPeatTime = 1/2 # how long does a chunk of peat topsoil last? https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2010WR009492
unitSaltTime = 0.5 # timescale for unit of salt https://link.springer.com/article/10.1023/A:1008470719913

alphaM = random.uniform(alphaM0,alphaM1,n)
alphaP = random.uniform(alphaP0,alphaP1,n)
alphaS = random.uniform(alphaS0,alphaS1,n)

# ----------------
# Beta parameters
# ----------------
# Define function to get 3 random points betwee 0,1, such that
# b1 + b2 + b3 =1
def pick3():
    d = random.uniform(0,1)
    c = random.uniform(0+d/2.,1-d/2.)
    b1 = c-d/2.
    b2 = d
    b3 = 1-b1-b2
    picks = [b1,b2,b3]
    random.shuffle(picks)
    return picks

# Must add to 1
# Mangrove gain
betaP = random.uniform(0,1,n) # from propagules & established saplings
betaG = 1-betaP # from endogenous growth of existing trees

# Mangrove (biomass) loss
betaL, betaD, betaS = np.zeros(n), np.zeros(n), np.zeros(n)

for k in range(n):
    betaL[k], betaD[k], betaS[k] = pick3()

# Peat gain
betaA = np.zeros(n)
betaR = np.zeros(n)
betaV = np.zeros(n)
for k in range(n):
    betaA[k], betaR[k], betaV[k] = pick3()

#Peat loss
betaE = random.uniform(0,1,n)
betaSB = 1 - betaE

# ----------------------
# Elasticity parameters
# ----------------------

# Default range for unknown but likely linear values
r0 = 0.5
r1 = 2


hydP = random.uniform(-1.0, 0, n)
# Mangroves
propM = random.uniform(r0, r1, n) 
propS = random.uniform(-3, 0.0, n)
growM = random.uniform(0, 1, n)
growS = random.uniform(-2,0.0,n)

evaptM = random.uniform(0,1,n)

precipBeta = random.uniform(0,0.3,n) #https://journals.ametsoc.org/doi/abs/10.1175/1520-0442%281993%29006%3C1077%3AEOCPR%3E2.0.CO%3B2

propPrecip = random.uniform(r0,3,n)
growPrecip = random.uniform(r0,r1,n)

drownHyd = random.uniform(0.0, 5.0, n)
drownM = random.uniform(r0, r1, n) 

stressM = random.uniform(r0, r1, n)
stressS = random.uniform(0.0, 5.0, n)

littM = random.uniform(r0, r1, n)

# Peat soils
accSed = random.uniform(r0, r1, n) 
sedHyd = random.uniform(0.5, 1, n)
accM = random.uniform(r0, r1, n)

retLitt = random.uniform(r0, r1, n)
retHyd = random.uniform(-2.0, 0.0, n)
    
volGrow = random.uniform(r0, r1, n)
volP = random.uniform(-2, -0.5, n)
volHyd = random.uniform(r0,r1,n)
volPrecip = random.uniform(0.5,1,n)

eroM = random.uniform(-3,0, n)

subsMort = random.uniform(r0, 4, n)
subsHyd = random.uniform(r0, r1, n)
subsP = random.uniform(0.5,1, n)

# Salinity
concS = random.uniform(r0, r1, n) #nonlinear values inferred from Teh 2008
concEvapt = random.uniform(r0,r1,n)

concHyd = random.uniform(0.5,3,n)
    
decrS = random.uniform(r0,r1,n)
decrPrecip = random.uniform(r0,r1,n)

evaptS = random.uniform(-2,0,n)


def stability(eigs):
    # take in vector of eigenvalues
    # tell if system stable or not (stable if all eigenvalues are less than zero)
    reals = np.real(eigs)
    if max(reals) < 0:
        result = 1
    else:
        result = 0
    return result

def ruthHurwitz3(coffs):
    a0, a1, a2 = coffs
    cond1 = a0 > 0
    cond2 = a2 > 0
    cond3 = a1*a2 > a0
    return (cond1 and cond2 and cond3)

def schurCohn(lam1,lam2,lam3):
    # https://www.hindawi.com/journals/ddns/2017/6186354/
    # test if Jacobian is locally asymptotically stable
    trace = lam1+lam2+lam3
    det = lam1*lam2*lam3
    minors = lam1*lam2+lam1*lam3+lam2*lam3
    cond1 = np.abs((trace+det))< (1 + minors)
    cond2 = np.abs(minors-trace*det)<(1-det**2)
    return int(cond1 and cond2)

# Construct dataframe to track parameters and associated eigenvalues
# Parameters that are varying
data = {'alphaM':alphaM,'alphaP':alphaP,'alphaS':alphaS,
        'betaG':betaG,'betaP':betaP,'betaD':betaD,'betaL':betaL,'betaS':betaS,
        'betaA':betaA,'betaR':betaR,'betaV':betaV,'betaE':betaE,'betaSB':betaSB,
        'hydP':hydP,'propM':propM,'propS':propS,\
        'growM':growM,'growS':growS,\
        'propPrecip':propPrecip,'growPrecip':growPrecip,\
        'drownHyd':drownHyd,'drownM':drownM,'stressM':stressM,\
        'stressS':stressS,'littM':littM,'accSed':accSed,\
        'sedHyd':sedHyd,'accM':accM,'retLitt':retLitt,'retHyd':retHyd,\
        'volGrow':volGrow,'volP':volP,'volHyd':volHyd,'volPrecip':volPrecip,'eroM':eroM,\
        'subsMort':subsMort,'subsHyd':subsHyd,'subsP':subsP,'concS':concS,\
        'concEvapt':concEvapt,'concHyd':concHyd,'evaptM':evaptM,'evaptS':evaptS,\
        'decrS':decrS,'decrPrecip':decrPrecip,'precipBeta':precipBeta}

eigs = [] #eigenvalue triplets
eigMax = [] # max Eigenvalue list
eigsV = [] #eigenvector triplets
stab = [] #stability (0,1)
determ = [] # determinant of Jacobian
traces = [] # trace of Jacobian
fixType = []
for j in tqdm(range(n)):
    
    dataJ = {k:v[j] for (k,v) in data.items()}
    jac = computeJac(dataJ)
   
    w, v = LA.eig(jac)
    det = LA.det(jac)
    tr = np.trace(jac)
    j2 = np.square(jac)
    lamCoff = -0.5*(tr**2-np.trace(j2))
    
    eigs.append(w)
    eigMax.append(np.real(np.max(w)))
    eigsV.append(v)
    stab.append(stability(w))
    
    determ.append(det)
    traces.append(tr)

    #coffs.append(lamCoff)
    
    #ruth = ruthHurwitz3((det,lamCoff,tr))
    #ruthStab.append(int(ruth))

#p1 = plt.figure(1)
#plt.scatter(traces,determ)
#plt.xlabel('Trace')
#plt.ylabel('Determinant')
#plt.show()

# Characterise types of systems from parameter range
typeList = list(map(typeNode,eigs))
typeCount = {}
for typ in set(typeList):
    typeCount[typ] = typeList.count(typ)



# Compute correlations
corrs = {k:(np.corrcoef(data[k],stab)[0,1]) for (k,v) in data.items() }

#Sort out only the big correlations
numPlot = 10 # number of variables to plot
remNum = len(corrs.items()) - numPlot #number of variables to remove
absCorrs = {k: np.abs(v) for (k,v) in corrs.items() if not isnan(v)}
corrsSorted = sorted(absCorrs.items(), key=lambda x: x[1], reverse=True)
delN = len(corrsSorted) - numPlot
del corrsSorted[-delN:]

bigCorrs = {k: corrs[k] for (k,v) in corrsSorted}


plt.bar(range(len(bigCorrs)), bigCorrs.values())

plt.ylabel('Pearson Correlation Coefficient')
plt.xlabel('Parameter')
plt.xticks(range(len(bigCorrs)), list(bigCorrs.keys()), rotation=60)
plt.show()



'''
posStable = set() # parameters for which higher values => greater stability
negStable = set() # parameters for which higher values => lower stability

for (param,corr) in bigCorrs.items():
    if corr > 0 :
        posStable.add(param)
    else:
        negStable.add(param)

# For max stability
maxStable = {}
minStable = {}

for param in posStable:
    maxStable[param] = max(data[param])
    minStable[param] = min(data[param])
    
for param in negStable:
    maxStable[param] = min(data[param])
    minStable[param] = max(data[param])
    
'''