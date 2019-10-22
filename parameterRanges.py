# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:53:3  2019

@author: Excalibur
Paramater ranges for elasticities
"""
# ----------------------
# Elasticity parameters
# ----------------------

alphaM = (1/3,2)
alphaP = (1, 6)
alphaS = (1/4, 1/2)

alphas = {'alphaM':alphaM, 'alphaP':alphaP, 'alphaS':alphaS}

betaA = (0,1)
betaR = (0,1)
betaV = (0,1)

betaE = (0,1)
betaSB = (0,1)

betaP = (0,1)
betaG = (0,1)

betaS = (0,1)
betaD = (0,1)
betaL = (0,1)
betas = {'betaG':betaG, 'betaP':betaP, 'betaD':betaD, 'betaS':betaS, 'betaL':betaL,
         'betaA':betaA, 'betaR':betaR, 'betaV':betaV, 'betaE':betaE, 'betaSB':betaSB}

# Default range for unknown but likely linear values
r0 = 0.5
r1 = 2


hydP = (-2.0,0)
# Mangroves
propM = (r0,r1)
propS = (-2.0,0)
growM = (r0, r1)
growS = (-1,0.0)

precipBeta = (0,1)

propPrecip = (r0,3)
growPrecip = (r0,r1)

drownHyd = (0.0, 5.0)
drownM = (r0, r1) 

stressM = (r0, r1)
stressS = (0.0, 5.0)

littM = (1, 2)

# Peat soils
accSed = (0, 1) 
sedHyd = (r0, r1)
accM = (0, r1)

retLitt = (r0, r1)
retHyd = (-2.0, 0.0)
    
volGrow = (r0, r1)
volP = (r0, r1)
volHyd = (r0,r1)
volPrecip = (0.5,1)

eroM = (-3,0)

subsMort = (r0, 4)
subsHyd = (r0, r1)
subsP = (r0,r1)

# Salinity
concS = (r0, r1) #nonlinear values inferred from Teh 2008
concEvapt = (r0,r1)

concHyd = (0.5,3)
    
decrS = (r0,r1)
decrPrecip = (r0,r1)

evaptS = (-3,-0.5)
evaptM = (r0,r1)


mangs = {'propM':propM, 'propS':propS, 'growM':growM,'growS':growS, 'drownHyd':drownHyd, \
         'drownM':drownM,'stressM':stressM, 'stressS':stressS, 'littM':littM,\
         'propPrecip':propPrecip,'growPrecip':growPrecip,\
         'evaptM':evaptM,'precipBeta':precipBeta}

peats = {'accSed':accSed, 'sedHyd':sedHyd, 'accM':accM,\
         'retLitt':retLitt, 'retHyd':retHyd, 'volGrow':volGrow,
         'volP':volP,'volPrecip':volPrecip, 'eroM':eroM, 'subsMort':subsMort,\
         'subsHyd':subsHyd, 'subsP':subsP, 'hydP':hydP,'volHyd':volHyd}

salts = {'concEvapt':concEvapt,'concHyd':concHyd, 'concS':concS, 'decrS':decrS,
         'decrPrecip':decrPrecip,'evaptS':evaptS}

ranges = {**alphas, **betas, **mangs, **peats, **salts}
