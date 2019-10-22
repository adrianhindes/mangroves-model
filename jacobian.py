# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:17:03 2019

@author: hindesa
Return Jacobian for mangroves model for given parameter set
"""
import numpy as np

def computeJac(data):
    alphaM = data['alphaM']
    alphaP = data['alphaP']
    alphaS = data['alphaS']
    
    betaG = data['betaG']
    betaP = data['betaP']
    
    betaD = data['betaD']
    betaS = data['betaS']
    betaL = data['betaL']
    
    betaA = data['betaA']
    betaR = data['betaR']
    betaV = data['betaV']
    
    betaE = data['betaE']
    betaSB = data['betaSB']
    
    hydP = data['hydP']
    
    propM = data['propM']
    propS = data['propS']
    growS = data['growS']
    growM = data['growM']
    
    propPrecip = data['propPrecip']
    growPrecip = data['growPrecip']
    
    drownHyd = data['drownHyd']
    drownM = data['drownM']
    stressM = data['stressM']
    stressS = data['stressS']
    
    littM = data['littM']
    accSed = data['accSed']
    sedHyd = data['sedHyd']
    accM = data['accM']
    retLitt = data['retLitt']
    retHyd = data['retHyd']
    
    volGrow = data['volGrow']
    volP = data['volP']
    volHyd = data['volHyd']
    volPrecip = data['volPrecip']
    
    eroM = data['eroM']
    subsMort = data['subsMort']
    subsHyd = data['subsHyd']
    subsP = data['subsP']
    
    
    concEvapt = data['concEvapt']
    evaptM = data['evaptM']
    concHyd = data['concHyd']
    
    decrPrecip = data['decrPrecip']
    
    precipBeta = 0 #data['precipBeta']
    evaptM = data['evaptM']
    evaptS = data['evaptS']
    

    

    
    # Define Jacobian matrix elements
    # Note syntax here is not strictly correct - dmdm == (dm/dt)/dm
    
    betaSB = 1-betaE
    betaV = 1-betaA-betaR
    
    mortD = betaD/(betaD+betaS)
    mortS = betaS/(betaD+betaS)

    dPropdM = propM+propPrecip*precipBeta*evaptM
    dGrowdM = growM+growPrecip*precipBeta*evaptM

    dmdm = betaP*dPropdM +betaG*dGrowdM-betaS*stressM -betaD*drownM -betaL*littM
         
    dmdp = -1*betaD*hydP*drownHyd

    dPropdS = propS+propPrecip*precipBeta*evaptS
    dGrowdS = growS+growPrecip*precipBeta*evaptS

    dmds = betaP*dPropdS + betaG*dGrowdS-betaS*stressS

    dVoldM = volGrow*(growM+growPrecip*precipBeta*evaptM)+volPrecip*precipBeta*evaptM
    dSubsdM = subsMort*(mortD*drownM+mortS*stressM)
    
    dpdm = betaA*accM +betaR*retLitt*littM + betaV*dVoldM - betaE*eroM -betaSB*dSubsdM

    dVoldP = volHyd*hydP+volP
    dSubsdP = subsHyd*hydP + subsP
            
    dpdp = hydP*(betaA*accSed*sedHyd + betaR*retHyd)+betaV*dVoldP-betaSB*dSubsdP


    dVoldS = volGrow*(growS+growPrecip*precipBeta*evaptS)+volPrecip*precipBeta*evaptS
    dSubsdS = subsMort*(mortS*stressS)

    dpds = betaV*dVoldS - betaSB*dSubsdS
    
    dsdm = evaptM*(concEvapt - decrPrecip*precipBeta)
    dsdp = concHyd*hydP
    dsds = concEvapt*evaptS - decrPrecip*precipBeta*evaptS
    
    # alpha paramater array
    alphas = np.array([ [alphaM, 0, 0], [0, alphaP, 0], [0, 0, alphaS]])
    alphas = alphas.astype(float)

    R1 = [dmdm, dmdp, dmds]
    R2 = [dpdm, dpdp, dpds]
    R3 = [dsdm, dsdp, dsds]

    jac0 =  np.array([R1,R2,R3])
    jac0 = jac0.astype(float)
    jac = np.matmul(alphas, jac0)
    
    return jac