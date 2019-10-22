# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 20:44:36 2019

@author: Excalibur
"""
import numpy as np

def typeNode(lams):
    lam1, lam2, lam3 = lams
    
    if any(lam == 0 for lam in lams):
        return "Non hyperbolic!"
    elif all(np.isreal(lam) for lam in lams):
        signs = list(map(np.sign,list(lams)))
        if all(lam<0 for lam in lams):
            return 'Attracting Node'
        elif all(lam > 0 for lam in lams):
            return 'Repelling Node'
        elif sum(signs) == -1:
            return 'Saddle'
        elif sum(signs) == 1:
            return 'Saddle'
    elif any(np.isreal(lam) for lam in lams):
            realLams = np.real(list(lams))
            signs = list(map(np.sign,list(realLams)))
            if sum(signs) == -3:
                return "Stable Focus-Node"
            elif sum(signs) == 3:
                return "Unstable Focus-Node"
            else:
                return "Saddle-Focus Point"
