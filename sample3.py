from numpy import random

def pick3():
    d = random.uniform(0,1)
    c = random.uniform(0+d/2.,1-d/2.)
    b1 = c-d/2.
    b2 = d
    b3 = 1-b1-b2
    picks = [b1,b2,b3]
    random.shuffle(picks)
    return picks