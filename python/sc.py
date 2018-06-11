# +
#
# PROGRAM: function Sc
#
#
# PURPOSE: finds the critical supersaturation from a dry diameter and
#          a set of kappa values, volume fractions, and solubilities. 
#
#
# AUTHOR: Markus Petters (petters@atmos.colostate.edu)
#         Department of Atmospheric Science
#         Colorado State University, Fort Collins, CO
#
#
# COMMENTS: This code is published in Petters and Kreidenweis,
#           ACPD, 2008. Default values for T are 298.15K and sigma =
#           0.072 J m-2 
#
#
# EXAMPLE: This code can be executed directly from a Python
#          interpreter (for example IDLE for Windows) or from the
#          command line in linux (e.g. "python sc.py"). This code was
#          tested on Python 2.5 (r25:51908) and Python 2.5.1
#          (r251:54863)
#
#
# MODIFICATION HISTORY:
#   
# Copyright (C) 2008, Markus Petters, Department of Atmospheric
# Science, Colorado State Univsersity
#
# This software may be used, copied, or redistributed as long as it is not
# sold and this copyright notice is reproduced on each copy made.  This
# routine is provided as is without any express or implied warranties
# whatsoever.
#
#-

from math import *
from numpy import *

def Sc(Dd, ki, ei, Ci, T=298.15, sigma=0.072, gmax=10, g=1+1e-11, inc=1.01):    
    A = 8.69251e-6*sigma/T
    xi = map(lambda x: x if x < 1 else 1, Ci*(g**3.0 - 1.0)/ei)
    k = dot(ki, ei*xi)
    S = lambda D, Dd, k, A: (D**3.0-Dd**3.0)/(D**3.0-Dd**3.0*(1.0-k))*exp(A/D)
    f = lambda x,y: x if x > y else y
    return 1 if g > gmax else reduce(f,[S(g*Dd,Dd,k,A), \
    Sc(Dd,ki,ei,Ci,T=T,sigma=sigma,gmax=gmax,g=g*inc,inc=inc)]) 

print Sc(100e-9, array([0.6,0.2]), array([0.5,0.5]), array([inf,0.1]))
print Sc(100e-9, array([0.6]), array([1]), array([inf]), T=273.15,inc=1.01)
