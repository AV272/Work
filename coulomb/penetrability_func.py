#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:51:56 2020

@author: lkst
"""
import math
import cmath
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpmath import *


def C_l(l, eta, z):
    F = coulombf(l, eta, z) 
    G = coulombg(l, eta, z)
    
    P_l = 1/(F**2 + G**2)
    delta_l = -cmath.atan(F/G)
    ch = gamma(l +1 + j*eta)
    phi_l = cmath.phase(ch)
    return(sqrt(P_l*exp(j*(delta_l + phi_l))))


def F2(x):
    return(coulombg(0,0,x))

def C_l_E(l,E):
    z1=z2=1
    el = sqrt(1/137)
    h = 1
    M1 = 1875.61 # MeV
    M2 = 1875.61
    M = M1*M2/(M1+M2)
    #eta = z1*z2*sqrt(M/E)*0.1574
    eta = 0.16*sqrt(1/E)
    z2 = 0.219*sqrt(E)*6.63
    return(cmath.phase(C_l(l,eta,z2)))

    
def C0(E):
    return(C_l_E(0, E))
def C1(E):
    return(C_l_E(1, E))
def C2(E):
    return(C_l_E(2, E))
def C3(E):
    return(C_l_E(3, E))

#plot(C_l_E, [0,1500], [0,1])
#plot(C_l_E, [0,1500])
plot([C0,C1,C2,C3], [0,1.5])