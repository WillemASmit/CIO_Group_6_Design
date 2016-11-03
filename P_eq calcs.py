# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 19:10:53 2016

@author: Wian
"""
import numpy as np
from scipy.optimize import fsolve
#mu = 0.000033
#e = 0.046
#rho = 800
#
#def P_eq(W,D):
#    
#    d = D*1000
#    def colebrook(e,D, Re):
#        def fn1(f):
#            return 1/(-2*np.log10(e/D/3.7+2.51/Re/np.sqrt(f)))**2
#        error = 10E-10
#        x = 0.02
#        x_n = 10
#        while abs(x-x_n)>error:
#            x_n = x
#            x = fn1(x_n)
#        return x
#
#    Re = 354*W/d/mu
#    f_prime = colebrook(e,D,Re)
#    Lp = 23
#    Lr = (48*30 + 135)*D
#    Le = Lp + Lr
#    del_P = (62544*f_prime*Le*W**2) / (rho*d**5) 
#    return del_P
    
#==============================================================================
# L-101, L-102 en L-103
#==============================================================================
print("L-101, L-102 en L-103")
print("======================================================================")
rho = 1027
e = 0.015 #316 Stainless
mu = 0.0041
W = 35000
V = W/rho
#Omdat V = 43.75 kies 'n del_P1m van 0.25-0.9. 
del_P1 = 0.25 #kPa
n = 0.6043
k = 0.05452

D = 50
D_n = 0.2
while abs(D-D_n)>1E-10:
    D = D_n
    d = D*1000
    u = W/rho/3600/(np.pi/4*D**2)
    Re = 8*(n/(6*n+2))**2*rho*u**(2-n)*D**(-n)
    f = 8*((6*n+2)/n)**n * k*u**(n-2)*D**(-n)/rho
    d = (62544*f*1*W**2/(rho*del_P1))**(1/5)
    D_n = d/1000

def Solv(D):
    u = W/rho/3600/(np.pi/4*D**2)
    return -del_P1 + ((6*n+2)/2)**n*4*k*u**n*D**(-1*(n+1))

#d = fsolve(Solv,0.2)
print("The first estimated diameter: ",d," mm")
print("""Choose a 3.5" AISI 316L S10 pipe with ID 3.26"(82.804mm)""")

d = 82.804
D = d/1000
u = W/rho/3600/(np.pi/4*D**2)
Re = 8*(n/(6*n+2))**2*rho*u**(2-n)*D**(-n)
f = 8*((6*n+2)/n)**n * k*u**(n-2)*D**(-n)/rho
Pf1m = (f/D + 0.5)*rho*u**2/2000
print("The checked delPf1m: ", Pf1m," kPa")

Lp = 23 + 5*u*60
Lr = (32*30 + 135 + 9*20)*D #Elbows, CheckValve, T
Le = Lp + Lr
del_Pf = (f*Le/D + 0.5)*rho*u**2/2000
print("The frictional pressure drop is: ",del_Pf, ' kPa')

#Thickness
Peq = 30 + 32 + 80 + 90
Pel = 2.8*rho*9.81/1000

M = 1.125
P = del_Pf + Peq + Pel 
P = P/1000
S = 0.5*(206.843)
E = 1
Y = 0.4
C = 0.58

t_min = M*(P*d/(2*(S*E+P*Y)) + C)
print("The minimum thickness is: ",t_min,' mm')
#==============================================================================
# L-104 en L-105
#==============================================================================
print()
print("L-104 en L-105")
print("======================================================================")
rho = 1000
e = 0.015 #316 Stainless
mu = 8.9E-4
W = 34000
V = W/rho
#Omdat V = 18 kies 'n del_P1m van 0.35-1.35. 
del_P1 = 0.35 #kPa

def Re_fun(D,W):
    d = D*1000
    return 354*W/d/mu

def colebrook(D, Re):
    D = D*1000
    def fn1(f):
        return 1/(-2*np.log10(e/D/3.7+2.51/Re/np.sqrt(f)))**2
    error = 10E-10
    x = 0.02
    x_n = 10
    while abs(x-x_n)>error:
        x_n = x
        x = fn1(x_n)
    return x

D = 97.3836
D_n = 0.2
while abs(D-D_n)>1E-2:
    D = D_n
    d = D*1000
    Re = Re_fun(D,W)
    f = colebrook(D,Re)
    d = (62544*f*1*W**2/(rho*del_P1))**(1/5)
    D_n = d/1000

print("The first estimated diameter: ",d," mm")
print("""Choose a 4" AISI 316L S5 pipe with ID 3.86"(97.3836mm)""")

d = 97.3836
D = d/1000

u = W/(np.pi/4*D**2)/3600/rho
Re = Re_fun(D,W)
f = colebrook(D,Re)
del_Pf1m = (f*1/D + 0.5)*rho*u**2/2000
print("The checked delPF1m: ",del_Pf1m, ' kPa')

Lp = 10
Lr = (24*30)*D
Le = Lp + Lr

del_Pf = (f*Le/D + 0.5)*rho*u**2/2000
print("The frictional pressure drop is: ",del_Pf, ' kPa')

#Thickness
Peq = 100
Pel = 1.2*1000*9.81/1000

#Note P will change to the max head of the pump

M = 1.125
P = del_Pf + Peq + Pel
P = P/1000
S = 0.5*(206.843)
E = 1
Y = 0.4
C = 0.58

t_min = M*(P*d/(2*(S*E+P*Y)) + C)
print("The minimum thickness is: ",t_min,' mm')