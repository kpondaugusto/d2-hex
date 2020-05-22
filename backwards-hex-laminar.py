#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
from scipy import interpolate
from scipy.integrate import odeint

mdot=0.004 # kg/s
mu=3.5e-5 # Pa*s
Re=329.
p=4*mdot/(Re*mu) # m
print('Flow perimeter %f m'%p)

rho=163.0 # kg/m^3
f=64./Re # laminar

dp=8.5 # Pa

L=10*0.0254 # m

#dp=(f*mdot**2/(8*rho))*L*p/A**3
#L*p/A**3=dp/(f*mdot**2/(8*rho))
alpha=dp/(f*mdot**2/(8*rho))

Cp=6565.4 # J/(kg*K)
Nu=4.8608 # annular or planar Nusselt number laminar
kt=0.104 # W/(m*K) a check on this number from

Ntu=3.24 # achieved in other simulation

#hc=Nu*kt*p/(4*A)
#Ntu=hc*Aw/(mdot*Cp)

#(p/A)*Aw=4*Ntu/(Nu*kt/(mdot*Cp))
beta=4*Ntu/(Nu*kt/(mdot*Cp))

#Aw/L*A**2=Pin/A**2=beta/alpha

# solve beta and alpha for A and Pin

Pin=(beta**3/(alpha*p**2*L**2))**(1./3.)
print('Perimeter needed for sufficient heat transfer %f m'%Pin)
a=p*Pin*L/beta
print('Area needed for flow %f m^2'%a)

# p = pi*Din + pi*Dout + 2*n*d
# Pin = pi*Din + 2*n*d
# if Din =~ Dout = HEXD
# then (p-Pin)/pi = HEXD

HEXD=(p-Pin)/pi
print('HEX diameter %f m'%HEXD)
if(HEXD<0):
    print('Error!!!')
    print('Perimeter needed for heat transfer exceeds required flow perimeter!!!')

twond=Pin-pi*HEXD
print('Total fin depth %f m'%twond)

n=124
fin_depth=twond/2/n
print('For %d fins, the depth of each is %f m'%(n,fin_depth))
