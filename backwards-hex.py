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
print(p)

rho=163.0 # kg/m^3
f=0.29

dp=8.5 # Pa

L=10*0.0254 # m

#dp=(f*mdot**2/(8*rho))*L*p/A**3
#L*p/A**3=dp/(f*mdot**2/(8*rho))
alpha=dp/(f*mdot**2/(8*rho))

Cp=6565.4 # J/(kg*K)
Nu=4.8608
kt=0.104 # W/(m*K) a check on this number from

Ntu=3.238987

#hc=Nu*kt*p/(4*A)
#Ntu=hc*Aw/(mdot*Cp)

#(p/A)*Aw=4*Ntu/(Nu*kt/(mdot*Cp))
beta=4*Ntu/(Nu*kt/(mdot*Cp))

#Aw/L*A**2=Pin/A**2=beta/alpha

# take many fins limit Pin=~p

p=(beta**3/(alpha*L**2))*(1./5.)
print(p)
a=p**2*L/beta
print(a)
