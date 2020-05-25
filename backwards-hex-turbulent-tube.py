#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
from scipy import interpolate
from scipy.integrate import odeint

# Reference for this idea
# https://uspas.fnal.gov/materials/10MIT/Lecture_3.2.pdf

mdot=0.004 # kg/s
mu=3.5e-5 # Pa*s
Re=10000.
p=4*mdot/(Re*mu) # m
print('Flow perimeter %f m'%p)
diameter=p/pi
print('Diameter of tube %f m'%diameter)

rho=163.0 # kg/m^3

Cp=6565.4 # J/(kg*K)

# heat xfer correlations
Pr=2.209513 # taken from d2-hex-drawing.py
B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
print('This is B1 %f.' %B1)
jh=0.023*Re**(-0.2)*B1
Nu=jh*Re*Pr**(1./3.)
print('According to correlations, Nu=%f'%Nu)


kt=0.104 # W/(m*K) a check on this number from

Ntu=3.0

hc=Nu*kt/diameter
print('The heat transfer coefficient is %f W/m^2-K'%hc)

#Ntu=hc*pi*D*L/(mdot*Cp)
L=Ntu/(Nu*kt*pi/(mdot*Cp))

print('Tube length needed for sufficient heat transfer %f m'%L)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

dp=(8*f*mdot**2/(rho*pi**2))*L/diameter**5 # tube
print('Pressure drop %f Pa'%dp)

dcoil=4.*0.0254 # m diameter of coiled tube
turns=L/(pi*dcoil)



print('Coiling around a Cu rod of diameter %f m would require %f turns'%(dcoil,turns))
print('This Cu rod would need to be at least %f m long (for zero Cu wall thickness)'%(turns*diameter))

print()

#optimizing diameter for preferred dp

dp2=10 #Pa
mdot=0.004 # kg/s
mu=3.5e-5 # Pa*s
Ntu=3.0
f=0.02 # following along in the example, have to guess an f 

D=((8*f*mdot**2*L)/(rho*pi**2*dp2))**(1/5)

print('this is D %f m'%D)

turns2=L/(pi*D)

print('Coiling around a HEX of diameter %f m would require %f turns'%(D,turns2))
print('This Cu rod would need to be at least %f m long (for zero Cu wall thickness)'%(turns2*D))


Re2=(4*mdot)/(pi*D*mu)

print('This is Re based on optimal diameter %f' %Re2)

Pr2=(mu*Cp)/(kt) # yes still dimensionless

jh2=0.023*Re2**(-0.2)*B1
Nu2=jh2*Re2*Pr2**(1./3.)
print('According to correlations, Nu=%f'%Nu2)

hc2=Nu*kt/D
print('The heat transfer coefficient is %f W/m^2-K'%hc2)



