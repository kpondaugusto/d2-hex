#!/usr/bin/python3

from math import pi

d2=4.76*0.0254 # (m) inner diameter of the outer tubular housing
d1=4.65*0.0254 # (m) diameter of the inner cold cylinder

a=pi*(d2**2-d1**2)/4 # (m^2) area for fluid flow
p=pi*(d2+d1) # (m) perimeter of flow region

dh=4*a/p # (m) hydraulic diameter

print('Flow area %f cm^2'%(a*(100)**2))
print('Flow perimeter %f cm'%(p*100))
print('Hydraulic diameter %f cm'%(dh*100))


mdot=0.004 # (kg/s) mass flow rate
G=mdot/a # (kg/(m^2*s)) mass flow rate per unit area

print('G is %f kg/(m^2*s)'%G)

# 350 micro-poise at 20.5 K.... this is at saturation, but probably
# close enough to 20 psia.
# https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote641.pdf 

mu=3.5e-5 # Pa*s

Re=dh*G/mu # should be dimensionless

print('The Reynolds number is %f'%Re)

Re2=4*mdot/(p*mu)

print('The Reynolds number is %f'%Re2)

import CoolProp.CoolProp as CP
import numpy as np
fluid='Deuterium'

T=20 # K
P=101325 # Pa
h=CP.PropsSI('H','T',T,'P',P,fluid)
s=CP.PropsSI('S','T',T,'P',P,fluid)
print(T,P,h,s)
