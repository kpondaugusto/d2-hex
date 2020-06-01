#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
from scipy import interpolate
from scipy.integrate import odeint


import CoolProp.CoolProp as CP
import numpy as np
fluid='Deuterium'

p_psi=20. # PSI
p=p_psi*6894.76 # Pa

Tin=23.4 # (K) inlet temp
Tw=20.7 # (K) temperature of cold wall

#Taking D from optimizing dp in backwards-hex-turbulent-tube.py


D=0.015949 # m
L=4.322831 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

a=pi*D**2/4 # m^2

acoil=

kt=0.104 # W/(m*K) a check on this number from

Dh=4*A/P # m


print('Flow area %f m^2'%A)
print('Flow perimeter %f m'%P)
print('Hydraulic diameter %f m'%Dh)
print()

mdot=0.004 # (kg/s) mass flow rate
G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


print('Mass flux (G) is %f kg/(m^2*s)'%G)
print()

mu=3.5e-5 # Pa*s

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
print('The Colburn factor for the turbulent flow is %f.' %jh)

T=Tin

Cp=CP.PropsSI('C','P',p,'T',T,fluid) # (kg/(m*K)) found from coolprop -
                                    # found via a table

Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

print('The Prandtl Number is %f.'%Pr)

Nuturb=jh*Re*Pr**(1./3.)
print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

#Ntu=hc*Aw/(mdot*Cp)
#print('The number of transfer units is %f'%Ntu)
print()


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
print('The density is %f kg/m^3'%rho)

dp=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop
# unit check:
# [L]=m
# [G**2]=kg^2/(s^2*m^4)
# [dh]=m
# [rho]=kg/m^3
# So [p]=(kg^2/(s^2*m^3))/(kg/m^2)=kg/(s^2*m)=(kg*m/s^2)/m^2=[force]/[area]
# =Pa (as expected)

print('The pressure drop is %f Pa'%dp)




