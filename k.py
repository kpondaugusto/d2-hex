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

#parameters

kt=0.104 # W/(m*K) a check on this number from

mdot=0.004 # (kg/s) mass flow rate

mu=3.5e-5 # Pa*s

p_psi=20. # PSI
p=p_psi*6894.76 # Pa

Tin=23.4 # (K) inlet temp
Tw=20.7 # (K) temperature of cold wall
T=Tin

Cp=CP.PropsSI('C','P',p,'T',T,fluid) # (kg/(m*K)) found from coolprop -
                                    # found via a table

#Pressure drops around system


#For the 1st constraction out of the HEX


D= 0.5*0.0254 # m
L= 1*0.0254 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


print('Flow area %f m^2'%A)
print('Flow perimeter %f m'%P)
print('Hydraulic diameter %f m'%Dh)
print()

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


print('Mass flux (G) is %f kg/(m^2*s)'%G)
print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
print('The Colburn factor for the turbulent flow is %f.' %jh)


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

print('The pressure drop for the first sudden expansion is %f Pa'%dp)


print()


#For the 2nd contraction out of the HEX


D= 0.5*0.0254 # m
L= 1*0.0254 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


print('Flow area %f m^2'%A)
print('Flow perimeter %f m'%P)
print('Hydraulic diameter %f m'%Dh)
print()

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


print('Mass flux (G) is %f kg/(m^2*s)'%G)
print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
print('The Colburn factor for the turbulent flow is %f.' %jh)


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

print('The pressure drop for the second sudden expansion is %f Pa'%dp)


print()


#45 deg









#45 deg





#sudden exp










#90











#valve ? sudden exp ? 





















# For 4 <= Re <= 512 from ASME

#for 90deg bend

Re=100

K90=(2.2**(2.19)+(88.98/Re)**2.19)**(1/2.19)

#print('This is the loss coefficient for a 90deg bend %f.' %K90)

#for 180 loop

K180=(1.92**(1.13)+(167.48/Re)**1.13)**(1/1.13)

#print('This is the loss coefficient for a 180deg loop %f.' %K180)

#Table 4.4 Vijayan nuclear government paper pg 21, turb flow i think  https://inis.iaea.org/collection/NCLCollectionStore/_Public/32/004/32004813.pdf

#exp and cont same as below

#90deg bend

K902=0.9

K1802=2.2

# For sudden expansion from hydrualic_Resistance.pdf for turb??

A1= 10 #m^2 area of pipe bf exp

A2= 20 #m^2 area of pipe after exp

Kexp=(1-A1/A2)**2

#print('This is the loss coefficient for a sudden expansion %f.' %Kexp)


# For sudden contraction from hydrualic_Resistance.pdf

A1= 10 #m^2 before cont

A3= 3 #m^2 after

Kcont=(1/2)*(1-A3/A1)**(3/4) #where A3 is area of pipe after contraction and a fair ways down so the stream is fully developed again see vena contracta effect for more

#print('This is the loss coefficient for a sudden contraction %f.' %Kcont)


#from pic looks like 3 cont, 4 90, 1 180, and 1 exp

Ktot=2*K902 + K1802 + Kexp + Kcont

#print('The loss coefficient for the whole loop is %f' %Ktot)








