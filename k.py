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

rho = 163

uc = 0.116

Cp=CP.PropsSI('C','P',p,'T',T,fluid) # (kg/(m*K)) found from coolprop -
                                    # found via a table

#Pressure drops around system

#2 sudden cont



# For sudden cont from hydrualic_Resistance.pdf for turb??

A1= (pi*0.127**2)/4 #m^2 area of pipe bf exp  L = 0.18641 D = 0.127

A3= (pi*0.03808**2)/4 #m^2 area of pipe after exp  L = 0.02093 D = 0.03808

Kcont=(1/2)*(1-A3/A1)**(3/4)

print('This is the loss coefficient for a sudden contraction %f.' %Kcont)





#pg 100 first text

dpcont1 = Kcont*(rho*uc**2)/2

print('The pressure drop due to the sudden contraction is %f Pa.'%dpcont1)


print()

# For sudden cont from hydrualic_Resistance.pdf for turb??

A12= (pi*0.03808**2)/4 #m^2 area of pipe bf exp  L = 0.02093 D = 0.03808

A32= (pi*0.0127**2)/4 #m^2 area of pipe after exp L =  D = 0.0127

Kcont2=(1/2)*(1-A3/A1)**(3/4)

print('This is the loss coefficient for a sudden contraction %f.' %Kcont2)





dpcont2 = Kcont*(rho*uc**2)/2

print('The pressure drop due to the sudden contraction is %f Pa.'%dpcont2)


print()

#1st cold pipe

D= 0.0134 # m
L= 1.937 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


#print('Flow area %f m^2'%A)
#print('Flow perimeter %f m'%P)
print('Hydraulic diameter %f m'%Dh)
#print()

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


#print('Mass flux (G) is %f kg/(m^2*s)'%G)
#print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
#print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
#print('The Colburn factor for the turbulent flow is %f.' %jh)


Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

#print('The Prandtl Number is %f.'%Pr)

Nuturb=jh*Re*Pr**(1./3.)
#print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
#print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
#print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

#Ntu=hc*Aw/(mdot*Cp)
#print('The number of transfer units is %f'%Ntu)
#print()


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
#print('The density is %f kg/m^3'%rho)

dp1=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop


print('The pressure drop for the first pipe is %f Pa'%dp1)


R = ((f*L)/D )*(1/A**2)

print('The R is %f'%R)

print()

#2 ~ 45 deg turns from newest heat exchanger text (liu etc)

#K45 = 0.00241*(B*45*)

K45 = 0.3

uc=G/rho


dp45 = K45*(rho*uc**2)/2

print('The pressure drop due to the 45deg turn is %f Pa.'%dp45)


print()



#2nd cold pipe


D= 0.0134 # m
L= 1.27 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


#print('Flow area %f m^2'%A)
#print('Flow perimeter %f m'%P)
#print('Hydraulic diameter %f m'%Dh)
#print()

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


#print('Mass flux (G) is %f kg/(m^2*s)'%G)
#print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
#print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
#print('The Colburn factor for the turbulent flow is %f.' %jh)


Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

#print('The Prandtl Number is %f.'%Pr)

Nuturb=jh*Re*Pr**(1./3.)
#print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
#print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
#print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

#Ntu=hc*Aw/(mdot*Cp)
#print('The number of transfer units is %f'%Ntu)
#print()


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
#print('The density is %f kg/m^3'%rho)

dp2=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop
# unit check:
# [L]=m
# [G**2]=kg^2/(s^2*m^4)
# [dh]=m
# [rho]=kg/m^3
# So [p]=(kg^2/(s^2*m^3))/(kg/m^2)=kg/(s^2*m)=(kg*m/s^2)/m^2=[force]/[area]
# =Pa (as expected)

print('The pressure drop for the second pipe is %f Pa'%dp2)

R = ((f*L)/D + K45 )*(1/A**2)

print('The R is %f'%R)

print()

#last cold pipe


D= 0.0135 # m
L= 0.8 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


#print('Flow area %f m^2'%A)
#print('Flow perimeter %f m'%P)
#print('Hydraulic diameter %f m'%Dh)
#print()

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


#print('Mass flux (G) is %f kg/(m^2*s)'%G)
#print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
#print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
#print('The Colburn factor for the turbulent flow is %f.' %jh)


Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

#print('The Prandtl Number is %f.'%Pr)

Nuturb=jh*Re*Pr**(1./3.)
#print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
#print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

#Ntu=hc*Aw/(mdot*Cp)
#print('The number of transfer units is %f'%Ntu)
#print()


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
#print('The density is %f kg/m^3'%rho)

dp3=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop
# unit check:
# [L]=m
# [G**2]=kg^2/(s^2*m^4)
# [dh]=m
# [rho]=kg/m^3
# So [p]=(kg^2/(s^2*m^3))/(kg/m^2)=kg/(s^2*m)=(kg*m/s^2)/m^2=[force]/[area]
# =Pa (as expected)

print('The pressure drop for the third pipe is %f Pa'%dp3)

R = ((f*L)/D )*(1/A**2)

print('The R is %f'%R)

print()


#sudden exp



# For sudden expansion from hydrualic_Resistance.pdf for turb??

Aexp1= (pi*0.0127**2)/4 #m^2 area of pipe bf exp

Aexp2= (pi*((159.5+29.5)*0.0254)**2) #m^2 area of pipe after exp

Kexp=(1-Aexp1/Aexp2)**2

print('This is the loss coefficient for a sudden expansion %f.' %Kexp)

uc=G/rho

dpexp = Kexp*(rho*uc**2)/2

print('The pressure drop due to the sudden expansion is %f Pa.'%dpexp)

R = (Kexp)*(1/Aexp2**2)

print('The R is %f'%R)

print()



# 180 deg loop


K1802=2.2

uc=G/rho

dp180 = K1802*(rho*uc**2)/2

print('The pressure drop due to the 180deg loop is %f Pa.'%dp180)

R = (K1802 )*(1/Aexp2**2)

print('The R is %f'%R)
print()

#sudden contract



# For sudden contraction from hydrualic_Resistance.pdf

A1=Aexp2 #m^2 before cont

Acont3= (pi*0.03175**2)/4 #m^2 after

Kcont3=(1/2)*(1-A3/A1)**(3/4) #where A3 is area of pipe after contraction and a fair ways down so the stream is fully developed again see vena contracta effect for more

print('This is the loss coefficient for a sudden contraction %f.' %Kcont)






dpcont3 = Kcont3*(rho*uc**2)/2

print('The pressure drop due to the sudden contraction is %f Pa.'%dpcont3)

R = (Kcont3)*(1/Acont3**2)

print('The R is %f'%R)

print()

#first hot pipe L4


D= 0.03175 # m
L= 1.5 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m

G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


#print('Mass flux (G) is %f kg/(m^2*s)'%G)
#print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
#print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
#print('The Colburn factor for the turbulent flow is %f.' %jh)


Pr=(mu*Cp)/(kt) # yes still dimensionless

Nuturb=jh*Re*Pr**(1./3.)
#print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
#print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
#print('The density is %f kg/m^3'%rho)

dp4=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop


print('The pressure drop for the fourth pipe is %f Pa'%dp4)

R = ((f*L)/D )*(1/A**2)

print('The R is %f'%R)

print()

#90 deg angle

K902=0.9

uh=G/rho

dp90 = K902*(rho*uh**2)/2

print('The pressure drop due to the 90deg turn is %f Pa.'%dp90)

R = (K902 )*(1/A**2)

print('The R is %f'%R)

print()

#second hot pipe L5


D= 0.0134 # m
L= 1.75 # m length i found i needed in backwards-hex-turbulent-tube.py

P=2*pi*D/2 # m

A=pi*D**2/4 # m^2



Dh=4*A/P # m


G=mdot/A # (kg/(m^2*s)) mass flow rate per unit area


#print('Mass flux (G) is %f kg/(m^2*s)'%G)
#print()

Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

f=0.316*Re**(-0.25)
print('The turbulent friction factor is %f.' %f)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
#print('This is B1 %f.' %B1)

jh=0.023*Re**(-0.2)*B1
#print('The Colburn factor for the turbulent flow is %f.' %jh)


Pr=(mu*Cp)/(kt) # yes still dimensionless

Nuturb=jh*Re*Pr**(1./3.)
#print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
#print()

hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to
print()


rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
#print('The density is %f kg/m^3'%rho)

dp5=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop

print('The pressure drop for the fifth pipe is %f Pa'%dp5)

R = ((f*L)/D )*(1/A**2)

print('The R is %f'%R)


print()





#valve ? sudden exp ?



A1= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kvalve=10

print('This is the loss coefficient for a valve %f.' %Kvalve)

uh=G/rho

dpvalve = Kvalve*(rho*uh**2)/2

print('The pressure drop due to the valve is %f Pa.'%dpvalve)

R = (Kvalve )*(1/A2**2)

print('The R is %f'%R)


print()


# For sudden expansion from hydrualic_Resistance.pdf for turb??

A1= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kexp2=(1-A1/A2)**2

print('This is the loss coefficient for a sudden expansion %f.' %Kexp2)


dpexp2 = Kexp2*(rho*uh**2)/2

print('The pressure drop due to the sudden expansion is %f Pa.'%dpexp2)

R = (Kexp2 )*(1/A2**2)

print('The R is %f'%R)


print()


dptot=dp1+dp2+dp3+dp4+dp5+dp90+dp180+dpcont1+dpcont2+dpexp+dpexp2+(2*dp45)+dpcont3+dpvalve

print('The total pressure drop is %f Pa' %dptot)

print()

#Table 4.4 Vijayan nuclear government paper pg 21, turb flow i think  https://inis.iaea.org/collection/NCLCollectionStore/_Public/32/004/32004813.pdf



