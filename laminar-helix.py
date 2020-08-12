#!/usr/bin/python3

from math import *

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import CoolProp.CoolProp as CP
import numpy as np
from numpy import *
fluid='Deuterium'



p_psi=20. # PSI
p=p_psi*6894.76 # Pa
kt=0.104 # W/(m*K) a

T=Tin=23.4 # (K) inlet temp
Tw=20.7 # (K) temperature of cold wall

mdot=0.004 # kg/s
mu=3.5e-5 # Pa*s viscosity

L=10*0.0254 #m length of tube

rho=163.0 # kg/m^3


Cp=6565.4 # J/(kg*K)

Ngrooves=1 # number of grooves

D=4.76*0.0254 # 0.015949 #m diameter of tube, 0.015949 from optimizing dp in backwards-hex-turbulent-tube.py

R=D/2

wprime= 0.015 #m width of groove

uprime= 0.01 # m width between grooves

depth=0.015 # m depth of groove

sinalpha=(Ngrooves*(wprime + uprime))/(pi*D) #pitch angle

alpha=arcsin(sinalpha)

print('The pitch angle is %f' %sinalpha)

print('Alpha is %f' %alpha)

Lprime=L/sinalpha #m length of wound groove

print('The length of the groove is %f m.' %Lprime)

turns=Lprime/(pi*D)

print('Coiling around a Cu rod of diameter %f m would require %f turns'%(D,turns))

#based off of sketch w/ jeff

w=wprime/sinalpha # m

print(w)

ahelix=Ngrooves*wprime*depth #Arect+2*Atri ?? m^2 area of one helical groove/fin thing

print('The area of the helical fins is %f m^2.'%ahelix)

phelix=Ngrooves*(2*depth+2*wprime) #m

print('The perimeter of the helical grooves is %f m.' %phelix)

Dh=4*ahelix/phelix #m

print('Hydraulic diameter %f m'%Dh)
print()

G=mdot/ahelix # (kg/(m^2*s)) mass flow rate per unit area

print('Mass flux (G) is %f kg/(m^2*s)'%G)


Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)
print()


rofc = 2 + ((2*R) + (R**2)/2)/(pi)

print(rofc)
Recritmin = 2100*(1 + 12/((D/Dh)**(0.5)))

print('The critical Re Number is %f.' %Recritmin)


print('R/a is %f' %(D/Dh))


print('a/R is %f' %(Dh/D))

#Dean number

De=Re*sqrt(Dh/D)

print('The Deans number is %f.' %De)

f = ((1 + (Dh/D)/3)**2*(De/(88.33)))**(0.5) #1 + 0.015*Re**(0.75)*(Dh/D)**(0.4) #0.1125*(De)**(0.5) #0.1033*De**(0.5)*((1 + 1.729/De)**(0.5)-(1.729/De)**(0.5))**(-3)


if Re < 2300 :
    fs=fRe/Re
    #f=64/Re #assuming cicular tube
    print('The laminar friction factor is %f.' %fs)
elif 3500 > Re > 2300 :
    fs=1.2036*Re**(-0.416) #from vijayan
    print('The friction factor is in between laminar and turbulent')
elif Re > 3500 :
    fs=0.316*Re**(-0.25)
    print('The turbulent friction factor is %f.' %fs)
    fc = 4*f*fs
    print('f is %f'%fc)



#print(Re*(D/Dh)**(-2))

muw=3.68e-5

Cpw=CP.PropsSI('C','P',p,'T',Tw,fluid) # (kg/(m*K))
#print(Cpw)
Pr=Prw=(muw*Cpw)/(kt)

Prb=(mu*Cp)/kt


print('This is fc(R/a)^-1/2 %f' %f)

Nu = 0.0266*((Re**(0.85)/(D/Dh)**(0.15)) + (0.225*(D/Dh)**(1.15)))*Pr**(0.4) #0.0575*Re**(0.33)*De**(0.42)*Pr**(0.43)*(Prb/Prw)
#Nu = ( (4.364+(4.636/(1 + 1342/(Pr*De**2))**2))**3 + (1.816*(De/(1 + (1.15/Pr)))**(3/2)) )**(1/3)
print('This is Nuc/Nus %f' %Nu)

B1=1.174*((3.7e-5)/(3.68e-5))**(0.14)

jh=0.023*Re**(-0.2)*B1

Pr=(mu*Cp)/(kt)

Nus = jh*Re*Pr**(1./3.)

print('Nuc is %f' %(Nu*Nus))

print()


u=G/rho #m/s
print(u)


dp = (fc*Lprime*rho*u**2)/(2*Dh)


print('the pressure drop is %f Pa' %dp)


#dp2 = (fc2*Lprime*rho*u**2)/(2*Dh)


#print('the pressure drop is %f Pa' %dp2)

#dp=(fc*Lprime*G**2)/(Dh*2*rho) # (Pa) pressure drop

#print('The pressure drop is %f Pa'%dp)


print()

h = Nus*Nu*kt/Dh # (mu*Cp*Pr**(-2/3)*Re*0.023*Re**(-0.2))/(Dh)


print('The heat transfer coeff is %f W/m^2K.' %h)
print()


Aw=Ngrooves*(wprime+2*depth)*Lprime
print('Area of cold wall %f m^2'%Aw)


Ntu=h*Aw/(mdot*Cp)
print('The number of transfer units is %f'%Ntu)
print()

T1=Tin
T2=T1-(T1-Tw)*(1-exp(-Ntu))
T2=Tw+(T1-Tw)*exp(-Ntu)

Qtotal=mdot*Cp*(T1-T2) # Eq. (6.43) of Barron

print('For inlet temperature %f K and wall temperature %f K'%(T1,Tw))
print('the outlet temperature is %f K'%T2)
print('and the total heat transfer rate is %f W'%Qtotal)
print()





