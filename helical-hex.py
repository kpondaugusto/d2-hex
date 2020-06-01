#!/usr/bin/python3

from math import *


from mpl_toolkits import mplot3d
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
mu=3.5e-5 # Pa*s

L=10*0.0254 #m length of tube

rho=163.0 # kg/m^3

Cp=6565.4 # J/(kg*K)

N=5 #number of turns

Ngrooves=1 # number of grooves

D=4.76*0.0254 # 0.015949 #m diameter of tube, 0.015949 from optimizing dp in backwards-hex-turbulent-tube.py


wprime= 0.01 #m width of groove

uprime= 0.01 # m width between grooves

depth=0.01 # m depth of groove

sinalpha=(Ngrooves*(wprime + uprime))/(pi*D) #pitch angle

alpha=arcsin(sinalpha)

print('The pitch angle is %f' %sinalpha)

print('Alpha is %f' %alpha)

Lprime=L/sinalpha #m length of wound groove

print('The length of the groove is %f m.' %Lprime)

turns=Lprime/(pi*D)



print('Coiling around a Cu rod of diameter %f m would require %f turns'%(D,turns))

#based off of sketch w/ jeff

w=wprime*tan(alpha) # m

ahelix=Ngrooves*wprime*depth #Arect+2*Atri ?? m^2 area of one helical groove/fin thing

print('The area of the helical fins is %f m^2.'%ahelix)

phelix=Ngrooves*(2*depth+2*wprime) #m

print('The perimeter of the helical grooves is %f m.' %phelix)


Dh=4*ahelix/phelix #m

print('Hydraulic diameter %f m'%Dh)
print()

#Aw=ptot*L
#print('Area of cold wall %f m^2'%Aw)

G=mdot/ahelix # (kg/(m^2*s)) mass flow rate per unit area

print('Mass flux (G) is %f kg/(m^2*s)'%G)


Re=Dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)

fRe=24.00*4

#Creating an elif for f based on Barron eq'ns

if Re < 2300 :
    f=fRe/Re
    #f=64/Re #assuming cicular tube
    print('The laminar friction factor is %f.' %f)
elif 3500 > Re > 2300 :
    f=1.2036*Re**(-0.416) #from vijayan
    print('The friction factor is in between laminar and turbulent')
elif Re > 3500 :
    f=0.316*Re**(-0.25)
    print('The turbulent friction factor is %f.' %f)


B1=1.174*((3.7e-5)/(3.68e-5))**(0.14) #viscosity taken from cams sheets
print('This is B1 %f.' %B1)


if Re < 3500 :
    print('It is laminar or in between')
elif Re > 3500 :
    jh=0.023*Re**(-0.2)*B1
    print('The Colburn factor for the turbulent flow is %f.' %jh)



print()

Cp=CP.PropsSI('C','P',p,'T',T,fluid) # (kg/(m*K)) found from coolprop -
                                    # found via a table

Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

print('The Prandtl Number is %f.'%Pr)

#If turb

if Re < 3500 :
    Nu=4.8608
    print('Nu=4.8608 because the flow is laminar')
elif Re > 3500 :
    Nuturb=jh*Re*Pr**(1./3.)
    print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
print()

if Re < 3500 :
    hc=Nu*kt/Dh # Barron eq'n 6.15
    print('The heat transfer coefficient for laminar flow is %f W/(m^2*K)'%hc)
elif Re > 3500 :
    hc=Nuturb*kt/Dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
    print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

#Ntu=hc*Aw/(mdot*Cp)
#print('The number of transfer units is %f'%Ntu)
#print()

#T1=Tin
#T2=T1-(T1-Tw)*(1-exp(-Ntu))
#T2=Tw+(T1-Tw)*exp(-Ntu)

#Qtotal=mdot*Cp*(T1-T2) # Eq. (6.43) of Barron

#print('For inlet temperature %f K and wall temperature %f K'%(T1,Tw))
#print('the outlet temperature is %f K'%T2)
#print('and the total heat transfer rate is %f W'%Qtotal)
#print()


dp=(f*L*G**2)/(Dh*2*rho) # (Pa) pressure drop
# unit check:
# [L]=m
# [G**2]=kg^2/(s^2*m^4)
# [dh]=m
# [rho]=kg/m^3
# So [p]=(kg^2/(s^2*m^3))/(kg/m^2)=kg/(s^2*m)=(kg*m/s^2)/m^2=[force]/[area]
# =Pa (as expected)

print('The pressure drop is %f Pa'%dp)


#finding optimal Ngrooves for optimal hc

n = np.arange(1,5,1)


hc = []

for i in range(len(n)):

    value=(Pr**(1/3)*B1*0.023*mdot**(0.8)*kt*2*(wprime+depth))/((2*mu*(wprime+depth)*n[i])**(0.8)*wprime*depth)
    
    hc.append(value)

plt.plot(n, hc)

plt.show()



dp = []

for i in range(len(n)):

    value=(0.316*mdot**(7/8)*L*(wprime + depth)**(5/4)*mu**(-3/4)*2**(1/4)*n[i]**(-3/4))/((wprime*depth)**2*rho)
    
    dp.append(value)

plt.plot(n, dp)

plt.show()


