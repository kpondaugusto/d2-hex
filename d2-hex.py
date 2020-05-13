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

#friction factor

f_lam=64/Re # for Re < 2300 (laminar flow), dimentionless

f_turb=0.316*Re**(-0.25) #for 3500 < Re < 20000 (turbulent flow), again dimentionless

print('The friction factor is %f.'%f_lam)

#Colburn J factor

jH=0.023*Re**(-0.2) # dimentionless

print('The J factor is %f.'%jH)

c=3e8 #(m/s) speed of light

#thermal conductivity

kt=0.104 #W/mK

#Prandtl Number

Pr=(mu*c)/(kt) #yes still dimentionless

print('The Prandtl Number is %f.'%Pr)


#Nusselt Number

Nu=jH*Re*Pr**(1/3) #more dimentionless numbers


print('The Nusselt Numebr is %f.'%Nu)

ro=165#(kg/m^3) denstiy

L=4.007 #(m) length of tube(s)

#pressure drop

p=(f_lam*L*G**2)/(dh*2*ro) #(Pa)


print('The pressure drop is %f Pa.'%p)

#convective heat transfer coefficient

hc=(dh)/(Nu*kt) #m^2k/W

print('The convective heat transfer coefficient is %f m^2k/W.'%hc)



