#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

d2=4.76*0.0254 # (m) inner diameter of the outer tubular housing
d1=4.75*0.0254 # (m) diameter of the inner cold cylinder before
               # cutting any grooves

fig,ax=plt.subplots()

ax.set_xlim([-d2/2,d2/2])
ax.set_ylim([-d2/2,d2/2])

# outer housing
circle1=plt.Circle((0,0),d2/2,color='r',fill=False)
ax.add_artist(circle1)

groove_depth=0.1*0.0254 # m
groove_width=0.06*0.0254 # m
ngrooves=124

# shorter names
r=d1/2
w=groove_width
d=groove_depth

for groove in range(ngrooves):
    theta=360./ngrooves
    center_angle=groove*theta

    alpha=center_angle*pi/180
    x=r*cos(alpha)
    y=r*sin(alpha)

    dalpha=asin(w/2/r)

    xsl=r*cos(alpha+dalpha)
    ysl=r*sin(alpha+dalpha)
    xel=xsl-d*cos(alpha)
    yel=ysl-d*sin(alpha)
    
    xsr=r*cos(alpha-dalpha)
    ysr=r*sin(alpha-dalpha)
    xer=xsr-d*cos(alpha)
    yer=ysr-d*sin(alpha)
    
    # draw line at edge of each groove
    line=plt.plot([xsl,xel],[ysl,yel],color='black')
    line=plt.plot([xsr,xer],[ysr,yer],color='black')
    # and bottom of groove
    line=plt.plot([xel,xer],[yel,yer],color='black')

    alpha_deg=alpha*180/pi
    dalpha_deg=dalpha*180/pi
    arc=patches.Arc((0,0),2*r,2*r,0,alpha_deg+dalpha_deg,alpha_deg+theta-dalpha_deg)
    ax.add_patch(arc)

# Calculation of perimeter of all those grooves
pgroove=2*d+w # inner "U" of a groove
parc=(theta-2*dalpha_deg)*pi/180*r # outer arc length between two grooves
perimeter=ngrooves*(pgroove+parc)

print('The length of the inner U of a groove is %f m'%pgroove)
print('The length of an arc between two grooves is %f m'%parc)
print('The length around all those fins and grooves is %f m'%perimeter)

print()

# Calculation of area for flow
annulus=pi*(d2**2-d1**2)/4
agroove=d*w # area of one groove # approximately
aeps=((2*dalpha)/(2*pi))*pi*r**2-2*0.5*(w/2)*(r*cos(dalpha)) # area
                                                             # between
                                                             # the arc
                                                             # and the
                                                             # area of
                                                             # the
                                                             # groove
                                                             # on the
                                                             # previous
                                                             # line
agroove_total=agroove+aeps
agrooves=agroove_total*ngrooves
area=annulus+agrooves


print('Area of annular region between cylinders %f m^2'%annulus)
print('Additional area cut out by each groove %f m^2'%(agroove+aeps))
print('Area of all the grooves %f m^2'%agrooves)
print('Total flow area %f m^2'%area)
print()

a=area # (m^2) area for fluid flow
p=perimeter+pi*d2 # (m) perimeter of flow region

dh=4*a/p # (m) hydraulic diameter

print('Flow area %f m^2'%a)
print('Flow perimeter %f m'%p)
print('Hydraulic diameter %f m'%dh)
print()

L=10.*0.0254 # (m) length of tube(s)
Aw=perimeter*L
print('Length of tubes %f m'%L)
print('Area of cold wall %f m^2'%Aw)
print()

mdot=0.004 # (kg/s) mass flow rate
G=mdot/a # (kg/(m^2*s)) mass flow rate per unit area

print('G is %f kg/(m^2*s)'%G)

# 350 micro-poise at 20.5 K.... this is at saturation, but probably
# close enough to 20 psia.
# https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote641.pdf 

mu=3.5e-5 # Pa*s

Re=dh*G/mu # should be dimensionless
print('The Reynolds number is %f'%Re)



# thermal conductivity

kt=0.104 # W/(m*K) a check on this number from
         # https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote641.pdf
         # gives the value of about 1.05 mW/(cm*K) = 0.105 W/(m*K) at
         # a temperature of about 22 K.  It says in this reference
         # that the data is uncertain at the 25% level, though.
#

# specific heat

import CoolProp.CoolProp as CP
import numpy as np
fluid='Deuterium'

p_psi=20. # PSI
p=p_psi*6894.76 # Pa

Tin=23.4 # (K) inlet temp
Tw=20.7 # (K) temperature of cold wall
#Ts=(Tw+Tin)/2 # (K) temperature of film, with which we will exchange heat.

T=Tin

Cp=CP.PropsSI('C','P',p,'T',T,fluid) # (kg/(m*K)) found from coolprop -
                                    # found via a table

Pr=(mu*Cp)/(kt) # yes still dimensionless
               # because (Pa*s)*(J/(kg*K))/(W/(m*K))
               # =((kg*m/(s^2*m^2))*s)*(W*s/(kg*K))*((m*K)/W) = 1

print('The Prandtl Number is %f.'%Pr)

#Nu_3_ii=4.8608 # Table 90 of Shah and London

Nu=3.657 # circular duct Barron Eq. (6.30)
#Nu=7.541 # Table 6.2 Barron -- parallel plate =~ thin annulus?

fRe=64.00 # circular duct Barron Eq. (6.27)
#fRe=96.00 # Table 6.2 Barron

# calculate heat transfer -- inner surface to fluid
hc=Nu*kt/dh
print('The heat transfer coefficient is %f W/(m^2*K)'%hc)

Ntu=hc*Aw/(mdot*Cp)
print('The number of transfer units is %f'%Ntu)
print()

T1=Tin
T2=T1-(T1-Tw)*(1-exp(-Ntu))
#T2=Tw+(T1-Tw)*exp(-Ntu)

Qtotal=mdot*Cp*(T1-T2) # Eq. (6.43) of Barron

print('For inlet temperature %f K and wall temperature %f'%(T1,Tw))
print('the outlet temperature is %f'%T2)
print('and the total heat transfer rate is %f'%Qtotal)
print()

# calculate pressure drop

f=fRe/Re
print('The friction factor is %f'%f)

rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
print('The density is %f kg/m^3'%rho)

dp=(f*L*G**2)/(dh*2*rho) #(Pa)
print('The pressure drop is %f Pa'%dp)

#plt.show()
