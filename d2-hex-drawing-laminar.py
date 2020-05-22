#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
from scipy import interpolate
from scipy.integrate import odeint

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
n=ngrooves=124

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

print(p,perimeter)

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

print('Mass flux (G) is %f kg/(m^2*s)'%G)

# 350 micro-poise at 20.5 K.... this is at saturation, but probably
# close enough to 20 psia.
# https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote641.pdf 

mu=3.5e-5 # Pa*s

Re=dh*G/mu # should be dimensionless
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




#all from Barron

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


#If turb

if Re < 3500 :
    print('Nu=4.8608 because the flow laminar')
elif Re > 3500 :
    Nuturb=jh*Re*Pr**(1./3.)
    print('This is the turbulent Nusselt Number %f.' %Nuturb)
    
print()

#Nu_3_ii=4.8608 # Table 90 of Shah and London

#Nu=3.657 # circular duct Barron Eq. (6.30)
#Nu=7.541 # Table 6.2 Barron -- parallel plate =~ thin annulus?
#NuT=7.541 # Table 86 of Shah and London, thin annulus
#NuT=7.541 # Table 138 of Shah and London, parallel plate
#NuT=4.861 # Table 138 of Shah and London, parallel plate one side insulated
Nu=4.8608 # Eq. (283) Shah and London, parallel plate one side insulated

#fRe=64.00 # circular duct Barron Eq. (6.27)
#fRe=96.00 # Table 6.2 Barron
fRe=24.00*4 # Table Table 86 of Shah and London, thin annulus

# calculate heat transfer -- inner surface to fluid

if Re < 3500 :
    hc=Nu*kt/dh # Barron eq'n 6.15
    print('The heat transfer coefficient for laminar flow is %f W/(m^2*K)'%hc)
elif Re > 3500 :
    hc=Nuturb*kt/dh # Barron eq'n 6.17 makes it incredibly tiny compared to eq'n 6.15 maybe should be using eq'n 6.40 ??
    print('The heat transfer coefficient for turbulent flow is %f W/(m^2*K)'%hc)

Ntu=hc*Aw/(mdot*Cp)
print('The number of transfer units is %f'%Ntu)
print()

T1=Tin
T2=T1-(T1-Tw)*(1-exp(-Ntu))
#T2=Tw+(T1-Tw)*exp(-Ntu)

Qtotal=mdot*Cp*(T1-T2) # Eq. (6.43) of Barron

print('For inlet temperature %f K and wall temperature %f K'%(T1,Tw))
print('the outlet temperature is %f K'%T2)
print('and the total heat transfer rate is %f W'%Qtotal)
print()

# calculate pressure drop

#f=fRe/Re
#print('The friction factor is %f'%f)

rho=CP.PropsSI('D','P',p,'T',T,fluid) # (kg/m^3)
print('The density is %f kg/m^3'%rho)

dp=(f*L*G**2)/(dh*2*rho) # (Pa) pressure drop
# unit check:
# [L]=m
# [G**2]=kg^2/(s^2*m^4)
# [dh]=m
# [rho]=kg/m^3
# So [p]=(kg^2/(s^2*m^3))/(kg/m^2)=kg/(s^2*m)=(kg*m/s^2)/m^2=[force]/[area]
# =Pa (as expected)

print('The pressure drop is %f Pa'%dp)

plt.title('The Cross Section of the Heat Exchanger.')
plt.show()




#Rectangular Channels in Parallel, Barron 6.3.3, and ex. 6.3

#for one

Af=w*d # (m^2)

Pf=2*w + 2*d # (m)

Dh=4*Af/Pf # (m)

print('The hydraulic diameter for one rectangular channel is %f m.' %Dh)

#all

Afall=d*(n*w) # (m^2)

#Pfall=2*d + 2*n*w # (m)
Pfall=2*n*d + 2*n*w # (m) # Jeff edited

Dhall=4*Afall/Pfall # (m)


print('The flow area for all the rectangular channels is %f m^2.' %Afall)
print()

print('The hydraulic diameter for all the rectangular channels is %f m.' %Dhall)

#Mass flux

Grect=mdot/Afall

Rerect=Dhall*Grect/mu # should be dimensionless

print('The Reynolds number for the rectangles is %f'%Rerect) #still laminar!                                                               #albeit larger than                                                           #the cylinder

#b/a

ba=w/d #dimentionless

print('The b/a is %f.' %ba) # b/c 0.600 I am using C1, C2, from Table 6.2

fRerect=C1=59.94

frect=fRerect/Rerect

print('The friction factor is %f.'%frect)

Nurect=C2=3.205 #Where I used Nu(T) for constant Temperature on a wall
                #(taking the one wall to be the outter constant, the
                #other C2=3.896=Nu(Q) for constant heat flux at all
                #four walls

                # Jeff says:  I don't fully get the comment above.
                
hcrect=(Nurect*kt)/Dhall # W/(m^2*K)

print('The hc of the rectangles is %f W/(m^2*K)' %hcrect)

#Wetted area of duct

#Awrect = 2*(w + d)*L # (m^2)
Aw_wet=(2*n*w+2*n*d)*L # (m^2) # Jeff edited # entire wetted wall area -- for friction losses
Aw_heat=(n*w+2*n*d)*L # (m^2) # Jeff edited # area for heat transfer

print('Wetted area of duct is %f' %Aw_wet)
print('Area for heat transfer is %f' %Aw_heat)

#transfer units

Nturect=(hcrect*Aw_heat)/(mdot*Cp)

print('The transfer units is %f' %Nturect)

#heat transfer

Qrect=mdot*Cp*(T1 - Tw)*(1 - exp(-Nturect)) #W

print('The heat transfer is %f W' %Qrect)

#pressure drop

dprect=(frect*L*Grect**2)/(Dhall*2*rho) #Pa

print('The pressure drop is %f Pa.' %dprect)

#heat conduction

dTdx = -(Qtotal/(kt*a))


print('The slope of the temperature curve is %f K/m.' %dTdx)

def ODE(T,x):

   dTdx = float(-Qtotal/(kt*a) )
   
   return dTdx


#inital conditions

T0 = 19 # K

x = np.arange(0, 10, 0.0001)
solnT = odeint(ODE, T0, x) #from scipy
z = interpolate.interp1d(x, solnT[:,0])

#making callable function

def T(x):
    
    T = z(x)
    
    return T

plt.plot(T(x),x)
plt.title('Temperature as a function of position')
plt.show()









