#!/usr/bin/python3

from math import pi,cos,sin,asin

d2=4.76*0.0254 # (m) inner diameter of the outer tubular housing
d1=4.65*0.0254 # (m) diameter of the inner cold cylinder

a=pi*(d2**2-d1**2)/4 # (m^2) area for fluid flow
p=pi*(d2+d1) # (m) perimeter of flow region

dh=4*a/p # (m) hydraulic diameter

print('Flow area %f cm^2'%(a*(100)**2))
print('Flow perimeter %f cm'%(p*100))
print('Hydraulic diameter %f cm'%(dh*100))


mdot=0.003 # (kg/s) mass flow rate
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

# Make drawing

import matplotlib.pyplot as plt

circle1=plt.Circle((0,0),d2/2,color='r',fill=False)

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.add_artist(circle1)

ax.set_xlim([-d2/2,d2/2])
ax.set_ylim([-d2/2,d2/2])

import matplotlib.patches as patches

#arc=patches.Arc((0,0),d1/2,d1/2,45,0,90)
#ax.add_patch(arc)

dtop=d1 # m
groove_depth=0.01 # m
groove_width=0.001 # m
ngrooves=200

for groove in range(ngrooves):
    theta=360./ngrooves
    center_angle=groove*theta
    print(groove,center_angle)
    r=dtop/2
    w=groove_width
    d=groove_depth
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

print(pgroove,parc,perimeter)

plt.show()
