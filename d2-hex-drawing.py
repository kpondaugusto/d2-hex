#!/usr/bin/python3

from math import pi,cos,sin,asin
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



plt.show()
