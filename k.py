#!/usr/bin/python3
from math import *

# For 4 <= Re <= 512 from ASME

#for 90deg bend

Re=100

K90=(2.2**(2.19)+(88.98/Re)**2.19)**(1/2.19)

print('This is the loss coefficient for a 90deg bend %f.' %K90)

#for 180 loop

K180=(1.92**(1.13)+(167.48/Re)**1.13)**(1/1.13)

print('This is the loss coefficient for a 180deg loop %f.' %K180)



# For sudden expansion from hydrualic_Resistance.pdf for turb??

A1= #m^2 area of pipe bf exp

A2= #m^2 area of pipe after exp

Kexp=(1-A1/A2)**2

print('This is the loss coefficient for a sudden expansion %f.' %Kexp)


# For sudden contraction from hydrualic_Resistance.pdf

A1= #m^2 before cont

A3= #m^2 after

Kcont=(1/2)*(1-A3/A1)**(3/4) #where A3 is area of pipe after contraction and a fair ways down so the stream is fully developed again see vena contracta effect for more

print('This is the loss coefficient for a sudden contraction %f.' %Kcont)



Ktot=K90 + K180 + Kexp + Kcont

print('The loss coefficient for the whole loop is %f' %Ktot)








