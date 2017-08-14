# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:06:10 2017
@author: Parth
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:19:22 2017
@author: Parth
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 17:10:21 2017
@author: parth
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')
import numpy as np
import math
import time
import os
from matplotlib import pyplot as plt

number = 82
dt = 1E-3


#%% Arrays initialisation
T = np.zeros(number)
Db = np.zeros(number)
A = np.zeros(number)
v = np.zeros(number) #bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)
N = np.zeros(number) #chack
M = np.zeros(number)
Pb = np.zeros(number)
par = np.zeros(number)
RF = np.zeros(number)


#%%Constants


k = 1.38E-23 #Boltzman constant J/K
g = 9.8 # gravitational acceleration
P0 = 1E6 # Atmosperic pressure 1 bar = 1E5 pascals
H = 0.6 # height of the pool
mass_of_vap_sod = 1E-3 # mass of vaporised sodium in kg
Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 300 # sodium pool temperature
cp = 1284.4 #kJ/kg
h = 1E1 #check
lg = 4197E3
dp = 2E-6 #particle diameters
Na = 6.022E23 #avogrado number
mug = 1.7784896E-5 #gas viscosity (nitrogen)
mul = 8.90E-4 # liquid viscosity (water)




#%% Input Parameter Needed




#M[0] = 1E-1 #initial mass of aerosols in bubble
Pb[0] = 0.98E6 #check
T[0]= 300 #check
#rhol[0]=949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
#Pb[0] =  P0 + rhol[0]*g*(H)
Db[0] = 0.6
v[0] = 0 #check value


#Db[:] = 0.6 #meter
v[0] = 0
#Pb[:] = 1E6

#Db[0] = 1E-1 #initial size of bubble
#Db[0]= (N[0]*k*T[0]/(3.14*Pb[0]))**(1/3.0)
#rhog[0] = M[0]/(0.144*3.14*Db[0]**3) #density of the gas in bubble

#rhog[:] = 5E-2

#rhog[:]= 20E-3/(1.5*0.001)
#N[:] = rhog[0]*Na/56
rhol[:] = 1000 #kg/m3
print('sodium density \t{:0.3e}'.format(rhol[0]))
print('bubble diameter \t{:0.3e}'.format(Db[0]*1E3))
print('pool height \t',H)
print('bubble diameter',Db[0]*100)
print('particle diameter',dp)
#T[:] = 300
v[0] = 2.4*((math.pi/6)*Db[0]**3)**1.0/6
#rhog[0] = 20E-3/(1.5*0.001)
rhog[:] = Pb[0]*14/(8.314*T[0]) #P*MM/(R*T)
#rhop = 10/(1/6* (math.pi*Db[0]**3))
M[0] = rhog[0]*(1/6* (math.pi*Db[0]**3))

rhop=10E-3/(1/6* (math.pi*Db[0]**3))

for i in range(number-1):

#    print(i)

    t[i] = i*dt

    y[i] = (i)*dt*v[i]

    Pb[i] = P0 + rhol[i]*g*(H-y[i]) # pressure equation

#    Db[0]= (N[0]*k*T[0]/(3.14*Pb[0]))**(1/3.0)

    Db[:] = 0.6

#    N[i] = rhog[i]*Na/23

#    if i > 0 :
#        Db[i] = (N[i]*k*T[i]/(3.14*Pb[i]))**(1.0/3.0) # diameter equation


    T[i+1] = T[i] + dt* (0-h*3.14*Db[i]**2 * (T[i]-Ts)) # temperature equation


#    rhol[i+1] = 949-0.223*(T[i+1]-273.15)-1.7E-5*(T[i+1]-273.15)**2 #density of sodium

    v[i+1] = v[i] + dt*(6/(rhol[i]*3.14*Db[i]**3))*(  (3.14*Db[i]**3 *g*(rhol[i]-rhog[i]))/6 - 12*3.14*Db[i]*v[i]) #velocity of bubble


    v[i] = 2.4*((3.14/6)*Db[i]**3)**1.0/6

    lembda = 82.06*T[i]/(float(math.sqrt(2))*3.14*Na*Pb[i])*1E-2 #mean free path of a gas molecule

    cunn = 1 + (2*lembda/dp) * (1.257 + 0.4 * math.exp(-0.55*dp/lembda)) #cunningham slip correction


#    cunn = 1
    theta = k*T[i]*cunn/(3*math.pi*mug*dp)

    tau = rhop*dp**2 *cunn /(18*mug)


#    tau = 10E-3/(3*math.pi*mug*dp)
#    tau = 20E-3/(3*math.pi*mul*dp)
    diff = 1.8*(8*theta/(v[i]*Db[i]**3))**(1.0/2) #diffusion coefficient

    sedimentation = 1.5*g*tau/(Db[i]*v[i]) #sedimentation coefficient

    inertial_impaction = 18*v[i]*tau/Db[i]**2 # inertial impaction coefficient

    c = diff + sedimentation + inertial_impaction

#    dy = dt*v[i]
#    par[i] = dy*(-c*M[i])
#    M[i+1] = M[i] + par[i]

    RF[i] = math.exp(-y[i]*c)
#    rhog[i+1] = M[i+1]/(0.144*3.14*Db[i]**3 )

#    N[i+1] = rhog[i+1]*Na/235
#    if i >0:
#        M[i] =  M[0] * math.exp(-c*y[i])
#        rhog[i] = M[i]/(0.144*3.14*Db[i]**3 )


    if y[i]>H:
        break


Rf_factor = 1/RF[i]

print('inertial impaction \t ',(inertial_impaction))
print('sedimentation \t',sedimentation)
print('diffusion \t ',(diff))

'''
plt.figure()
plt.plot(t[0:i],Db[0:i],label='diameter',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('diameter (meter)',color='blue',size=14)
plt.legend()
plt.savefig('diameter')
plt.figure()
plt.plot(t[0:i],np.log(1/(RF[0:i])),label='Decontamination fraction',color='blue')
plt.xlabel('Time',color='blue',size=14)
plt.ylabel('Decontamination  fraction',color='blue',size=14)
plt.legend()
plt.savefig('Decontamination')

plt.figure()
plt.plot(t[0:i],(RF[0:i]),label='release fraction',color='blue')
plt.xlabel('Time',color='blue',size=14)
plt.ylabel('Release fraction',color='blue',size=14)
plt.legend()
plt.savefig('Release')
'''
'''
plt.figure()
plt.plot(t[:-4],T[:-4],label = 'Temperature',color='blue')
plt.xlabel('time',color='blue',size=14)
plt.ylabel('Temp (K)',color='blue',size=14)
plt.legend()
plt.savefig('temperature')
plt.figure()
plt.plot(t[:-4],Pb[:-4],label='pressure',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('Pressure (Pa)',color='blue',size=14)
plt.legend()
plt.savefig('pressure')
plt.figure()
plt.plot(t[:-4],v[:-4],label='velocity',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('velocity(m/s)',color='blue',size=14)

plt.legend()
plt.savefig('velocity')

'''



'''

#convert meter to cm
Db =(Db*1E2)
y = y*1E2-250
#convert array to list
Db = list(Db)
y = list(y)
import turtle

def draw_circle(turtle, color, size, x, y):
    turtle.penup()
    turtle.color(color)

    turtle.goto(x,y)
    turtle.pendown()

    turtle.circle(size)


turtle.shape("circle")
#turtle.resizemode("user")
turtle.shapesize(1E-2,1E-2)
#tommy.speed(7)
scale = 100
# Draw three circles:
for i in range(int(len(y)/scale)):
    draw_circle(turtle, "green", Db[scale*i], 25, y[scale*i])

'''





