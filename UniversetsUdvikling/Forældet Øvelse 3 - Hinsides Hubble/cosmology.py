#python-program to work out the angular-distance, luminosity distance etc.
#world models with Lamda different from 0. I follow Longair's book :
#Galaxy-formation

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

plt.close('all')

def drdz(z):
    return cH/np.sqrt((1+z)**2*(Omega*z+1)-OL*z*(z+2))

def r(z):
    I = spi.quad(drdz,0,z)
    return I[0]

def D(z):
   if (kappa < 0.): 
      Rtilt = 1/np.sqrt(-kappa)
      res=Rtilt*np.sinh(r(z)/Rtilt)
   if kappa > 0.:
      Rtilt = 1/sqrt(kappa)
      res = Rtilt*np.sin(r(z)/Rtilt)
   if kappa == 0.: 
      res = r(z)
   return res

def AngDist(z):
   return D(z)/(1+z)

def LumDist(z):
   return D(z)*(1+z)

def dV(z):
   return D(z)**2*drdz(z)

def drdzprop(z):
   return drdz(z)/(1+z)

def time(z):
   I = spi.quad(drdzprop,np.inf,z)
   return -I[0]/(3.0e10/3.0857e24)/(3600.*24.*365.25)

def recession(z):
   return np.sqrt(Omega*(1+z)+OL/(1+z)**3)*D(z)/cH

def recessionnow(z):
   return np.sqrt(Omega+OL)*D(z)/cH

OL=0.69
Omega=1.-OL
cH = 2998./0.68
kappa = (Omega+OL-1.)/cH**2.

print('Time back to Big Bang:',time(0)/1.e9,'Gyr')
print('Time back to the Big Bang at z=1.72:',time(1.72)/1.e9, 'Gyr')
print('Look back time to z=1.72:',(time(0)-time(1.72))/1.e9, 'Gyr')
print('Time from z=3.6 to 5.0:',(time(3.6)-time(5.0))/1.e9, 'Gyr')
print('Time from z=1.6 to 2.0:',(time(1.6)-time(2.0))/1.e9, 'Gyr')
print('Time from z=0 to 0.42:',(time(0)-time(0.42))/1.e9, 'Gyr')
print('Recession velocity at z=10',recession(10))
print('Recession velocity at z=10 now',recessionnow(10))
print('Kpc per arcsec at z=1',AngDist(1.)*2.*np.pi/(3600.*360.)*1000.)
