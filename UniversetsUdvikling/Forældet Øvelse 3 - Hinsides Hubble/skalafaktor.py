#python-program work out the scale factor for different multi-component
#universes
#models are defined by Omega_m, Omega_L and Omage_0

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

plt.close('all')

def dtda(a):
    return 1./np.sqrt(omegam/a+omegal*a**2+(1.-omegan))

#Benchmark
omegam = 0.3
omegal = 0.7
omegan = omegam+omegal
a = 10**(np.arange(0,100,1)/40.-2.)
t = [0 for x in range(100)]
for n in range(100):
    I = spi.quad(dtda,np.min(a),a[n])
    t[n] = I[0]
tbm = t

#Flat, matterdominated
omegam = 1.0
omegal = 0.0
omegan = omegam+omegal
a = 10**(np.arange(0,100,1)/40.-2.)
t = [0 for x in range(100)]
for n in range(100):
    I = spi.quad(dtda,np.min(a),a[n])
    t[n] = I[0]
tf = t

#Loittering universe
omegam = 0.3370
omegal = 1.77010
omegan = omegam+omegal
a = 10**(np.arange(0,100,1)/40.-2.)
t = [0 for x in range(100)]
for n in range(100):
    I = spi.quad(dtda,np.min(a),a[n])
    t[n] = I[0]
tl = t

fig, ax = plt.subplots()
ax.plot(tbm, a, 'b', label='Benchmark')
ax.plot(tf, a, 'r', label='Einstein–de Sitter')
ax.plot(tl, a, 'g', label='Loitering')
ax.set(xlabel='H0*t', ylabel='Skalafaktor')
ax.grid()
ax.legend()
fig.savefig("test.png")
plt.show()

