from PyMoosh import *
import numpy as np
import matplotlib.pyplot as plt

def permittivity_glass(wl):
    #epsilon=2.978645+0.008777808/(wl**2*1e-6-0.010609)+84.06224/(wl**2*1e-6-96)
    epsilon_glass = (1.5130 - 3.169e-9*wl**2 + 3.962e3/wl**2)**2
    return epsilon_glass
Glass = Material(permittivity_glass)

Wavelength = 600

h_gap = 5
h_au = 10

Materials = ['Silver', 1, 'Gold', Glass]
Stack = [0, 1, 2, 3]
Thicknesses = [0, h_gap, h_au, 0]     
     
S = Structure(Materials, Stack, Thicknesses)

'''mode = guided_modes(S, Wavelength, 1, 5+1j, 10+1j, initial_points=100)
print(mode)'''

n_eff = 7.110943204110561+0.7127984342605717j

S.thickness[0] = 20
S.thickness[3] = 20
x, prof = profile(S, n_eff, Wavelength, 1, pixel_size=0.1)
#print(x, prof)

plt.plot(x, np.real(prof), linewidth = 2)
plt.ylabel('Re(H)')
plt.xlabel('z (nm)')

Sum = 0
for i in S.thickness[:-1]: 
    #plt.text(S, min(y), str(Materials[Stack[i]]) +'\n('+ str(Thicknesses[Stack[i]]) +' nm)', rotation='vertical', size='x-small') 
    Sum += i
    plt.axvline(Sum, color='red', alpha=0.2)
plt.show()