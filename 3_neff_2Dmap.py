from PyMoosh import *
import numpy as np
import math
import xarray as xr
import matplotlib.pyplot as plt

def permittivity_glass(wl):
    #epsilon=2.978645+0.008777808/(wl**2*1e-6-0.010609)+84.06224/(wl**2*1e-6-96)
    epsilon_glass = (1.5130 - 3.169e-9*wl**2 + 3.962e3/wl**2)**2
    return epsilon_glass
Glass = Material(permittivity_glass)

Wavelength = 600

h_gap = np.linspace(10, 1, 10)
h_au = np.linspace(10, 2, 25)

#h_au = np.linspace(50, 2, 49) 
#h_au = np.linspace(50, 2, 193) 
#h_au = np.append(np.linspace(50, 3.5, 94), np.linspace(3.4, 2, 15))
#h_au = np.linspace(50, 2, 385) 

#base = 1.1
#h_au = np.logspace(math.log(50, base), math.log(2, base), 200, base=base)
#print(h_au)

n_eff = xr.DataArray(data = np.zeros ((len(h_gap), len(h_au)), dtype = complex), 
                     coords = {'h_Gap': h_gap, 
                               'h_Au': h_au})
 

Materials = ['Silver', 1, 'Gold', Glass]
Stack = [0, 1, 2, 3]
Thicknesses = [0, 999, 999, 0]          # 999 => value is gonna be set correctly below  
S = Structure(Materials, Stack, Thicknesses)
                  
for i in range(len(h_gap)):
    S.thickness[1] = h_gap[i]           # float(n_eff.coords['h_gap'][i])
                   
    for j in range(len(h_au)): 
        S.thickness[2] = h_au[j] 
        #S.plot_stack()
        
        if   j != 0: n_start = complex(n_eff[ i , j-1])
        elif i != 0: n_start = complex(n_eff[i-1,  0 ])
        else: n_start = 5+0j
        
        n_eff[i, j] = steepest(n_start, 1e-12, 1000, S, Wavelength, 1)
        print('Steepest Done', (h_gap[i], h_au[j]))


print(complex(n_eff.isel(h_Gap=5, h_Au=0))) #n_eff for h_gap = 5nm and h_au = 10nm

#print(n_eff)
n_eff = n_eff.real
#print(n_eff)

n_eff.plot()
plt.title('$n_{eff}$')
plt.xlabel('$h_{Au}$')
plt.ylabel('$h_{Gap}$')
plt.show()

n_eff.plot.line(x = 'h_Au')
plt.xlabel('$h_{Au}$')
plt.ylabel('$n_{eff}$')
plt.show()