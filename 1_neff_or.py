import PyMoosh as pm
import numpy as np
import matplotlib.pyplot as plt
import time
plt.rcParams['figure.dpi'] = 150

def permittivity_glass(wl):
    #epsilon=2.978645+0.008777808/(wl**2*1e-6-0.010609)+84.06224/(wl**2*1e-6-96)
    epsilon_glass = (1.5130 - 3.169e-9*wl**2 + 3.962e3/wl**2)**2
    return epsilon_glass

Glass = pm.Material(permittivity_glass)

Materials = ['Silver', 1, 'Gold', Glass]
Stack = [0, 1, 2, 3]

wavelength = 600

h_gap = np.linspace(10, 1, 10)
h_au = [49, 48, 47]

print(h_gap, h_au)


for j in h_au:
    n_eff = []
    n_start=3.1 + 0.19j
    
    for i in h_gap:
        Thicknesses = [0, i, j, 0]     # On retrouve pas la forme du bon mode pour < 48
        structure = pm.Structure(Materials, Stack, Thicknesses, verbose=False)   
        
        res = pm.steepest(n_start, 1e-12, 10000, structure, wavelength, 1)
        n_start = res
        n_eff.append(res)
       
    plt.plot(h_gap, np.real(n_eff), label='$h_{au} = $' + str(j) + ' nm')


plt.ylabel('$Re(n_{eff})$')
plt.xlabel('$h_{gap} (nm)$')
plt.legend()
plt.show()

print(n_eff)    #h_gap = 10nm => n_eff = 3.1 + 0.19j