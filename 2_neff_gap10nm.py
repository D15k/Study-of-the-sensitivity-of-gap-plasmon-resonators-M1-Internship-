from PyMoosh import *
import numpy as np
import matplotlib.pyplot as plt
import time
plt.rcParams['figure.dpi'] = 150


def steepest(start, tol, step_max, struct, wl, pol, steepest_step=1e-2):
    """ Steepest descent to find a zero of the `dispersion`
    function. The advantage of looking for a zero is that you
    know when the algorithm can stop (when the value of the function
    is smaller than `tol`).

    Args:
        start (complex): effective index where the descent starts
        tol (real): when dispersion is smaller than tol, the descent stops
        step_max (integer): maximum number of steps allowed (0 -> infinite)
        struct (Structure): the object describing the multilayer
        wl (float): wavelength in vacuum
        pol: 0 for TE, 1 for TM
        steepest_step: step of the steepest descent, 1e-2 if unscpecified

    Returns:

        (float) : last effective index reached at the end of the descent

    """


    k_0 = 2 * np.pi / wl
    z = start*k_0
    delta = abs(z) * 0.001
    dz = steepest_step * delta
    step = 0
    current = dispersion(z,struct,wl,pol)

    while (current > tol) and ((step < step_max) or (step_max == 0)):

        grad = (
        dispersion(z+dz,struct,wl,pol)
        -current
#        -dispersion(z-dz,struct,wl,pol)
        +1j*(dispersion(z+1j*dz,struct,wl,pol)
#        -dispersion(z-1j*dz,struct,wl,pol))
        -current)
        )/(dz)

        if abs(grad)!=0 :
            z_new = z - delta * grad / abs(grad)
        else:
            # We have a finishing condition not linked to the gradient
            # So if we meet a gradient of 0, we just divide the step by two
            delta = delta/2.
            z_new = z

        value_new = dispersion(z_new,struct,wl,pol)
        if (value_new > current):
            # The path not taken
            delta = delta / 2.
            dz = dz / 2.
        else:
            current = value_new
            z = z_new
    #        print("Step", step, z,current)
        step = step + 1

    #print("End of the loop")
    if step == step_max:
        print("Warning: maximum number of steps reached")

    return z/k_0


def permittivity_glass(wl):
    #epsilon=2.978645+0.008777808/(wl**2*1e-6-0.010609)+84.06224/(wl**2*1e-6-96)
    epsilon_glass = (1.5130 - 3.169e-9*wl**2 + 3.962e3/wl**2)**2
    return epsilon_glass

Glass = Material(permittivity_glass)

Materials = ['Silver', 1, 'Gold', Glass]
Stack = [0, 1, 2, 3]

wavelength = 600

nbr_pts = [50, 100, 150]
    
for j in nbr_pts:
    h = np.linspace(50, 1, j)
    
    print(h)
    
    n_eff = []
    n_start=2
    # Faut commencer par les valeurs grandes pour attraper le bon mode pr les prochaines descentes de gradiant
    # Plus le nombre de points est grand, moins on a de chance de changer de mode  .
    
    for i in h:
        Thicknesses = [0, 9, i, 0]     
        structure = Structure(Materials, Stack, Thicknesses, verbose=False)   
        
        res = steepest(n_start, 1e-12, 10000, structure, wavelength, 1, 1e-2)
        n_start = res
        n_eff.append(res)

    plt.plot(h, np.real(n_eff), label=str(j) + ' points')


plt.ylabel('$Re(n_{eff})$')
plt.xlabel('$h_{Au} (nm)$')
plt.title('$h_{Gap} = 9 \ nm$')
plt.legend()
plt.show()

print(n_eff)    