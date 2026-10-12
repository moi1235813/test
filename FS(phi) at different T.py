from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from si import compute_adiabatic_T


# Simulation parameters
p = 6e6 # pressure [Pa]

Tfuel = 300
Tox = 900

eq_ratios = 0.15 #np.linspace(0.2, 5, 10)
air = "O2:1,N2:3.76"

width = 0.005  # m
loglevel = 2  # amount of diagnostic output (0 to 8)

#gas = ct.Solution('h2o2.yaml')
gas = ct.Solution('output2.yaml')
#gas = ct.Solution('Mevel2017.yaml')
#gas = ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml')

plt.clf()
plt.subplot(1, 1, 1)




gas.set_equivalence_ratio(eq_ratios, fuel="H2", oxidizer=air)

T_mix, X_target = compute_adiabatic_T(gas, T_oxidiser=Tox, P_pa=p, T_fuel=Tfuel)
gas.TPX = T_mix, p, X_target


# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.1, curve=0.2, prune=0.05)

f.transport_model = 'multicomponent'
f.solve(loglevel=loglevel, auto=True)

print("flame speed: " + str(f.velocity[0]))

#Fuel.append(eq_ratios)
#speeds.append(f.velocity[0])

grid = f.grid * 1000

plt.plot(grid, f.T, color='C3')

plt.title("Spatial temperature distribution")
plt.xlabel('x [mm]')
plt.ylabel('T [K]')
plt.grid()
plt.tight_layout()
plt.show()


dTdx = np.gradient(f.T, f.grid)
flameThic = (f.T[-1]-f.T[0]) / (np.max(dTdx))

print(f"Flame thickness: {flameThic}")

mass_fractions = f.Y

o2_idx = gas.species_index('O2')
h2o_idx = gas.species_index('H2O')
h2_idx = gas.species_index('H2')
n2_idx = gas.species_index('N2')

plt.figure(figsize=(8, 5))
plt.plot(grid, mass_fractions[o2_idx, :], label='$Y_{O2}$', lw=2)
plt.plot(grid, mass_fractions[h2o_idx, :], label='$Y_{H2O}$', lw=2)
plt.plot(grid, mass_fractions[h2_idx, :], label='$Y_{H2}$', lw=2)
plt.plot(grid, mass_fractions[n2_idx, :], label='$Y_{N2}$', lw=2)
#plt.plot(Fuel, speeds, label=f"T={Tin} K")


#print(Fuel, speeds)

plt.title("Specie mass fractions")
plt.xlabel('x [mm]')
plt.ylabel('Y')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
