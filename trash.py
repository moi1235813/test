from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 296  # unburned gas temperature [K]
n_states = 20


eq_ratios = np.linspace(0.2, 6, 20)
air = "O2:1,N2:3.76"

width = 0.03  # m
loglevel = 0  # amount of diagnostic output (0 to 8)

gas = ct.Solution('h2o2.yaml')
"""gas = ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml')"""

Fuel = []
speeds = []

for i in eq_ratios:
    gas.set_equivalence_ratio(i, fuel="H2", oxidizer=air)
    gas.TP = Tin, p

    # Set up flame object
    f = ct.FreeFlame(gas)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.10)

    # Solve with mixture-averaged transport model
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel=loglevel, auto=True)

    print("flame speed: " + str(f.velocity[0]))

    Fuel.append(i)
    speeds.append(f.velocity[0])


plt.clf()
plt.subplot(1, 1, 1)

plt.plot(Fuel, speeds)

plt.xlabel('φ')
plt.ylabel('Flame speed (m/s)')

plt.tight_layout()
plt.show()