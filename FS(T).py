
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


# Simulation parameters
p = [1, 40, 80] # pressure [Pa]
 # unburned gas temperature [K]

dt = 20
n = 10

air = "O2:1,N2:3.76"

width = 0.03  # m
loglevel = 0  # amount of diagnostic output (0 to 8)

#gas = ct.Solution('h2o2.yaml')
gas = ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml')

plt.clf()
plt.subplot(1, 1, 1)



for pressure in p:
    
    Tin = 700
    i = 0

    Fuel = []
    speeds = []

    while i < n:
        gas.set_equivalence_ratio(1, fuel="NH3", oxidizer=air)
        gas.TP = Tin, pressure * 10**5

        # Set up flame object
        f = ct.FreeFlame(gas)
        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.10)

        # Solve with mixture-averaged transport model
        f.transport_model = 'mixture-averaged'
        f.solve(loglevel=loglevel, auto=True)

        print("flame speed: " + str(f.velocity[0]))

        Fuel.append(Tin)
        speeds.append(f.velocity[0])

        Tin += dt
        i += 1

    plt.plot(Fuel, speeds, label=f"p = {pressure} bar")





plt.title("Laminar flame speed by equivalence ratio \n p = 1atm")
plt.xlabel('Temperature (K)')
plt.ylabel('Laminar flame speed (m/s)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()