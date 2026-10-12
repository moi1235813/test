import matplotlib.pyplot as plt
import cantera as ct
from si import compute_adiabatic_T


p = 6e6  # pressure
Tfuel = 300
Tox = 900
eq_ratios = 1 #np.linspace(0.2, 5, 10)
air = "O2:1,N2:3.76"
tsurf = 900

# each mdot value will be solved to convergence, with grid refinement, and
# then that solution will be used for the next mdot
mdot = [150, 200]  # kg/m^2/s

rxnmech = 'output2.yaml'  # reaction mechanism file


# The solution domain is chosen to be 20 cm
width = 0.005  # m

loglevel = 1  # amount of diagnostic output (0 to 5)

# Grid refinement parameters
ratio = 3
slope = 0.1
curve = 0.2
prune = 0.06

gas = ct.Solution(rxnmech)

gas.set_equivalence_ratio(eq_ratios, fuel="H2", oxidizer=air)

T_mix, X_target = compute_adiabatic_T(gas, T_oxidiser=Tox, P_pa=p, T_fuel=Tfuel)
gas.TPX = T_mix, p, X_target


# Create the stagnation flow object with a non-reactive surface.  (To make the
# surface reactive, supply a surface reaction mechanism. See example
# catalytic_combustion.py for how to do this.)
sim = ct.ImpingingJet(gas=gas, width=width)

# set the mass flow rate at the inlet
sim.inlet.mdot = mdot[0]

# set the surface state
sim.surface.T = tsurf

sim.set_grid_min(1e-6)
sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
sim.show()

sim.solve(loglevel, auto=True)

results = []
for m, md in enumerate(mdot):
    sim.inlet.mdot = md
    sim.solve(loglevel)
    results.append(sim.to_array())


sim.show_stats()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
styles = ['-', '--']
for n, i in enumerate([0, -1]):
    ax1.plot(results[i].grid, results[i].velocity,
             linestyle=styles[n], color='C4', label=rf'$\dot m = {mdot[i]:.2f}$ kg/s/m²')
    ax2.plot(results[i].grid, results[i].T,
             linestyle=styles[n], color='C3', label=rf'$\dot m = {mdot[i]:.2f}$ kg/s/m²')

ax1.set_ylabel('heat release rate [MW/m³]', color='C4')
ax1.set(xlabel='distance from inlet [m]')
ax2.set_ylabel('temperature [K]', color='C3')
ax2.legend()
plt.show()