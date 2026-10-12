import cantera as ct
import matplotlib.pyplot as plt
from simthing import ignition_delay

# 'h2o2.yaml'

fuels = ['example_data/ammonia-CO-H2-Alzueta-2023.yaml']

def statemachine(temp = 1000, pressure = ct.one_atm, i=1, fuel='h2o2.yaml'):
    air = "O2:1,N2:3.76"
    gas = ct.Solution(fuel)

    if fuel == 'h2o2.yaml':
        gas.set_equivalence_ratio(i, fuel="H2", oxidizer=air)
    else:
        gas.set_equivalence_ratio(i, fuel="H2", oxidizer=air)

    gas.TP = temp, pressure
    r = ct.IdealGasReactor(gas)

    sim = ct.ReactorNet([r])
    sim.verbose = True

    # limit advance when temperature difference is exceeded
    delta_T_max = 20.
    r.set_advance_limit('temperature', delta_T_max)

    dt_max = 1.e-4
    t_end = 80000 * dt_max
    states = ct.SolutionArray(gas, extra=['t'])

    print('{:10s} {:10s} {:10s} {:14s}'.format(
        't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
    while sim.time < t_end:
        sim.advance(sim.time + dt_max)
        states.append(r.thermo.state, t=sim.time * 1e3)
        print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
            sim.time, r.T, r.thermo.P, r.thermo.u))

    return states

n_temp_states = 10
d_t = 30

plt.subplot(1,1,1)

for thing in fuels:
    temp = 900
    pressure_bar = 60
    states = []
    delays = []
    initialtemps = []

    no = 0
    while no < n_temp_states:
        initialtemps.append(temp)
        s = statemachine(temp, pressure_bar * 10 ** 5, 1, thing)
        states.append(s)
        delays.append(ignition_delay(s))
        temp += d_t
        no += 1

    print(delays)

    if thing == "h2o2.yaml":
        plt.plot(initialtemps, delays, label="H2", linestyle="dashed")
    else:
        plt.plot(initialtemps, delays, label="NH3")


plt.xlabel('Temperature (K)')
plt.ylabel('ID (ms)')
plt.yscale("log")
plt.legend()
plt.grid(visible=True)

plt.title(f"Comparison between hydrogen's and ammonia's ID by temperature\n p = {pressure_bar} bar")

plt.tight_layout()
plt.show()