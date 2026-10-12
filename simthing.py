import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

Mevel   = ct.Solution('Mevel2017.yaml')
h202    = ct.Solution('h2o2.yaml')
Alzueta = ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml')
Burke   = ct.Solution('Burke2012.yaml')

def statemachine(temp = 1000, pressure = ct.one_atm, i=1, solution = Burke):
    air = "O2:1,N2:3.76"

    gas = solution

    gas.set_equivalence_ratio(i, fuel="H2", oxidizer=air)
    gas.TP = temp, pressure
    r = ct.IdealGasReactor(gas)

    sim = ct.ReactorNet([r])
    sim.verbose = True

    # limit advance when temperature difference is exceeded
    delta_T_max = 20.
    r.set_advance_limit('temperature', delta_T_max)

    dt_max = 1.e-8
    states = ct.SolutionArray(gas, extra=['t'])

    print('{:10s} {:10s} {:10s} {:14s}'.format(
        't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))

    HOT = True
    ignition = False

    while HOT:
        sim.advance(sim.time + dt_max)
        states.append(r.thermo.state, t=sim.time * 1.e3)
        print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(sim.time, r.T, r.thermo.P, r.thermo.u))

        if sim.time > dt_max * 10 and np.max(np.gradient(states.T[-10:], states.t[-10:])) > 1000:
            ignition = True

        if ignition and np.max(np.gradient(states.T[-50:], states.t[-50:])) < 0.001:
            HOT = False

    return states

def ignition_delay(states):
    ID = states.t[np.argmax(np.gradient(states.T, states.t))]
    print(f"Ignition delay: {ID} ms")

    fig, ax = plt.subplots()
    ax.plot(states.t, states.T)
    ax.axvline(ID, color="red")
    #plt.show()

    return ID


    
#max(states.T) + min(states.T)) / 2

#a = ignition_delay(statemachine(1000, 6.e6))
b = ignition_delay(statemachine(1250, 6.e6, 1, Burke))
c = ignition_delay(statemachine(1250, 6.e6, 1, Mevel))
#c = ignition_delay(statemachine(1200, 6.e6))

print(b, c)
