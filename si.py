import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

Mevel   = ct.Solution('Mevel2017.yaml')
h202    = ct.Solution('h2o2.yaml')
Alzueta = ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml')
Burke   = ct.Solution('Burke2012.yaml')
Burke2  = ct.Solution('output2.yaml')

T_FUEL = 300.0  # K - cold H2 fuel inlet temperature (ref: Yip et al.)


def compute_adiabatic_T(gas, T_oxidiser, P_pa, T_fuel=T_FUEL):
    """
    Compute the adiabatic mixture temperature for a PURE H2 fuel case.

    WHY THIS EXISTS:
    In the real engine, H2 is injected COLD (T_fuel ~380 K) into HOT
    compressed air (T_oxidiser). The instant they mix - before any
    chemistry happens - the local temperature is not  T_oxidiser. It depends on how much cold fuel vs hot air is
    present, because the fuel physically cools whatever air it mixes
    into. This mixes ENTHALPY (not temperature directly) between two
    pure streams, since enthalpy - not T - is what's actually
    conserved when two streams combine adiabatically.

    ASSUMES: gas's current mole-fraction composition (gas.X) already
    represents the FULLY MIXED target composition - i.e. call this
    right after gas.set_equivalence_ratio(...), before overwriting
    gas.TP with a uniform temperature.

    Parameters
    ----------
    gas         : ct.Solution - Cantera gas object, already set to the
                  target fully-mixed composition via
                  set_equivalence_ratio()
    T_oxidiser  : float - hot ambient/air temperature [K] (this is
                  what the script previously called 'temp')
    P_pa        : float - pressure [Pa]
    T_fuel      : float - cold H2 fuel inlet temperature [K]

    Returns
    -------
    T_mix   : float - adiabatic mixture temperature [K]
    X_target: dict/array - the original fully-mixed mole-fraction
              composition (returned so the caller can restore it,
              since this function temporarily overwrites gas.TPX
              several times while probing the pure streams)
    """

    # Save the fully-mixed TARGET composition before we start
    # overwriting gas.TPX to probe the pure fuel/air streams below.
    X_target = gas.X.copy()

    # ------------------------------------------------------------
    # STEP 1: Pure fuel stream (100% H2) enthalpy at cold T_fuel
    # ------------------------------------------------------------
    gas.TPX = T_fuel, P_pa, {'H2': 1.0}
    h_fuel  = gas.enthalpy_mass   # specific enthalpy of pure cold H2 [J/kg]
    Y_fuel  = gas.Y.copy()

    # ------------------------------------------------------------
    # STEP 2: Pure oxidiser stream (air) enthalpy at hot T_oxidiser
    # ------------------------------------------------------------
    gas.TPX = T_oxidiser, P_pa, {'O2': 0.21, 'N2': 3.76 / 4.76}
    h_ox    = gas.enthalpy_mass   # specific enthalpy of pure hot air [J/kg]
    Y_ox    = gas.Y.copy()

    # ------------------------------------------------------------
    # STEP 3: Mixture fraction xi = mass fraction of the final
    # mixture that originated from the fuel stream.
    #
    # Set gas to the fully-mixed TARGET composition (at an arbitrary
    # reference T - mole/mass fraction ratios don't depend on T) and
    # read off the H2 mass fraction directly. That is xi, since fuel
    # here is pure H2.
    # ------------------------------------------------------------
    gas.TPX    = T_oxidiser, P_pa, X_target
    Y_mix_mole = gas.Y.copy()

    idx_H2 = gas.species_index('H2')
    xi     = float(np.clip(Y_mix_mole[idx_H2], 0.0, 1.0))

    # ------------------------------------------------------------
    # STEP 4: Adiabatic enthalpy mixing - mix h and Y linearly in xi.
    # This is the physical core: enthalpy mixes linearly, temperature
    # does not (heat capacity differs between streams).
    # ------------------------------------------------------------
    h_mix = xi * h_fuel + (1.0 - xi) * h_ox
    Y_mix = xi * Y_fuel + (1.0 - xi) * Y_ox
    Y_mix = Y_mix / Y_mix.sum()   # renormalize against float drift

    # ------------------------------------------------------------
    # STEP 5: Solve for the temperature consistent with h_mix at
    # this pressure and composition. This T_mix is the physically
    # correct "just-mixed, not-yet-reacted" starting temperature.
    # ------------------------------------------------------------
    gas.HPY = h_mix, P_pa, Y_mix
    T_mix   = gas.T

    return float(T_mix), X_target


def statemachine(temp=1000, pressure=ct.one_atm, i=1, solution=Burke):
    air = "O2:1,N2:3.76"

    gas = solution

    gas.set_equivalence_ratio(i, fuel="H2", oxidizer=air)

    # ------------------------------------------------------------
    # ADIABATIC MIXING LINE INITIALISATION
    #
    # Previously: gas.TP = temp, pressure  <- assumed fuel+air were
    # ALREADY mixed and equilibrated to a single uniform temperature
    # 'temp'. That ignores the cooling effect of injecting H2 cold
    # (T_FUEL) into hot air (temp)
    #
    # Now: 'temp' is treated as T_oxidiser (the hot air temperature),
    # and the reactor is initialised at T_mix, the physically correct
    # adiabatically-mixed temperature for this composition.
    # ------------------------------------------------------------
    T_mix, X_target = compute_adiabatic_T(gas, T_oxidiser=temp, P_pa=pressure, T_fuel=T_FUEL)
    gas.TPX = T_mix, pressure, X_target

    r = ct.IdealGasReactor(gas)

    sim = ct.ReactorNet([r])
    sim.rtol = 1e-10
    sim.atol = 1e-16
    sim.verbose = True

    # limit advance when temperature difference is exceeded
    delta_T_max = 20.
    r.set_advance_limit('temperature', delta_T_max)

    dt_max = 1.e-2
    states = ct.SolutionArray(gas, extra=['t'])

    print('{:10s} {:10s} {:10s} {:14s}'.format(
        't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
    print(f'(T_oxidiser={temp:.1f} K -> T_mix={T_mix:.1f} K, T_fuel={T_FUEL:.1f} K)')

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

"""a = ignition_delay(statemachine(900, 6.e6, 0.1, Burke2))
b = ignition_delay(statemachine(900, 6.e6, 0.2, Burke2))
c = ignition_delay(statemachine(900, 6.e6, 0.3, Burke2))"""
d = ignition_delay(statemachine(900, 6.e6, 0.15, Burke2))
#e = ignition_delay(statemachine(900, 6.e6, 0.7, Burke2))
#f = ignition_delay(statemachine(900, 6.e6, 1, Burke2))
#c = ignition_delay(statemachine(1200, 6.e6))

print(d)