import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

MECH   = 'output2.yaml'        # mechanism used for the sweep (was Burke2)
AIR    = "O2:1,N2:3.76"
T_FUEL = 300.0                 # K - cold H2 fuel inlet temperature (flame-speed script used 400 K)


def compute_adiabatic_T(gas, T_oxidiser, P_pa, T_fuel=T_FUEL):
    """
    Adiabatic (frozen, unreacted) mixing temperature of cold pure H2 with hot air.
    gas must already hold the fully-mixed target composition
    (call right after gas.set_equivalence_ratio(...)).
    Returns (T_mix, X_target).
    """
    X_target = gas.X.copy()

    # pure fuel stream
    gas.TPX = T_fuel, P_pa, {'H2': 1.0}
    h_fuel, Y_fuel = gas.enthalpy_mass, gas.Y.copy()

    # pure oxidiser stream
    gas.TPX = T_oxidiser, P_pa, {'O2': 0.21, 'N2': 3.76 / 4.76}
    h_ox, Y_ox = gas.enthalpy_mass, gas.Y.copy()

    # mixture fraction = H2 mass fraction of the mixed state (fuel is pure H2)
    gas.TPX = T_oxidiser, P_pa, X_target
    xi = float(np.clip(gas.Y[gas.species_index('H2')], 0.0, 1.0))

    # enthalpy mixes linearly in xi, temperature does not
    h_mix = xi * h_fuel + (1.0 - xi) * h_ox
    Y_mix = xi * Y_fuel + (1.0 - xi) * Y_ox
    Y_mix /= Y_mix.sum()

    gas.HPY = h_mix, P_pa, Y_mix
    return float(gas.T), X_target


def _run_reactor(gas, T_mix, X_target, pressure, dt_max, t_max, T_ign):
    """
    One constant-volume reactor run with a fixed output step dt_max [s].
    'Ignited' = T has passed T_ign (halfway between T_mix and equilibrium T).
    """
    gas.TPX = T_mix, pressure, X_target   # reset to the unreacted mixed state

    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])
    sim.rtol = 1e-10
    sim.atol = 1e-16

    r.set_advance_limit('temperature', 20.0)   # resolve the ignition event

    states = ct.SolutionArray(gas, extra=['t'])   # t stored in ms
    ignited = False

    while sim.time < t_max:
        sim.advance(sim.time + dt_max)
        states.append(r.thermo.state, t=sim.time * 1.e3)

        if not ignited and r.T > T_ign:
            ignited = True

        if ignited and len(states.t) > 50:
            if np.max(np.abs(np.gradient(states.T[-50:], states.t[-50:]))) < 1e-3:
                break

    return states, ignited


def statemachine(temp=1000, pressure=ct.one_atm, fraction=0.02, solution=None,
                 T_fuel=T_FUEL, t_max=20.0, dt_max=1.e-2,
                 min_pts=50, max_refine=6, verbose=False):
    """
    Constant-volume ignition run started from the adiabatic mixing temperature.
    temp is the hot oxidiser temperature.

    dt_max is reduced automatically: if fewer than min_pts output points fall
    before the ignition delay, the run is repeated with
        dt_max = min(tau / min_pts, dt_max / 2)
    up to max_refine times.

    Returns (states, ignited, T_mix).
    """
    gas = solution if solution is not None else ct.Solution(MECH)
    gas.set_mixture_fraction(fraction, fuel="H2", oxidizer=AIR)

    T_mix, X_target = compute_adiabatic_T(gas, T_oxidiser=temp, P_pa=pressure, T_fuel=T_fuel)

    # Ignition threshold relative to the available temperature rise, so very
    # lean mixtures (small rise) are not missed by a fixed +400 K criterion.
    gas.TPX = T_mix, pressure, X_target
    gas.equilibrate('UV')                 # constant-volume reactor -> UV equilibrium
    T_eq = gas.T
    T_ign = T_mix + 0.5 * (T_eq - T_mix)

    for attempt in range(max_refine + 1):
        states, ignited = _run_reactor(gas, T_mix, X_target, pressure, dt_max, t_max, T_ign)

        if not ignited:
            break   # nothing to resolve - a finer step won't help

        tau_ms = states.t[np.argmax(np.gradient(states.T, states.t))]
        n_before = int(np.sum(states.t < tau_ms))

        if verbose:
            print(f'  Z={fraction:.3f} attempt {attempt}: dt_max={dt_max:.2e} s, '
                  f'tau={tau_ms:.4e} ms, points before ignition={n_before}')

        if n_before >= min_pts:
            break

        dt_max = min(tau_ms * 1e-3 / min_pts, dt_max / 2.0)
    else:
        if verbose:
            print(f'  Z={fraction:.3f}: max_refine reached, tau may be under-resolved')

    if verbose:
        print(f'Z={fraction:.3f}: T_ox={temp:.1f} K -> T_mix={T_mix:.1f} K, '
              f'ignited={ignited}, dt_max={dt_max:.2e} s')

    return states, ignited, T_mix


def ignition_delay(states, ignited=True, plot=False):
    """Ignition delay [ms] = time of maximum dT/dt. NaN if no ignition."""
    if not ignited:
        return np.nan

    ID = states.t[np.argmax(np.gradient(states.T, states.t))]

    if plot:
        fig, ax = plt.subplots()
        ax.plot(states.t, states.T)
        ax.axvline(ID, color="red")
        ax.set(xlabel="t [ms]", ylabel="T [K]")

    return ID


if __name__ == "__main__":
    # Guarded so that `from si import compute_adiabatic_T` in the flame-speed
    # script no longer launches this sweep on import.

    p    = 6e6                         # Pa
    T_ox = 900.0                       # K
    # log-spaced so the lean side (where the shortest delays are) is resolved
    fractions = np.geomspace(0.0005, 0.08, 50)

    gas = ct.Solution(MECH)

    # stoichiometric mixture fraction, computed directly
    gas.set_equivalence_ratio(1.0, fuel="H2", oxidizer=AIR)
    Z_st = gas.Y[gas.species_index('H2')]

    # fuel temperatures to compare [K] - one curve each
    T_fuels = [300.0, 400.0, 500.0]

    results = {}   # T_fuel -> array with columns: phi, Z, T_mix, tau_ign [ms]
    for T_fuel in T_fuels:
        data = []
        for fraction in fractions:
            states, ignited, T_mix = statemachine(T_ox, p, fraction, gas,
                                                  T_fuel=T_fuel, verbose=True)
            tau = ignition_delay(states, ignited)

            #gas.set_equivalence_ratio(phi, fuel="H2", oxidizer=AIR)
            #Z = gas.Y[gas.species_index('H2')]

            data.append([fraction, T_mix, tau])
            with open("ignition_sweep.txt", "a") as file:
                file.write(f"[{T_fuel}, {fraction}, {T_mix}, {tau}]\n")

        data = np.array(data)
        results[T_fuel] = data

        no_ign = np.isnan(data[:, 2])
        if no_ign.any():
            print(f"T_fuel={T_fuel:.0f} K: no ignition within t_max for "
                  f"phi = {np.round(data[no_ign, 0], 2)}")

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()

    for k, (T_fuel, data) in enumerate(results.items()):
        ax1.semilogy(data[:, 0], data[:, 2], '-', color=f'C{k}',
                     label=rf"$T_{{fuel}}$={T_fuel:.0f} K")
        ax2.plot(data[:,0], data[:,1], linestyle="dashed")

    ax1.axvline(x=Z_st, color="gray", linestyle="dotted")
    ax1.set(xlabel="φ", ylabel=r"$\tau_{ign}$ [ms]")
    ax1.legend(loc='upper right', fontsize=8)

    ax1.grid(which='both', alpha=0.4)
    plt.title(f"Ignition delay vs H2 mass fraction (T_air = {T_ox:.0f} K, p = 60bar)")
    plt.tight_layout()
    plt.show()