import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

T_fuel = 400
T_air = 900

U = 8
D = 1e-4

p = 6e6

def compute_adiabatic_T(gas, T_oxidiser, P_pa, T_fuel):
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

air = "O2:1,N2:3.76"

Zs = np.linspace(0, 1, 100)

gas = ct.Solution('output2.yaml')
gas.set_equivalence_ratio(1, fuel="H2", oxidizer=air)

Z_st = Z_st = gas.Y[gas.species_index("H2")]

data = []

for Z in Zs:
    gas = ct.Solution('output2.yaml')
    gas.set_mixture_fraction(Z, fuel="H2", oxidizer=air)

    T_mix, X_target = compute_adiabatic_T(gas, T_oxidiser=T_air, P_pa=p, T_fuel=T_fuel)
    gas.TPX = T_mix, p, X_target

    Re = gas.density*U*D/gas.viscosity
    Ma = U / gas.sound_speed

    data.append([Z, Re, Ma])

    print(f"Z: {Z:.2f}, Re = {Re:.2f}, Ma = {Ma:.4f}")

data = np.array(data)

plt.clf()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(Zs, data[:, 1], color="black", label="Reynolds number")
ax2.plot(Zs, data[:, 2], label="Mach number", color="red")


ax1.axvline(x=Z_st, color="gray", linestyle="dotted")

plt.title(f"Re and Ma vs mixture fraction")
ax1.set_ylabel('Re')
ax1.set(xlabel="Z")
ax2.set_ylabel('Ma')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

ax1.grid()

plt.show()