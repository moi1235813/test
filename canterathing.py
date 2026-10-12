import matplotlib.pyplot as plt
from simthing import statemachine
from simthing import ignition_delay

plt.subplot(1,1,1)

d_t = 30
n_temp_states = 15
pressure_list = [40, 80, 150]



for p in pressure_list:
    temp = 900
    states = []
    delays = []
    initialtemps = []

    i = 0
    while i<n_temp_states:
        initialtemps.append(temp)
        s = statemachine(temp, p*10**5, 1)
        states.append(s)
        delays.append(ignition_delay(s))
        temp += d_t
        i += 1

    print(delays)


    plt.plot(initialtemps, delays, label=f"P={p} bar")


plt.xlabel('Temperature (K)')
plt.ylabel('Ignition delay (ms)')
plt.yscale("log")
plt.grid()
plt.legend()
plt.title("Hydrogen-air ignition delay \n φ=1")
plt.tight_layout()
plt.show()