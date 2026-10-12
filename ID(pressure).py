
import matplotlib.pyplot as plt
from simthing import statemachine
from simthing import ignition_delay

plt.subplot(1,1,1)

d_p = 0.25
n_pressure_states = 35

temps = [950, 1000, 1050]

for temp in temps:
    pressure_bar = 1
    states = []
    delays = []
    initialp = []

    no = 0
    while no<n_pressure_states:
        initialp.append(pressure_bar)
        s = statemachine(temp, pressure_bar*10**5)
        states.append(s)
        delays.append(ignition_delay(s))
        pressure_bar += d_p
        no += 1

    print(delays)



    plt.plot(initialp, delays, label=f"T = {temp}K")



plt.xlabel('Pressure (bar)')
plt.ylabel('ID (s)')
#plt.yscale("log")
plt.legend()
plt.grid(visible=True)

plt.title(f"ID of hydrogen-air mixture as a function of pressure")

plt.tight_layout()
plt.show()