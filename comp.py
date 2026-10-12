
import matplotlib.pyplot as plt
from simthing import statemachine
from simthing import ignition_delay

plt.subplot(1,1,1)

d_t = 20
n_temp_states = 20

ER = [0.25, 1, 5]

for comp in ER:
    temp = 850
    pressure_bar = 60
    states = []
    delays = []
    initialtemps = []

    no = 0
    while no<n_temp_states:
        initialtemps.append(temp)
        s = statemachine(temp, pressure_bar*10**5, comp)
        states.append(s)
        delays.append(ignition_delay(s))
        temp += d_t
        no += 1

    print(delays)


    if comp == 1:
        plt.plot(initialtemps, delays, label=f"φ = {comp}", linestyle="dashed")
    else:
        plt.plot(initialtemps, delays, label=f"φ = {comp}")


plt.xlabel('Temperature (K)')
plt.ylabel('τ (ms)')
plt.yscale("log")
plt.legend()
plt.grid(visible=True)

plt.title(f"ID of hydrogen-air mixture as a function of the initial \n temperature at different equivalence ratios \n p = {pressure_bar} bar")

plt.tight_layout()
plt.show()