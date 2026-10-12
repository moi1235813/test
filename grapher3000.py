import matplotlib.pyplot as plt
from simthing import statemachine
from simthing import ignition_delay
import cantera as ct


temp = 900
d_t = 30
no = 5
states = []
delays = []
initialtemps = []


i = 0
while i<no:
    initialtemps.append(temp)
    s = statemachine(temp, ct.one_atm*60)
    states.append(s)
    delays.append(ignition_delay(s))
    temp += d_t
    i += 1


print(delays)
plt.clf()
plt.subplot(1, 2, 1)

for i in states:
    plt.plot(i.t, i.T)

plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')

plt.subplot(1,2,2)
plt.plot(initialtemps, delays)
plt.xlabel('Temperature (K)')
plt.ylabel('Ignition delay (ms)')

plt.tight_layout()
plt.show()