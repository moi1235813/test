import matplotlib.pyplot as plt
from simthing import statemachine
from simthing import ignition_delay
import numpy as np


temp = 1000
pressure = 1*10**5
eq_ratios = np.linspace(0.2, 6, 20)
delays = []

for phi in eq_ratios:
    s = statemachine(temp, pressure, phi)
    delays.append(ignition_delay(s))



print(delays)


plt.subplot(1,1, 1)
plt.plot(eq_ratios, delays)
plt.xlabel('ER')
plt.ylabel('Ignition delay (ms)')

plt.tight_layout()
plt.show()