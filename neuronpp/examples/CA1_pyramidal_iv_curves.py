import numpy as np
import matplotlib.pyplot as plt
from neuronpp.cells.CA1_pyramidal import CA1PyramidalCell
import neuron


# Create cell
neuron.h.CVode()
neuron.h.CVode().active(True)

cell = CA1PyramidalCell(name="cell1")
soma = cell.filter_secs("soma")

fig, ax = plt.subplots(1, 1)
injections = [-.2, 0, .200, .400]
neuron.h.celsius = 34

for inj in injections:
    print(inj)
    ic = neuron.h.IClamp(soma.hoc(0.5))
    rec_v = neuron.h.Vector().record(soma.hoc(0.5)._ref_v)
    time = neuron.h.Vector().record(neuron.h._ref_t)
    ic.delay = 300
    ic.dur = 1000
    ic.amp = inj
    neuron.h.finitialize(-70)
    neuron.h.fcurrent()
    neuron.run(1500)
    ax.plot(time, rec_v, label="%4.3fnA" % inj)


 
ax.set_xlabel("time (s)")
ax.set_ylabel("V (mV)")
ax.legend()
plt.show()
