import neuron
import numpy as np
import matplotlib.pyplot as plt
from neuronpp.cells.combe2018_cell import Combe2018Cell

# Create cell
cell = Combe2018Cell(name="cell")

soma = cell.filter_secs("soma")

fig, ax = plt.subplots(1, 1)
injections = [-.100, 0, .100, .200]
out = []
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

f1 = open("combe_orig_apic.txt", "w")
f2 = open("combe_orig_trunk.txt", "w")
f3 = open("combe_orig_dend.txt", "w")
f4 = open("combe_orig_sa.txt", "w")
for sec in cell.secs:
    if "trunk" in sec.hoc.name():
        f = f2
    elif "apic" in sec.hoc.name():
        f = f1
    elif "dend" in sec.hoc.name():
        f = f3
    else:
        f = f4

    mechs = sec.hoc.psection()

    for key in sorted(mechs.keys()):
        if isinstance(mechs[key], dict):
            for new_key in  sorted(mechs[key].keys()):
                f.write("%s %s %s " % (sec.hoc.name(), key, new_key) + str(mechs[key][new_key]) + "\n")
        else:
            f.write("%s %s " % (sec.hoc.name(), key) + str(mechs[key]) + "\n")

ax.set_xlabel("time (s)")
ax.set_ylabel("V (mV)")
ax.legend()
plt.show()
