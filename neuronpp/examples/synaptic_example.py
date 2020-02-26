from neuronpp.utils.run_sim import RunSim

from neuronpp.cells.cell import Cell
from neuronpp.utils.record import Record
from neuronpp.core.cells.netstim_cell import NetStimCell


# Prepare cell
cell = Cell(name="cell", compile_paths="../commons/mods/ebner2019")
cell.load_morpho(filepath='../commons/morphologies/swc/my.swc')
cell.add_sec("dend[1]", diam=10, l=10, nseg=10)
cell.connect_secs(source="dend[1]", target="soma")
cell.insert("pas")
cell.insert("hh")

# Two examples of synapses with NetStim:
stim_cell = NetStimCell("stim_cell")
stim = stim_cell.make_netstim(start=250, number=30, interval=1)
soma = cell.filter_secs("soma")

# 1) Hoc-style synapse
pp = cell.add_point_process(mod_name="ExpSyn", seg=soma(0.5))
cell.add_netcon(source=stim, point_process=pp, netcon_weight=0.01, delay=1)

# 2) Recommended synapse
syn1 = cell.add_sypanse(source=stim, seg=soma(0.5), weight=0.01, mod_name="Syn4P", delay=1)

# 3) Event synapse
syn2 = cell.add_sypanse(source=None, seg=soma(0.5), weight=0.01, mod_name="Syn4P", delay=1)

# prepare plots
rec_v = Record(soma(0.5), variables="v")

# run
sim = RunSim(init_v=-55, warmup=20)

# Making external events to the synapse
syn2.make_event(10)
syn2.make_event(20)
syn2.make_event(30)
syn2.make_event(40)

# Example of online update of the graph, however sim.run(runtime=1000) can be call in a single run
for i in range(1000):
    sim.run(runtime=1)
    rec_v.plot(animate=True)
