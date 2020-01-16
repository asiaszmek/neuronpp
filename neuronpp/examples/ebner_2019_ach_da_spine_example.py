from neuron import h
from neuron.units import mV
import numpy as np
import matplotlib.pyplot as plt

from neuronpp.cells.ebner2019_ach_da_cell import Ebner2019AChDACell
from neuronpp.core.cells.netstim_cell import NetStimCell
from neuronpp.core.cells.spine_cell import SpineCell
from neuronpp.core.cells.vecstim_cell import VecStimCell
from neuronpp.core.utils.record import Record
from neuronpp.core.utils.run_sim import RunSim


class Ebner2019AChDaSpineCell(Ebner2019AChDACell, SpineCell):
    def __init__(self, name):
        SpineCell.__init__(self, name)
        Ebner2019AChDACell.__init__(self, name)


WEIGHT = 0.0035		# µS, conductance of (single) synaptic potentials
WARMUP = 200


if __name__ == '__main__':
    h.load_file('stdrun.hoc')

    # define cell
    cell = Ebner2019AChDaSpineCell(name="cell")
    cell.load_morpho(filepath='morphologies/swc/my.swc', seg_per_L_um=1, add_const_segs=11)
    cell.add_spines(spine_number=10, head_nseg=10, neck_nseg=10, name='dend')
    cell.add_soma_mechanisms()
    cell.add_apical_mechanisms(sections='dend head neck')
    cell.add_4p_ach_da_synapse(point_process_name="head", loc=1)  # add synapse at the top of each spine's head

    # Create stims
    ns_cell = NetStimCell("netstim_cell")
    stim1 = ns_cell.add_netstim("stim1", start=WARMUP + 1, number=1)

    vs_cell = VecStimCell("vecstim_cell")
    stim2 = vs_cell.add_vecstim("stim3", np.array([WARMUP+50]))

    # stimulation
    cell.add_netcons(source=stim1, weight=WEIGHT, delay=1, mod_name="SynACh", name="head[0][0]")
    cell.add_netcons(source=stim2, weight=WEIGHT, delay=1, mod_name="SynDa", name="head[0][0]")
    # empty source for events
    ncons = cell.add_netcons(source=None, weight=WEIGHT, delay=1, mod_name="Syn4PAChDa", name="head[0][0]")

    # create plots
    rec_4psyn = Record(cell.filter_point_processes(mod_name="Syn4PAChDa", name="head[0][0]"), variables="w")

    # init and run
    h.finitialize(-70 * mV)
    sim = RunSim(warmup=WARMUP)
    sim.run(runtime=200)

    # Event delivery
    ncons[0].hoc.event(h.t+10)

    sim.run(runtime=200)

    # plot
    rec_4psyn.plot()
    plt.show()
