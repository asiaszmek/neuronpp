from neuronpp.core.cells.netstim_cell import NetStimCell

from neuronpp.cells.hoc_cell import HocCell
from neuronpp.utils.record import Record
from neuronpp.utils.iclamp import IClamp
from neuronpp.utils.run_sim import RunSim

REPS = 5  # Number of pre-and postsynaptic spikes
DT = 0.025
AMP = 2.7
DUR = 5.0
WARMUP = 2300
COOL_DOWN = 100
WEIGHT = 0.0035

freq = 50  # freqs.append(10, 20, 40, 50)
interval = 1000 / freq
delta_t = 10  # LTP

if __name__ == '__main__':
    # Cell definition
    cell = HocCell("cell", compile_paths="../commons/mods/hay2011 ../commons/mods/ebner2019")
    cell.load_hoc("../commons/hocmodels/ebner2019/load_model.hoc", cell_template_name="L5PCtemplate")
    soma = cell.filter_secs("soma")[0]

    # Netstim to synapse
    stim = NetStimCell("stim")
    st = stim.make_netstim(start=WARMUP, number=REPS, interval=interval)
    cell.make_netcons(source=st, weight=WEIGHT, mod_name="Syn4P")

    # IClamp to soma
    iclamp = IClamp(segment=soma.hoc(0.5))
    start_t = WARMUP + delta_t
    for i in range(REPS):
        start_t += i * interval
        iclamp.stim(delay=start_t, dur=DUR, amp=AMP)

    # Record
    rec = Record(cell.filter_secs("apic[1],apic[50]"), loc=0.5)

    # Run
    sim = RunSim(init_v=-70, warmup=WARMUP, dt=DT)
    total_time = REPS * interval + COOL_DOWN
    sim.run(total_time)

    # Plot
    rec.plot(position="merge")