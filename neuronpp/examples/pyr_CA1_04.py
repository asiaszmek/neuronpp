import os

from neuronpp.utils.simulation import Simulation
from neuronpp.utils.record import Record
from neuronpp.cells.migliore2018 import Migliore2018_CA1

path = os.path.dirname(os.path.abspath(__file__))
morpho_path = os.path.join(path, "..", "commons", "morphologies",
                           "asc", "mpg141208_B_idA.asc")

if __name__ == '__main__':
    # define cell
    cell = Migliore2018_CA1(name="CA1_PC_cAC_sig5")
    cell.load_morpho(filepath=morpho_path)
    for sec in cell.secs:
        print(sec.hoc.psection()["morphology"])
