import os

from neuronpp.cells.cell import Cell
from neuronpp.core.cells.core_hoc_cell import CoreHocCell

path = os.path.dirname(os.path.abspath(__file__))
f_path = os.path.join(path, "..", "commons/hocmodels/migliore2018")

class CA1_04(Cell, CoreHocCell):
    def __init__(self, name=None, model_folder=f_path, spine_number=0,
                 spine_secs_names="apic",
                 spine_seed: int = None):
        """
        :param name:
            The name of the cell
        :param model_folder:
            The folder where the main folder of Combe et al. 2018 model is located
        :param spine_number:
            The number of spines added to the model with random_uniform distribution to the sections
             specified by 'spine_sec' param.
        :param spine_secs_names:
            The section or sections where to put spines. It can be:
              * a string - as a filter name, so you can set "apic" to add spies to all apical
                dendrites

              * a regex, which need to be prefixed with 'regex:' string before
                eg. 'regex:(apic)|(basal)'
              will return all sections wich have a name containing 'apic' or 'basal' string

              * a list of existing sections in the cell
        :param spine_seed:
            Seed value for the random_uniform spike distribution. Default is None
            meaning - there is no seed
        """
        Cell.__init__(self, name, model_folder)
        CoreHocCell.__init__(self, name)

        main_file = os.path.join(model_folder, "load_cell_pyr_04.hoc")
        self.load_hoc(main_file, hoc_template_name="CA1_PC_cAC_sig5")

        
        # Add spines with AMPA and NMDA synapses


class CA1_08(Cell, CoreHocCell):
    def __init__(self, name=None, model_folder=f_path, spine_number=0,
                 spine_secs_names="apic",
                 spine_seed: int = None):
        """
        :param name:
            The name of the cell
        :param model_folder:
            The folder where the main folder of Combe et al. 2018 model is located
        :param spine_number:
            The number of spines added to the model with random_uniform distribution to the sections
             specified by 'spine_sec' param.
        :param spine_secs_names:
            The section or sections where to put spines. It can be:
              * a string - as a filter name, so you can set "apic" to add spies to all apical
                dendrites

              * a regex, which need to be prefixed with 'regex:' string before
                eg. 'regex:(apic)|(basal)'
              will return all sections wich have a name containing 'apic' or 'basal' string

              * a list of existing sections in the cell
        :param spine_seed:
            Seed value for the random_uniform spike distribution. Default is None
            meaning - there is no seed
        """
        Cell.__init__(self, name, model_folder)
        CoreHocCell.__init__(self, name)

        main_file = os.path.join(model_folder, "load_cell_pyr_08.hoc")
        self.load_hoc(main_file, hoc_template_name="CA1_PC_cAC_sig6")

        
        # Add spines with AMPA and NMDA synapses
