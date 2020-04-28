from neuron.rxd import rxd

from neuronpp.core.cells.rxd_tools import RxDTool
from neuronpp.core.cells.spine_cell import SpineCell
from neuronpp.core.hocwrappers.rxd import RxD
from neuron.units import nM, ms, sec


class CaCell(SpineCell):
    def _define_regions(self, secs, dx_3d_size):
        """
        Set up cytosol, membrane and the extracellular compartment in secs. 

        :param secs:
            List of sections (hoc objects)
        :param dx_3d_size:
            Size of a single voxel
        """
        cyt = rxd.Region(secs=secs, name="cytosol", nrn_region='i',
                              dx=dx_3d_size)
        mem = rxd.Region(secs=secs, name="membrane",
                              geometry=rxd.membrane())
        extc = rxd.Region(secs=secs, name='extracellular', nrn_region='o',
                               geometry=rxd.Shell(1, 2))
   
        return cyt, mem, extc 


    def __init__(self, name=None, compile_paths=None, regions=None,
                 dx_3d_size=None):
        """
        :param name:
            Name of the cell
        """
        SpineCell.__init__(self, name, compile_paths=compile_paths)
        if regions is None:
            secs = [sec.hoc for sec in self.secs]
        elif isinstance(regions, str):
            regions = self.filter_secs(regions)
        self.cyt, self.mem, self.extc = self._define_regions(self, secs,
                                                             dx_3d_size)
        self.ca = rxd.Species(regions=self.cyt, initial=75 * nM,
                              name='Ca', charge=2, d=0.016)
         
    ### Add pumps and buffers
