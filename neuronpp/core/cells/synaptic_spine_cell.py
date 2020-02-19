from neuronpp.core.cells.spine_cell import SpineCell
from neuronpp.core.cells.complex_synaptic_cell import ComplexSynapticCell


class SynapticSpineCell(SpineCell, ComplexSynapticCell):
    def __init__(self, name=None, compile_paths=None):
        ComplexSynapticCell.__init__(self, name, compile_paths=compile_paths)
        SpineCell.__init__(self, name)

    def add_synapses_with_spine(self, source, mod_name: str, secs, weight=1, rand_weight=False,
                                number=1, delay=0, head_nseg=2, neck_nseg=2, tag: str = None, **synaptic_params):
        """

        :param source:
            Can be only: hocwrappers.NetStim, hocwrappers.VecStim, hocwrappers.Sec or None. If it is Sec also loc param need to be defined.
            If None it will create NetConn with no source, which can be use as external event source
        :param weight:
        :param rand_weight:
            if True, will find rand weight [0,1) and multiply this by weight.
        :param number:
        :param tag:
        :param mod_name:
        :param delay:
        :param secs:
        :param head_nseg:
        :param neck_nseg:
        :param synaptic_params:
        :return:
        """
        heads, _ = self.make_spines(spine_number=number, secs=secs, head_nseg=head_nseg, neck_nseg=neck_nseg)

        # loc=1.0 put synase on the top of the spine's head
        syns = []
        for h in heads:
            h_segment = h(1.0)
            syn = self.add_sypanse(source=source, seg=h_segment, mod_name=mod_name, weight=weight,
                                   rand_weight=rand_weight, delay=delay, tag=tag, **synaptic_params)
            syns.append(syn)

        return syns, heads
