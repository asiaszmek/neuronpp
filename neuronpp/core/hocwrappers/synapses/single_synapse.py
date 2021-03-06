from typing import Union, List, Optional, Iterable

from neuronpp.core.hocwrappers.netcon import NetCon
from neuronpp.core.decorators import non_removable_fields
from neuronpp.core.hocwrappers.hoc_wrapper import HocWrapper
from neuronpp.core.hocwrappers.synapses.synapse import Synapse
from neuronpp.core.hocwrappers.point_process import PointProcess


@non_removable_fields("sources", "target")
class SingleSynapse(HocWrapper, Synapse):
    def __init__(self, source, point_process: PointProcess, name,
                 netcon: Optional[Union[List[NetCon], NetCon]] = None, tag=None):
        self.point_process = point_process
        self.point_process_name = point_process.mod_name

        name = "%s[%s]" % (self.point_process_name, name)
        HocWrapper.__init__(self, hoc_obj=point_process.hoc, parent=point_process.parent, name=name)

        self.tag = tag

        self.sources = []
        if source is not None:
            self.sources.append(source)
        self.target = point_process.parent

        if netcon is None:
            self.netcons = []
        else:
            if not isinstance(netcon, Iterable):
                netcon = [netcon]
            self.netcons = netcon

    def make_event(self, time, use_global_sim_time=True):
        """
        Currently it makes events to all NetConns connected to the synapse.

        :param time:
            time in ms of the simulational event
        :param use_global_sim_time:
            If true it will use global time of hoc simulation
            (don't need to add h.t or sim.time the the event time)
        """
        for nc in self.netcons:
            nc.make_event(time, use_global_sim_time)

    def add_netcon(self, source: Optional[HocWrapper], weight: float = 1.0, delay: float = 1.0,
                   threshold: float = 10) -> NetCon:
        """
        Currently it allows to add single new source
        :param source:
            Can be only: hocwrappers.NetStim, hocwrappers.VecStim, hocwrappers.Sec or None.
            If it is Sec also loc param need to be defined.
            If remain None it will create NetConn with no source, which can be use as external
            event source
        :param weight:
        :param delay:
            in ms
        :param threshold:
            threshold for NetConn, default=10
        :return:
            created NetCon
        """
        nc = self.point_process.cell.add_netcon(source=source, point_process=self.point_process,
                                                netcon_weight=weight,
                                                delay=delay, threshold=threshold)
        self.netcons.append(nc)
        if source is not None:
            self.sources.append(source)

        return nc

    def __repr__(self):
        ncs = '_'.join([str(nc) for nc in self.netcons])
        return "{}[{}][{}][{}]".format(self.parent, self.__class__.__name__, self.name, ncs)
