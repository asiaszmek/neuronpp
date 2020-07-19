import os
from numpy import exp
from neuron import h
from neuronpp.cells.cell import Cell
from neuronpp.cells.morphology_points_short_axon import axon_points
from neuronpp.cells.morphology_points import trunk_points
from neuronpp.cells.morphology_points import points_apic, points_apic_continued
from neuronpp.cells.morphology_points import points_dend, points_dend_continued
import neuronpp.cells.CA1_pyramidal_parameters as params


path = os.path.dirname(os.path.abspath(__file__))
f_path = os.path.join(path, "..", "commons/mods/CA1_pyramidal")
maximum_segment_length = 75

def dist_e_pas(x):
    return -65.168 - 5*x/150

def dist_kad(x):
    return 15./(1. + exp((300-x)/50))*0.011237

def dist_hd(x):
    return (1. + 3./100. * x)*1.9199e-05


class CA1PyramidalCell(Cell):
    def distribute_channel(self, x, mech, mech_param, val):
        mech_obj = getattr(x, mech)
        setattr(mech_obj, mech_param, val)

    def add_pas(self):
        soma = self.filter_secs("soma", as_list=True)[0].hoc
        for section in self.secs:
            sec = section.hoc
            dist = h.distance(sec(0.5)) - h.distance(soma(0.5))
            if "soma" in sec.name():
                sec.insert("pas")
                sec.g_pas = 1/params.Rm_soma
                sec.e_pas = dist_e_pas(dist)
                sec.Ra = params.Ra_soma
                sec.cm = params.Cm_soma

            if "axon" in sec.name():
                sec.insert("pas")
                sec.g_pas = 1/params.Rm_axon
                sec.e_pas = params.axon_e_pas
                sec.Ra = params.Ra_axon
                sec.cm = params.Cm_axon

            if "trunk" in sec.name():
                sec.insert("pas")
                sec.g_pas = 1/params.Rm_trunk
                sec.e_pas = dist_e_pas(dist)
                sec.Ra = params.Ra_trunk
                sec.cm = params.Cm_trunk
            if "apic" in sec.name():
                sec.insert("pas")
                sec.g_pas = 1/params.Rm_non_trunk
                sec.e_pas = dist_e_pas(dist)
                sec.Ra = params.Ra_non_trunk
                sec.cm = params.Cm_non_trunk
            if "dend" in sec.name():
                sec.insert("pas")
                sec.g_pas = 1/params.Rm_basal
                sec.e_pas = dist_e_pas(dist)
                sec.Ra = params.Ra_basal
                sec.cm = params.Cm_basal

    def add_sodium_channels(self):
        for sec in self.secs:
            sec.hoc.insert("nax")
            sec.hoc.ena = params.potNa
            if "axon" in sec.hoc.name():
                sec.hoc.gbar_nax = params.axon_nax
            elif  "trunk" in sec.hoc.name():
                sec.hoc.insert("na3")
                sec.gbar_na3 = params.gnadend
                sec.hoc.gbar_nax = params.nax

            elif "apic" in sec.hoc.name() :
                sec.hoc.insert("na3")
                sec.gbar_na3 = params.gnadend
                sec.hoc.gbar_nax = params.nax
            elif "dend" in sec.hoc.name():
                sec.hoc.insert("na3dend")
                sec.hoc.insert("nap")
                sec.hoc.gbar_nax = params.nax
            elif "soma" in sec.hoc.name():
                sec.hoc.gbar_nax = params.nax
                sec.hoc.insert("na3")
                sec.gbar_na3 = params.gna
                sec.ena = params.potNa

    def add_potassium_channels(self):
        soma = self.filter_secs("soma", as_list=True)[0].hoc

        for sec in self.secs:
            sec.hoc.insert("kdr")
            xdist = h.distance(sec.hoc(0.5))
            if "axon" in sec.hoc.name():
                sec.hoc.insert("kmb")
                sec.hoc.insert("kap")
                sec.hoc.ek = params.potK
                sec.hoc.gbar_kmb = params.axon_km
                sec.hoc.gkabar_kap = params.axon_kap
                sec.hoc.gkdrbar_kdr = params.axon_kdr
            elif "apic" in sec.hoc.name() or "trunk" in sec.hoc.name():
                sec.hoc.insert("hd")
                sec.hoc.insert("kad")
                sec.hoc.ek = params.potK
                sec.hoc.gkdrbar_kdr = params.apic_kdr
            elif "dend" in sec.hoc.name():
                sec.hoc.insert("kad")
                sec.hoc.insert("hd")
                sec.hoc.ek = params.potK
                sec.hoc.gkdrbar_kdr = params.dend_kdr

            elif "soma" in sec.hoc.name():
                sec.hoc.insert("kmb")
                sec.hoc.insert("kap")
                sec.hoc.insert("hd")
                sec.hoc.ek = params.potK
                sec.hoc.gkabar_kap = params.soma_kap
                sec.hoc.gbar_kmb = params.soma_km
                sec.hoc.gkdrbar_kdr = params.soma_kdr

            if "kad" in sec.hoc.psection()["density_mechs"]:
                for seg in sec.hoc:
                    xdist = h.distance(seg) - h.distance(soma(0.5))
                    val = dist_kad(xdist)
                    self.distribute_channel(seg, "kad", "gkabar", val)

            if "hd" in sec.hoc.psection()["density_mechs"]:
                for seg in sec.hoc:
                    xdist = h.distance(seg) - h.distance(soma(0.5))
                    val = dist_hd(xdist)
                    self.distribute_channel(seg, "hd", "ghdbar", val)

    def _shorten_axon(self):
        L_target = 60
        nseg0 = 5
        nseg_total = 2*nseg0
        nSec = len(self.axon)
        if not nSec:
            raise "No axon"
        i1 = self.axon[0].v
        try:
            i2 = self.axon[1].v
        except IndexError:
            i2 = i1
        count = 0
        diams = self.axon[0].diam
        # delete old axon
        axon = self.remove_secs("axon")
        del self.axon
        self.axon = []
        
        axon_0 = self.add_sec("axon_0")
        axon_1 = self.add_sec("axon_1")
        soma = self.filter_secs("soma")
        self.connect_secs(axon_0, soma)
        self.connect_secs(axon_1, axon_0)
        self.axon.append(axon_0.hoc)
        self.axon.append(axon_1.hoc)
        for i, axon in enumerate(self.axon):
            axon.L = L_target/2
            axon.nseg = int(nseg_total/2)
            axon.diam = diams
            if i == 0:
                axon.v = i1
            else:
                axon.v = i2

    def make_axon(self):
        # axon
        self.axon.hoc.pt3dclear()
        for points in axon_points:
            h.pt3dadd(*points, sec=self.axon.hoc)

    def make_soma(self):
        self.soma.hoc.pt3dclear()
        h.pt3dadd(10, 0, 30, 20, sec=self.soma.hoc)
        h.pt3dadd(21.5, 0.4, 30, 2.1, sec=self.soma.hoc)

    def make_trunk(self):
        for i, sec in enumerate(self.trunk):
            sec.hoc.pt3dclear()
            for points in trunk_points[i]:
                h.pt3dadd(*points, sec=sec.hoc)

    def make_apic(self):
        len_1 = len(points_apic)
        for i, sec in enumerate(self.apic):
            sec.hoc.pt3dclear()
            if i < len_1:
                for points in points_apic[i]:
                    h.pt3dadd(*points, sec=sec.hoc)
            else:
                for points in points_apic_continued[i - len_1]:
                    h.pt3dadd(*points, sec=sec.hoc)

    def make_dend(self):
        len_1 = len(points_dend)
        for i, sec in enumerate(self.dend):
            sec.hoc.pt3dclear()
            if i < len_1:
                for points in points_dend[i]:
                    h.pt3dadd(*points, sec=sec.hoc)
            else:
                for points in points_dend_continued[i - len_1]:
                    h.pt3dadd(*points, sec=sec.hoc)

    def make_morphology(self):
        self.soma = self.add_sec("soma")
        self.axon = self.add_sec("axon")
        self.trunk = []
        for i in range(19):
            self.trunk.append(self.add_sec("trunk_%d" % i))
        self.apic = []
        for i in range(72):
            self.apic.append(self.add_sec("apic_%i" % i))
        self.dend = []
        for i in range(51):
            self.dend.append(self.add_sec("dend_%i" % i))

        self.connect_secs( self.axon, self.soma)
        for i in range(19):
            if i == 0:
                self.connect_secs(self.trunk[i], self.soma)
            else:
                self.connect_secs(self.trunk[i], self.trunk[i - 1])

        for i in range(42, 48):
            self.connect_secs(self.apic[i], self.apic[i - 1])
        self.connect_secs(self.apic[48], self.apic[46])
        for i in range(49, 51):
            self.connect_secs(self.apic[i], self.apic[48])
        self.connect_secs(self.apic[51], self.apic[45])
        self.connect_secs(self.apic[52], self.apic[44])
        for i in range(53, 55):
            self.connect_secs(self.apic[i], self.apic[52])
        self.connect_secs(self.apic[55], self.trunk[17])
        self.connect_secs(self.apic[56], self.trunk[18])
        self.connect_secs(self.apic[57], self.apic[56])
        self.connect_secs(self.apic[58], self.apic[56])
        self.connect_secs(self.apic[59], self.apic[55])
        for i in [60, 61]:
            self.connect_secs(self.apic[i], self.apic[i - 1])
        self.connect_secs(self.apic[62], self.apic[60])
        self.connect_secs(self.apic[63], self.apic[59])
        for i in [64, 65]:
            self.connect_secs(self.apic[i], self.apic[i - 1])
        self.connect_secs(self.apic[66], self.apic[64])
        self.connect_secs(self.apic[67], self.apic[63])
        for i in [68, 69]:
            self.connect_secs(self.apic[i], self.apic[i - 1])
        self.connect_secs(self.apic[70], self.apic[68])
        self.connect_secs(self.apic[71], self.apic[67])
        self.connect_secs(self.apic[0], self.trunk[0])
        for i in [1, 2]:
            self.connect_secs(self.apic[i], self.apic[0])
        for i in [3, 4]:
            self.connect_secs(self.apic[i], self.apic[2])
        for i in [5, 6]:
            self.connect_secs(self.apic[i], self.trunk[i-4])
        for i in [7, 8]:
            self.connect_secs(self.apic[i], self.apic[6])
        for i in [9, 10]:
            self.connect_secs(self.apic[i], self.apic[8])
        for i in [11, 12]:
            self.connect_secs(self.apic[i], self.apic[10])
        self.connect_secs(self.apic[13], self.trunk[3])
        for i in [14, 15]:
            self.connect_secs(self.apic[i], self.apic[13])
        for i in [16, 17]:
            self.connect_secs(self.apic[i], self.trunk[i - 12])
        for i in [18, 19]:
            self.connect_secs(self.apic[i], self.apic[17])
        self.connect_secs(self.apic[20], self.trunk[6])
        for i in [21, 22]:
            self.connect_secs(self.apic[i], self.apic[20])
        for i in [23, 24]:
            self.connect_secs(self.apic[i], self.apic[22])
        for i in [25, 26, 27]:
            self.connect_secs(self.apic[i], self.trunk[i - 18])
        for i in [28, 29]:
            self.connect_secs(self.apic[i], self.apic[27])
        for i in[30, 31, 32]:
            self.connect_secs(self.apic[i], self.trunk[i-20])
        self.connect_secs(self.apic[33], self.apic[42])
        self.connect_secs(self.apic[34], self.apic[43])
        self.connect_secs(self.apic[35], self.trunk[14])
        self.connect_secs(self.apic[36], self.trunk[15])
        self.connect_secs(self.apic[37], self.trunk[16])
        for i in [38, 39]:
            self.connect_secs(self.apic[i], self.apic[i - 1])
        self.connect_secs(self.apic[40], self.apic[38])
        self.connect_secs(self.apic[41], self.apic[37])
        self.connect_secs(self.dend[0], self.soma, 0, 0)
        for i in range(1, 4):
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[4], self.dend[2])
        self.connect_secs(self.dend[5], self.dend[1])
        for i in [6, 7]:
            self.connect_secs(self.dend[i], self.dend[5])
        for i in [8, 9]:
            self.connect_secs(self.dend[i], self.dend[7])
        for i in [10, 11]:
            self.connect_secs(self.dend[i], self.dend[9])
        self.connect_secs(self.dend[12], self.dend[0])
        for i in range(13, 16):
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[16], self.dend[14])
        self.connect_secs(self.dend[17], self.dend[13])
        for i in [18, 19]:
            self.connect_secs(self.dend[i], self.dend[17])
        self.connect_secs(self.dend[20], self.dend[12])
        for i in [21, 22]:
            self.connect_secs(self.dend[i], self.dend[20])
        self.connect_secs(self.dend[23], self.soma, 0, 0)
        for i in range(24, 27):
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[27], self.dend[25])
        self.connect_secs(self.dend[28], self.dend[24])
        self.connect_secs(self.dend[29], self.dend[23])
        for i in [30, 31]:
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[32], self.dend[30])
        for i in [33, 34]:
            self.connect_secs(self.dend[i], self.dend[32])
        self.connect_secs(self.dend[35], self.dend[29])
        for i in range(36, 39):
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[39], self.dend[37])
        self.connect_secs(self.dend[40], self.dend[36])
        self.connect_secs(self.dend[41], self.dend[35])
        self.connect_secs(self.dend[42], self.soma, child_loc=0, parent_loc=0)
        for i in [43, 44]:
            self.connect_secs(self.dend[i], self.dend[i - 1])
        self.connect_secs(self.dend[45], self.dend[43])
        for i in [46, 47]:
            self.connect_secs(self.dend[i], self.dend[45])
        self.connect_secs(self.dend[48], self.dend[42])
        for i in [49, 50]:
            self.connect_secs(self.dend[i], self.dend[48])

        self.make_axon()
        self.make_soma()
        self.make_trunk()
        self.make_apic()
        self.make_dend()

    def add_soma_mechanisms(self):
        soma = self.filter_secs(name="soma", as_list=True)
        for s in soma:
            sec = s.hoc

            sec.insert("cal")
            sec.gcalbar_cal = params.soma_caL/10
            
            sec.insert("cat")
            sec.gcatbar_cat = params.soma_caT
        
            sec.insert("car")
            sec.gcabar_car = params.gsomacar
            
            sec.insert("SK_channel")
            sec.cac_SK_channel = params.cac_SK_channel
            sec.gbar_SK_channel = params.gbar_SK_channel
            
            sec.insert("BK_channel")  # K(Ca) fAHP potassium type current
            sec.gkbar_BK_channel = params.gkbar_BK_channel


    def add_trunk_mechanisms(self):
        trunk = self.filter_secs("trunk")
        for s in trunk:
            sec = s.hoc
            sec.insert("car")
            sec.gcabar_car = 0.1*params.soma_car
            sec.insert("cat")
            sec.insert("SK_channel")
            sec.cac_SK_channel = params.cac_SK_channel
            sec.insert("calH")
            sec.insert("BK_channel")
            for i, seg in enumerate(sec):
                if i == sec.nseg - 1:
                    xdist = h.distance(sec(1.0))
                else:
                    xdist = h.distance(seg)

                fr = xdist/params.caT_distal_distance

                if xdist > 50:
                    self.distribute_channel(seg, "calH", "gcalbar",
                                             2*params.soma_calH)
                else:
                    self.distribute_channel(seg, "calH", "gcalbar",
                                             0.1*params.soma_calH)

                if xdist < 100:
                    self.distribute_channel(seg, "cat", "gcatbar", 0)
                else:
                    val = params.caT_distal_maxfactor*params.soma_caT*fr
                    self.distribute_channel(seg, "cat", "gcatbar", val)

                if xdist < params.SK_channel_distal_distance and xdist > 50:
                    self.distribute_channel(seg, "SK_channel", "gbar",
                                             5*params.soma_SK_channel)
                    self.distribute_channel(seg, "BK_channel", "gkbar",
                                             2*params.BK_channel_init)

                else:
                    self.distribute_channel(seg, "SK_channel", "gbar",
                                             0.5*params.soma_SK_channel)
                    self.distribute_channel(seg, "BK_channel", "gkbar",
                                             0.5*params.BK_channel_init)

    def add_apical_mechanisms(self):
        apic = self.filter_secs("apic", as_list=True)
        for s in apic:
            sec = s.hoc
            sec.insert("car")
            sec.insert("calH")
            sec.gcabar_car = 0.1*params.soma_car
            sec.insert("cat")
            sec.insert("SK_channel")
            sec.cac_SK_channel = params.cac_SK_channel

            sec.insert("BK_channel")
            for i, seg in enumerate(sec):
                if i == sec.nseg - 1:
                    xdist = h.distance(sec(1.0))
                else:
                    xdist = h.distance(seg)
                fr = xdist/params.caT_distal_distance
                if xdist > 50:
                    new_calH = 2*params.soma_calH
                else:
                    new_calH = 0.1*params.soma_calH
                if xdist < 100:
                    new_caT = 0
                else:
                    new_caT = params.caT_distal_maxfactor*params.soma_caT*fr

                if xdist < params.SK_channel_distal_distance and xdist > 50:
                    new_SK = 5*params.soma_SK_channel
                    new_BK = 2*params.BK_channel_init
                else:
                    new_SK = 0.5*params.soma_SK_channel
                    new_BK = 0.5*params.BK_channel_init
                self.distribute_channel(seg, "calH", "gcalbar",
                                         new_calH)
                self.distribute_channel(seg, "cat", "gcatbar", new_caT)
                self.distribute_channel(seg, "SK_channel", "gbar",
                                         new_SK)
                self.distribute_channel(seg, "BK_channel", "gkbar",
                                         new_BK)

            
    def add_calcium(self, decay=True):
        if decay:
            ca_sections = self.filter_secs("soma", as_list=True) + self.filter_secs("trunk", as_list=True) + self.filter_secs("apic", as_list=True)
            for sec in ca_sections:
                section = sec.hoc
                section.insert("cad")
                section.taur_cad = params.taur_cad
                section.eca = 140#params.potCa
        else:
            print("Unimplemented mechanism")

    def __init__(self, name=None, compile_paths=f_path, morpho="migliore"):
        """
        :param name:
            The name of the cell
        :param compile_paths:
            Folder with channels
        """
        Cell.__init__(self, name=name, compile_paths=compile_paths)
        if morpho == "migliore":
            # adjust segment_number
            morpho_path = os.path.join(path, "..", "commons", "morphologies",
                                       "asc", "mpg141208_B_idA.asc")
            self.load_morpho(filepath=morpho_path)
            self._shorten_axon()
            for sec in self.secs:
                sec.hoc.nseg = 1 + int(sec.hoc.L/30)
        else:
            self.make_morphology()
            for sec in self.secs:
                sec.hoc.nseg = 1 + int(sec.hoc.L/maximum_segment_length)
        h.distance()
        self.add_pas()
        self.add_sodium_channels()
        self.add_potassium_channels()
        # self.add_soma_mechanisms()
        # self.add_axon_mechanisms()
        # self.add_trunk_mechanisms()
        # self.add_apical_mechanisms()
        # h.celsius = 34
        # self.add_calcium()
        # for sec in self.secs:
        #     if h.ismembrane("ca_ion", sec=sec.hoc):
        #         sec.hoc.eca = 140
        #         h.ion_style("ca_ion",0, 1, 0, 0, 0, sec=sec.hoc)
