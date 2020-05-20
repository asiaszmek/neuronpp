
Rm_trunk     = 36900 # Non-oblique dendritic specific membrane resistance.
Rm_non_trunk = 36900 # Apical oblique specific membrane resistance
Rm_basal     = 15484.67 # Basal specific membrane resistance
Rm_tip       = 36900 # Tip specific membrane resistance
Rm_soma      = 20000 # Somatic specific membrane resistance
Rm_axon      = 28000 # Axonal specific membrane resistance

Ra_basal     = 150 # Basal specific axial resistance
Ra_trunk     = 150 # Somatic specific axial resistance
Ra_non_trunk = 150 # Somatic specific axial resistance
Ra_soma      = 150 # Somatic specific axial resistance
Ra_tip       = 150 # Apical tip specific axial resistance
Ra_axon      = 50  # Axonal specific axial resistance

Cm_default   = 1          # Default specific capacitance
Cm_axon      = Cm_default # Axonal specific capacitance
Cm_soma      = 1          # Somatic specific capacitance
Cm_soma      = 1.5        # Somatic specific capacitance
Cm_trunk     = 1.192      # Trunk specific capacitance
Cm_trunk     = 1.5        # Trunk specific capacitance
Cm_non_trunk = 1.192      # Oblique specific capacitance
Cm_non_trunk = 1.5        # Oblique specific capacitance
Cm_basal     = 1.144      # Basal specific capacitance
Cm_tip       = 1.192      # Apical tip specific capacitance
soma_caL = 0.00006/10
soma_car = 0.00003
gsomacar = 0.00008
soma_caLH = 0.0001 #  0.00017
soma_caT = 0.0003
soma_km = 0*0.001
potNa = 50 
mykca_init = 0.9*1.5*0.03 # 0.03 flag
soma_nap_gnabar = 0*.5*0.000014
K_nap = 4.5
vhalf_nap = -60.4
soma_kca = 0.7*4.5*0.0001# 0.003 flag
soma_kap = 7*0.0005
soma_hbar = 1.8e-6
soma_kad = 7*0.0005
gna = 1.2*0.035 #  0.035
gkdr = 2.2*0.015 # flag
gkd = 0.005
gnadend = 0.015
gkdrdend = 2.2*0.015 # flag
gnanotrunk = 0*0.01 # 0.035
gkdrnotrunk = 2.2*0.015
AXKdr = 1
AXNa = 1
e_pas = -70
potK = -80
soma_K_h = 8.8
soma_vhalf_h = -82
cac_kca=0.00075 //0.0005
gbar_kca = 0.5*soma_kca
gkbar_mykca = 5.5*mykca_init