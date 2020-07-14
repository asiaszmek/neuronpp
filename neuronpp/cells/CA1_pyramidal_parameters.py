Rm_soma = 14000 # Somatic specific membrane resistance
Rm_non_trunk = Rm_soma  # Apical oblique specific membrane resistance
Rm_basal = Rm_soma  # Basal specific membrane resistance
Rm_axon = 1/0.00018 # Axonal specific membrane resistance
Rm_trunk = Rm_soma  # 36900  # Non-oblique dendritic specific Rm
axon_e_pas = -76.5
Ra_basal = 108  # Basal specific axial resistance
Ra_trunk = 108  # Somatic specific axial resistance
Ra_non_trunk = 108  # Somatic specific axial resistance
Ra_soma = 108  # Somatic specific axial resistance
Ra_axon = 74  # Axonal specific axial resistance

Cm_default = 1      # Default specific capacitance
Cm_axon = Cm_default  # Axonal specific capacitance
Cm_soma = 1.5     # Somatic specific capacitance
Cm_trunk = 1.5     # Trunk specific capacitance
Cm_non_trunk = 1.5     # Oblique specific capacitance
Cm_basal = 1.144    # Basal specific capacitance

soma_caL = .00006
soma_car = 0.00003
gsomacar = 0.00008
soma_calH = 0.0001  # 0.00017
soma_caT = 0.00003
soma_km = 0*0.001
potNa = 50
BK_channel_init = 0.9*1.5*0.03  # 0.03 flag
soma_nap_gnabar = 0*.5*0.000014
soma_K_nap = 4.5
soma_vhalf_nap = -60.4
soma_SK_channel = 0.7*4.5*0.0001  # 0.003 flag
soma_kap = 7*0.0005
soma_hbar = 1.8e-6
soma_kad = 7*0.0005
gna = 1.2*0.035  # 0.035
gkdr = 2.2*0.015  # flag
gkd = 0.005
gnadend = 0.015
gkdrdend = 2.2*0.015  # flag
gnanotrunk = 0*0.01  # 0.035
gkdrnotrunk = 2.2*0.015
AXKdr = 1
AXNa = 1
e_pas = -70
potK = -80
soma_K_h = 8.8
soma_vhalf_h = -82
cac_SK_channel = 0.00075
gbar_SK_channel = 0.5*soma_SK_channel
gkbar_BK_channel = 5.5*BK_channel_init


caT_distal_maxfactor = 4
caT_distal_distance = 350

SK_channel_distal_maxfactor = 1   # ORIG>> maximum cond. factor in dendrites
SK_channel_distal_distance = 200  # ORIG>> distance in dendrites for maximum cond.

dend_kap = 0.0025036

taur_cad = 20
potCa = 140
