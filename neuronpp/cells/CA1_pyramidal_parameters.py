Rm_soma = 14000 # Somatic specific membrane resistance
Rm_non_trunk = Rm_soma  # Apical oblique specific membrane resistance
Rm_basal = Rm_soma  # Basal specific membrane resistance
Rm_axon = 1/0.00018 # Axonal specific membrane resistance
Rm_trunk = Rm_soma  # 36900  # Non-oblique dendritic specific Rm
axon_e_pas = -76.5
Ra_basal = 108.4  # Basal specific axial resistance
Ra_trunk = 108.4  # Somatic specific axial resistance
Ra_non_trunk = 108.4  # Somatic specific axial resistance
Ra_soma = 108.4  # Somatic specific axial resistance
Ra_axon = 74  # Axonal specific axial resistance

Cm_default = 1      # Default specific capacitance
Cm_axon = Cm_default  # Axonal specific capacitance
Cm_soma = 1     # Somatic specific capacitance
Cm_trunk = 1     # Trunk specific capacitance
Cm_non_trunk = 1     # Oblique specific capacitance
Cm_basal = 1    # Basal specific capacitance

soma_caL = .00006
gsomacar = 0.00003
soma_calH = 0.00001  # 0.00017
soma_caT = 0.000003

soma_km = 0.002
axon_km = 0.024
soma_kap = 0.091
axon_kap = 0.16

potNa = 50
gna = 0.035  # 0.035

nax = 0.03
axon_nax = 0.06

gnadend = 0.015

AXNa = 1

BK_channel_init = 0.9*1.5*0.03  # 0.03 flag
gkbar_BK_channel = 5.5*BK_channel_init


soma_SK_channel = 0.7*4.5*0.0001  # 0.003 flag
cac_SK_channel = 0.00075
gbar_SK_channel = 0.5*soma_SK_channel
SK_channel_distal_maxfactor = 1   # ORIG>> maximum cond. factor in dendrites
SK_channel_distal_distance = 200  # ORIG>> distance in dendrites for maximum cond.


axon_kdr = 0.04
apic_kdr = 0.0073
dend_kdr = 0.0073
soma_kdr = 0.0073





potK = -90



caT_distal_maxfactor = 4
caT_distal_distance = 350


taur_cad = 20
potCa = 140
