
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex

#Test Membrane Transport Mechanisms
from biocrnpyler import Membrane_Signaling_Pathway_MM

class test_membrane_signaling_MM():
    tcs = Membrane_Signaling_Pathway_MM()
    MSP = Species("MSP1")
    MSP.ATP=2
    RP = Species("RP1")
    sub_assign = Species("S1")
    sub_signal = Species("S2")
    product= Species('RP_active')
    energy = Species("E1")
    waste = Species("W1")
    c1 =  Complex([sub_signal, MSP])
    c2 =  Complex([MSP.ATP*[energy], c1])
    c3 =  Complex([c1, MSP.ATP*[waste], sub_assign])
    c4 =  Complex([c1, sub_assign])
    c5 = Complex([c4, RP])
    c6 = Complex([c1, RP, sub_assign])

    c_fake = Species("C")
    
    #Test Update Species
    assert len(tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)) == 7
    assert c1 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)
    # assert c2 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)
    # assert c3 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)
    # assert c4 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)
    # assert c5 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)
    # assert c6 in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)

    assert c_fake in tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste, complex = c_fake)
    
    #Test Update Reactions
    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, product, energy, waste, kb_sigMS = 2e-3, ku_sigMS = 2e-10, 
        kb_autoPhos = 2e-3, ku_autoPhos = 2e-10, ku_waste = 1e-1, kb_phosRP = 2e-3, ku_phosRP = 2e-10, 
        k_phosph = 1e-1, ku_activeRP = 2e-1,ku_dephos = 2e-10)) == 8

    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, product, energy, waste, kb_sigMS = 2e-3, ku_sigMS = 2e-10,
        kb_autoPhos = 2e-3, ku_autoPhos = 2e-10, ku_waste = 1e-1, kb_phosRP = 2e-3, ku_phosRP = 2e-10, 
        k_phosph = 1e-1, ku_activeRP = 2e-1, ku_dephos = 2e-10, complex_species = c_fake,)) == 8