#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

# Test Membrane Transport Mechanisms
from biocrnpyler import (
    Complex,
    Component,
    Diffusion_Facilitated_Carrier,
    Integration_MembraneProtein,
    ParameterKey,
    Transport_PrimaryActive_ABCexporter,
    Transport_SecondaryActive_Symporter,
    Transport_SecondaryActive_Antiporter,
    Diffusion_Simple,
    Diffusion_Facilitated_Channel,
    Species,
)


def contains(element, nested_array):
    """Recursively checks if an element is in a nested list."""
    return any(
        contains(element, sublist)
        if isinstance(sublist, list)
        else element == sublist
        for sublist in nested_array
    )


def total_length(nested_array):
    """Recursively counts the total number of elements in a nested list."""
    count = 0
    for item in nested_array:
        if isinstance(item, list):
            count += total_length(item)  # Recursively count sublist elements
        else:
            count += 1  # Count individual elements
    return count


class test_diffusion_simple:
    sd = Diffusion_Simple()
    substrate = Species('DMi')
    product = Species('DMo')
    c_fake = Species('C')

    # Test Update Species
    assert len(sd.update_species(substrate, product)) == 2
    assert substrate in sd.update_species(substrate, product)
    assert product in sd.update_species(substrate, product)

    # Test Update Reactions
    assert (
        len(
            sd.update_reactions(
                substrate,
                product,
                k_diff=1.0,
            )
        )
        == 1
    )
    assert (
        len(
            sd.update_reactions(
                substrate,
                product,
                k_diff=1.0,
                complex_species=c_fake,
            )
        )
        == 1
    )


class test_membrane_integration:
    mpi = Integration_MembraneProtein()
    MP = Species('MP1')
    MP.size = 2
    IMP = Species('IMP1')
    c1 = Complex([MP] * MP.size)
    c_fake = Species('C')

    # Test Update Species
    assert len(mpi.update_species(MP, IMP)) == 3
    assert c1 in mpi.update_species(MP, IMP)
    assert c_fake in mpi.update_species(MP, IMP, complex=c_fake)

    # Test Update Reactions
    # Define sensor parameter dictionary and component
    insertion_param_dict = {
        ParameterKey(
            mechanism='integration_membraneprotein',
            part_id=None,
            name='kb_oligomer',
        ): 2e-3,
        ParameterKey(
            mechanism='integration_membraneprotein',
            part_id=None,
            name='ku_oligomer',
        ): 2e-10,
        ParameterKey(
            mechanism='integration_membraneprotein', part_id=None, name='kex'
        ): 2e-3,
        ParameterKey(
            mechanism='integration_membraneprotein',
            part_id=None,
            name='kcat',
        ): 2e-10,
    }
    insertion_params = Component(
        'insertion_params', parameters=insertion_param_dict
    )

    assert len(mpi.update_reactions(MP, IMP, component=insertion_params)) == 2
    assert (
        len(
            mpi.update_reactions(
                MP,
                IMP,
                component=insertion_params,
                complex_species=c_fake,
            )
        )
        == 2
    )


class test_diffusion_facilitated_channel:
    st = Diffusion_Facilitated_Channel()
    MC = Species('MC1')
    substrate = Species('S1')
    product = Species('P1')
    c_fake = Species('C')

    # Test Update Species
    assert len(st.update_species(MC, substrate, product)) == 3

    # Test Update Reactions
    assert len(st.update_reactions(MC, substrate, product, k_diff=1.0)) == 1
    assert (
        len(
            st.update_reactions(
                MC,
                substrate,
                product,
                k_diff=1.0,
                complex_species=c_fake,
            )
        )
        == 1
    )


class test_diffusion_facilitated_carrier:
    ft = Diffusion_Facilitated_Carrier()
    MC = Species('MC1')
    substrate = Species('S1')
    product = Species('P1')
    carrier_in = Species(MC.name, material_type='protein', compartment=MC.compartment, attributes=['in'])
    c1 = Complex([product, MC])
    c2 = Complex([substrate, carrier_in])
    c_fake = Species('C')

    # Test Update Species
    assert total_length(ft.update_species(MC, substrate, product)) == 5

    # Test Update Reactions
    # Define sensor parameter dictionary and component
    transport_param_dict = {
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='kb_subMC',
        ): 2e-3,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='ku_subMC',
        ): 2e-10,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='kf_trnspMC',
        ): 2e-3,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='kr_trnspMC',
        ): 2e-3,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='kb_prodMC',
        ): 2e-3,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='ku_prodMC',
        ): 2e-10,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='k_out',
        ): 2e-3,
        ParameterKey(
            mechanism='diffusion_facilitated_carrier',
            part_id=None,
            name='k_in',
        ): 2e-10,
    }
    transport_params = Component(
        'transport_params', parameters=transport_param_dict
    )

    # Test Update Reactions
    assert (
        len(
            ft.update_reactions(
                MC, substrate, product, component=transport_params
            )
        )
        == 6
    )
    assert (
        len(
            ft.update_reactions(
                MC,
                substrate,
                product,
                component=transport_params,
                complex_dict=None,
            )
        )
        == 6
    )


class test_transport_primaryactive_abcexporter:
    pat = Transport_PrimaryActive_ABCexporter()
    MP = Species('MC1')
    MP.ATP = 2
    MP.attributes = ['exporter']
    substrate = Species('S1')
    product = Species('P1')
    energy = Species('E1')
    waste = Species('W1')
    c1 = Complex([substrate, MP])
    c2 = Complex([MP.ATP * [energy], c1])
    c3 = Complex([MP.ATP * [energy], product, MP])
    c4 = Complex([MP.ATP * [waste], MP])
    c_fake = Species('C')

    # Test Update Species
    assert (
        total_length(
            pat.update_species(MP, substrate, product, energy, waste)
        )
        == 9
    )

    # Test Update Reactions
    # Define sensor parameter dictionary and component
    transport_param_dict = {
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kb_subMP',
        ): 2e-3,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='ku_subMP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kb_subMPnATP',
        ): 2e-3,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='ku_subMPnATP',
        ): 2e-1,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kf_trnspMP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kr_trnspMP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kb_prodMP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='ku_prodMP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='kb_MP',
        ): 2e-10,
        ParameterKey(
            mechanism='transport_primaryactive_abcexporter',
            part_id=None,
            name='ku_MP',
        ): 2e-10,
    }
    transport_params = Component(
        'transport_params', parameters=transport_param_dict
    )

    assert (
        len(
            pat.update_reactions(
                MP,
                substrate,
                product,
                energy,
                waste,
                component=transport_params,
            )
        )
        == 8
    )


class test_transport_secondaryactive_symporter:
    st = Transport_SecondaryActive_Symporter(driving_ion={'Ion1': '1:1'})
    MC = Species('MC1')
    substrate = Species('S1')
    product = Species('P1')

    # Test Update Species
    assert total_length(st.update_species(MC, substrate, product)) == 9

    # Test Update Reactions
    symporter_param_dict = {
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kb_ionMC_out'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='ku_ionMC_out'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kb_subMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='ku_subMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kf_trnspMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kr_trnspMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kb_prodMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='ku_prodMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='kb_ionMC_in'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='ku_ionMC_in'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='k_out'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_symporter', part_id=None, name='k_in'): 2e-10,
    }
    transport_params = Component('transport_params', parameters=symporter_param_dict)

    assert len(st.update_reactions(MC, substrate, product, component=transport_params)) == 6


class test_transport_secondaryactive_antiporter:
    at = Transport_SecondaryActive_Antiporter(driving_ion={'Ion1': '1:1'})
    MC = Species('MC1')
    substrate = Species('S1')
    product = Species('P1')

    # Test Update Species
    assert total_length(at.update_species(MC, substrate, product)) == 9

    # Test Update Reactions
    antiporter_param_dict = {
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kb_ionMC_out'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='ku_ionMC_out'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kf_ionX'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kr_ionX'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kb_subMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='ku_subMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kf_trnspMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kr_trnspMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kb_prodMC'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='ku_prodMC'): 2e-10,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='kb_ionMC_in'): 2e-3,
        ParameterKey(mechanism='transport_secondaryactive_antiporter', part_id=None, name='ku_ionMC_in'): 2e-10,
    }
    transport_params = Component('transport_params', parameters=antiporter_param_dict)

    assert len(at.update_reactions(MC, substrate, product, component=transport_params)) == 6