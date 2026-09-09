#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

# Test Membrane Signaling Mechanisms
from biocrnpyler import (
    Complex,
    Component,
    ParameterKey,
    Species,
    # Sensor_TwoComponentSignaling,
)
# Assuming the signaling mechanisms are imported from their respective module
from biocrnpyler.mechanisms.signaling import Sensor_TwoComponentSignaling


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


class test_sensor_twocomponentsignaling:
    # Initialize mechanism
    tcs = Sensor_TwoComponentSignaling()

    # Define species
    MS = Species('EnvZ')
    MS.ATP = 2
    RP = Species('OmpR')
    Pi = Species('Phosphate')
    Sig = Species('Osmolarity')
    Prod = Species('OmpR_active')
    ATP = Species('ATP')
    ADP = Species('ADP')

    # Pre-define some expected complexes to test for
    c_activated_ms = Complex([Sig, MS])
    c_atp_activated_ms = Complex([MS.ATP * [ATP], c_activated_ms])

    # Test Update Species
    def test_update_species(self):
        # The method returns 6 base species + a list of 7 complexes = 13 total species
        generated_species = self.tcs.update_species(
            self.MS, self.RP, self.Pi, self.Sig, self.Prod, self.ATP, self.ADP
        )

        assert total_length(generated_species) == 13
        assert contains(self.MS, generated_species)
        assert contains(self.c_activated_ms, generated_species)
        assert contains(self.c_atp_activated_ms, generated_species)

    # Test Update Reactions
    def test_update_reactions(self):
        # Define parameter dictionary and component
        signaling_param_dict = {
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='kb_sigMS'): 1.0,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_sigMS'): 0.1,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='kb_autoPhos'): 1.0,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_autoPhos'): 0.1,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='k_hydro'): 0.5,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_waste'): 1.0,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='kb_phosRP'): 1.0,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_phosRP'): 0.1,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='k_phosph'): 0.5,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_activeRP'): 1.0,
            ParameterKey(mechanism='membrane_sensor', part_id=None, name='ku_dephos'): 0.01,
        }
        signaling_params = Component(
            'signaling_params', parameters=signaling_param_dict
        )

        generated_reactions = self.tcs.update_reactions(
            self.MS,
            self.RP,
            self.Pi,
            self.Sig,
            self.Prod,
            self.ATP,
            self.ADP,
            component=signaling_params
        )

        # The mechanism generates exactly 9 reactions for the signaling cascade
        assert (len(generated_reactions) == 9)
