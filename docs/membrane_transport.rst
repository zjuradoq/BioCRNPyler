==============================
Membrane-associated Classes and Mechanisms
==============================

The following Jupyter notebooks provide a set of tutorial
introductions to BioCRNpyler's membrane-associate features.

-------------
Introduction
-------------

BioCRNpyler supports modeling genetic circuits with membrane-associated 
features using computational tools. This functionality enables users to 
construct and simulate genetic circuits that incorporate membrane components 
and mechanisms within simplified cell-free system models. 
 

BioCRNpyler's membrane-associated features allow for the representation of 
inducer dynamics, including diffusion across membranes or transport via 
channels and transporters. This facilitates the design and refinement of models 
that connect mechanistic biology with computational analysis, bridging the gap 
between conceptual design and experimental implementation in synthetic biology.


.. image:: figures/membrane_transport.png
   :alt: Types of membrane transport components and mechanisms
   :align: center
   :width: 400px

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Membrane Components/Mechanisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following membrane-associated mechanisms that are available in BioCRNpyler:

- **:ref:`Simple Diffusion <simple-diffusion>`**: Models the passive movement of small, nonpolar molecules across the membrane, driven 
  by concentration gradients, without the need for membrane proteins or energy input.

- Membrane protein-mediated mechanisms
    - **:ref:`Membrane Protein Integration <Membrane-Protein-Integration>`**: Models the the insertion and proper orientation of proteins into the membrane, ensuring 
    their structural and functional integration within the lipid bilayer.

    - **:ref:`Simple Transport <simple-transport>`**: Models the passive movement of substrates through membrane pores/channels along concentration
   gradients, without requiring energy input.

    - **:ref:`Facilitated Transport <facilitated-transport>`**: Models the passive movement of substrates along concentration gradients by binding to carrier 
  proteins that undergo conformational changes, without requiring energy input.
  
    - **:ref:`Primary Active Transport <primary-active-transport>`**: Models the active movement of substrates against concentration gradients by binding to membrane 
  pumps, which undergo conformational changes driven by energy input (e.g., ATP).

    - **:ref:`Two-Component Signaling <two-component-signaling>`**: Models the environmental sensing through a signaling pathway involving a sensor kinase and 
  phosphorylation of a response regulator protein, enabling adaptive cellular responses.

.. - Multicellular communication
.. - Examples

The figure presented below illustrates the various options available for modeling transport 
and two-component signaling within BioCRNpyler. It specifically highlights the membrane 
components, indicated by purple boxes, along with their corresponding mechanisms, represented 
by the blue box.

.. image:: figures/membrane_model_flowchart.png
   :alt: Flowchart for membrane modeling options
   :align: center
   :width: 400px

.. _simple-diffusion:

----------------
Simple Diffusion
----------------

Simple diffusion allows molecules to passively cross membranes down their concentration 
gradient. This is the most basic mechanism by which molecules can traverse a membrane, commonly 
referred to as passive diffusion. In this process, a molecule can dissolve in the lipid bilayer, 
diffuse across it, and reach the other side. This mechanism does not require the assistance of 
membrane proteins, and the transport direction is determined by the concentration gradient, 
moving from areas of high concentration to areas of low concentration.

~~~~~~~~~~
Component: ``DiffusibleMolecule()``
~~~~~~~~~~

A Diffusible Molecule refers to a class of molecules that can pass through cell membranes 
without assistance. Examples of such molecules include gases like oxygen (O<sub>2</sub>) and 
carbon dioxide (CO<sub>2</sub>), as well as small polar but uncharged molecules. In contrast, 
larger uncharged molecules and charged molecules require membrane proteins for transport across 
the membrane.

The following code defines a diffusible molecule called ``S``:

.. code-block:: python

    S = DiffusibleMolecule('name')

Unless otherwise specified, the species ``S`` will reside in the ``internal`` compartment.  
The membrane component ``DiffusibleMolecule(Component)`` will then create a species ``product``,  
which is a copy of ``S`` but located in the ``external`` compartment.

To access more information about this component, use:

.. code-block:: python

    help(DiffusibleMolecule)

~~~~~~~~~~
Mechanism: ``Simple_Diffusion()``
~~~~~~~~~~

In BioCRNpyler, the ``DiffusibleMolecule`` component uses the mechanism ``Simple_Diffusion``, which can be defined as:

.. code-block:: python

    mech_tra = Simple_Diffusion()
    transport_mechanisms = {mech_tra.mechanism_type: mech_tra}

~~~~~~~~~~
Creating a Mixture
~~~~~~~~~~

We can now create a mixture that uses this mechanism:

.. code-block:: python

    M = Mixture(
        "DiffusibleMolecule",
        components=[S],
        parameter_file="membrane_toolbox_parameters.txt",
        mechanisms=transport_mechanisms
    )

~~~~~~~~~~
Compiling the CRN
~~~~~~~~~~

Finally, compile the chemical reaction network:

.. code-block:: python

    CRN = M.compile_crn()

~~~~~~~~~~
Example: 
~~~~~~~~~~

Consider the following diffusion step:

1. **Diffusion of small molecules (i.e. Nitrate):**

.. math::

    NO3_{internal} \rightleftharpoons NO3_{external}

Define diffusible molecule:
--------------------------

.. code-block:: python

    # Define diffusible molecule
    NO3 = DiffusibleMolecule('NO3')
    # nitrate=

    # Mechanisms
    mech_tra = Simple_Diffusion()
    transport_mechanisms = {mech_tra.mechanism_type: mech_tra}

    # Create mixture
    M0 = Mixture("Diffusible_Molecule", components=[NO3],
                 parameter_file="membrane_toolbox_parameters.txt",
                 mechanisms=transport_mechanisms)

    # Compile the CRN with Mixture.compile_crn
    CRN = M0.compile_crn()

    # Print the CRN to see what you created
    print(CRN.pretty_print())

Console Output:
--------------

.. code-block:: text

    Species(N = 2) = {
    NO3 (@ 0),  NO3 (@ 0),  
    }

    Reactions (1) = [
    0. NO3 <--> NO3
     Kf=k_forward * NO3_Internal
     Kr=k_reverse * NO3_External
      k_forward=0.0002
      found_key=(mech=simple_diffusion, partid=None, name=k_diff).
      search_key=(mech=simple_diffusion, partid=NO3, name=k_diff).
      k_reverse=0.0002
      found_key=(mech=simple_diffusion, partid=None, name=k_diff).
      search_key=(mech=simple_diffusion, partid=NO3, name=k_diff).

    ]