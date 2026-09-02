.. ACPYPE documentation master file, created by
   sphinx-quickstart on Tue Dec 28 23:14:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ACPYPE's documentation!
==================================


**ACPYPE - AnteChamber PYthon Parser interfacE**

A tool based in **Python** to use **Antechamber** to generate topologies for chemical
compounds and to interface with others python applications like CCPN and ARIA.
Topologies files to be generated so far: CNS/XPLOR, GROMACS, CHARMM and AMBER.

Quick start
-----------

``acpype`` has two modes, and both are a single command.

**From a small molecule** -- a ``.mol2``, ``.pdb``, ``.mdl``/``.mol`` file, or a SMILES
string -- it drives the whole AmberTools pipeline:

.. code-block:: bash

   acpype -i molecule.mol2 -b MOL -c bcc -n 0 -a gaff2

That single command does the work of three separate AmberTools runs -- ``antechamber``
to assign atom types and charges, ``parmchk2`` to fill in missing parameters, and
``tleap`` to build the topology -- and then also writes GROMACS, CNS/XPLOR and CHARMM
versions into one ``MOL.acpype/`` folder.

**From existing AMBER files**, converting a LEaP ``prmtop``/``inpcrd`` pair to GROMACS
without needing AmberTools at all:

.. code-block:: bash

   acpype -p FFF_AC.prmtop -x FFF_AC.inpcrd

Run ``acpype -h`` for every option and a description of each output file. The
`README <https://github.com/alanwilter/acpype#readme>`_ covers installation and
worked examples for GROMACS, CNS/XPLOR and NAMD.

Using ACPYPE as a library
-------------------------

``acpype.acs_api.acpype_api`` returns the generated topologies as a JSON string, and
``acpype.topol.ACTopol`` exposes the same pipeline as an object:

.. code-block:: python

   from acpype.topol import ACTopol

   molecule = ACTopol("molecule.mol2", chargeType="bcc", atomType="gaff2", basename="MOL")
   molecule.createACTopol()
   molecule.createMolTopol()

The full API is documented below.

API reference
-------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

   acpype

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
