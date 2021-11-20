default_force_field = """
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p
source leaprc.gaff2
"""

leap_template = """
{force_fields}
{more_force_fields}
x = loadpdb {input_pdb}
{box_info}
{more_leap_cmds}
set default PBradii mbondi3
set default nocenter on
saveAmberParm x {prmtop} {rst7}
quit
"""

water_str = "HETATM  607  O   HOH A 167      12.812  27.696   7.842  1.00 56.38           O"
