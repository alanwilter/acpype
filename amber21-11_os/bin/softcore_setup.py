#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
#
# Postprocesses a prmtop/rst for use with softcore TI in AMBER 10 (or later, hopefully)
#
# Thomas T. Joseph <thomas.joseph@mssm.edu>
import sys
import re
import os
import subprocess

class Error(Exception):
    pass

def validate(condition, msg="Argh!"):
    """If condition is False, aborts the program with the specified message."""
    if not condition:
        raise Error(msg)

class AmberSystem:
    # POINTERS block indices, stolen from the prmtop format spec on ambermd.org
    NATOM  = 0 # total number of atoms 
    NTYPES = 1 # total number of distinct atom types
    NBONH  = 2 # number of bonds containing hydrogen
    MBONA  = 3 # number of bonds not containing hydrogen
    NTHETH = 4 # number of angles containing hydrogen
    MTHETA = 5 # number of angles not containing hydrogen
    NPHIH  = 6 # number of dihedrals containing hydrogen
    MPHIA  = 7 # number of dihedrals not containing hydrogen
    NHPARM = 8 # currently not used
    NPARM  = 9 # currently not used
    NEXT   = 10 # number of excluded atoms
    NRES   = 11 # number of residues
    NBONA  = 12 # MBONA + number of constraint bonds
    NTHETA = 13 # MTHETA + number of constraint angles
    NPHIA  = 14 # MPHIA + number of constraint dihedrals
    NUMBND = 15 # number of unique bond types
    NUMANG = 16 # number of unique angle types
    NPTRA  = 17 # number of unique dihedral types
    NATYP  = 18 # number of atom types in parameter file, see SOLTY below
    NPHB   = 19 # number of distinct 10-12 hydrogen bond pair types
    IFPERT = 20 # set to 1 if perturbation info is to be read in
    NBPER  = 21 # number of bonds to be perturbed
    NGPER  = 22 # number of angles to be perturbed
    NDPER  = 23 # number of dihedrals to be perturbed
    MBPER  = 24 # number of bonds with atoms completely in perturbed group
    MGPER  = 25 # number of angles with atoms completely in perturbed group
    MDPER  = 26 # number of dihedrals with atoms completely in perturbed groups
    IFBOX  = 27 # set to 1 if standard periodic box, 2 when truncated octahedral
    NMXRS  = 28 # number of atoms in the largest residue
    IFCAP  = 29 # set to 1 if the CAP option from edit was specified
    
    def __init__(self, prmtop_filename, rst_filename=None):
        """Initializes this instance with AMBER prmtop and rst files."""
        self.blocks = {}
        self.block_list = []
        self.formats = {}
        self.format_strings = {}
        self.x, self.y, self.z = [], [], []
        self.box = None
        
        self.load_prmtop(prmtop_filename)
        if rst_filename is not None:
            self.load_rst(rst_filename)
                
    def load_prmtop(self, filename):
        """Loads an AMBER prmtop file into this instance."""
        
        # These regular expressions are used in parsing FLAG and FORMAT lines
        block_name_re = re.compile('%FLAG\s+(\w+)')
        self.format_re = re.compile('%FORMAT\((\d+)(.)([\d\.]+)')
        
        # Read header line, TODO check version of file
        fp = open(filename, "r")
        self.header = fp.readline()
        
        # Read first line first
        line = fp.readline()
        while True:
            if line == "": break
            name = block_name_re.match(line).groups()[0]
            line = fp.readline()
            validate(line != "", "%s truncated after %%FLAG card." % filename)
            format = self.format_re.match(line).groups()
            tokens_per_line = int(format[0])
            token_type = format[1]
            token_length = float(format[2])
            self.formats[name] = (tokens_per_line, token_type, token_length)
            self.format_strings[name] = line
            self.new_block(name, line)
            
            # Read lines until we get one that starts with %, which signifies
            # a new block
            line = fp.readline()
            if line == "": break
            while line[0] is not "%":
                # Subtract 1 from length to ignore newline
                for i in range(0, len(line) - 1, int(token_length)):
                    self.blocks[name].append(line[i:i + int(token_length)])
                line = fp.readline()
                if line == "": break
        
            # Converts tokens to the correct type and saves them in arrays
            # I: integer; a: alphanumeric; E: float
            if token_type == "I":
                self.blocks[name] = [int(x) for x in self.blocks[name]]
            elif token_type == "E":
                self.blocks[name] = [float(x) for x in self.blocks[name]]
        
    def load_rst(self, filename):
        """Loads an AMBER restart file into this instance."""
        
        assert 'POINTERS' in self.blocks, "POINTERS block missing - valid prmtop not loaded"
        
        fp = open(filename, "r")
        fp.readline() # Eat header line
        
        coords_left = self.blocks['POINTERS'][AmberSystem.NATOM] * 3
        validate(coords_left == int(fp.readline()) * 3, \
            "Number of atoms %d in %s differs from that specified in prmtop" \
                % (coords_left, filename))
        
        # There are 6 coordinates, 12 characters each per line
        while coords_left > 0:
            line = fp.readline()
            self.x.append(float(line[0:12]))
            self.y.append(float(line[12:24]))
            self.z.append(float(line[24:36]))
            coords_left -= 3
            if coords_left >= 3:
                self.x.append(float(line[36:48]))
                self.y.append(float(line[48:60]))
                self.z.append(float(line[60:72]))
                coords_left -= 3
                
        validate(coords_left == 0, \
            "Number of coordinates read is not a multiple of 3: is %s corrupt?" % \
            filename)
        
        # If a box is specified in the prmtop, read it in
        if self.blocks['POINTERS'][AmberSystem.IFBOX] != 0:
            self.box = []
            line = fp.readline()
            self.box.append(float(line[0:12]))
            self.box.append(float(line[12:24]))
            self.box.append(float(line[24:36]))
            self.box.append(float(line[36:48]))
            self.box.append(float(line[48:60]))
            self.box.append(float(line[60:72]))
        
    def save_prmtop(self, filename):
        """Saves the prmtop part of this AmberSystem into a prmtop file."""
        
        fp = open(filename, "w")
        fp.write(self.header)
        
        # Eases converting from Fortran to Python/C format specifiers
        type_table = {'I': 'd', 'E': 'E', 'a': 's'}

        # The following ugly code deals with printing tokens in the correct format
        for name in self.block_list:
            fp.write("%%FLAG %s\n" % name)
            fp.write(self.format_strings[name])
            format = self.format_re.match(self.format_strings[name]).groups()
            tokens_per_line = int(format[0])
            token_type = type_table[format[1]]
            token_length = float(format[2])
            token_length_int = int(token_length)
            # Print x.y vs x length specifiers correctly
            if float(token_length_int) == token_length:
                token_length_dec = ""
            else:
                token_length_dec = "%.1f" % (token_length - token_length_int)
                token_length_dec = token_length_dec[1:]
            s = "%%%d%s%s" % (token_length_int, token_length_dec, token_type)

            count = 0
            for x in self.blocks[name]:
                fp.write(s % x)
                count += 1
                if count % tokens_per_line == 0:
                    fp.write("\n")
            
            if count % tokens_per_line != 0:
                fp.write("\n")
        
        fp.close()
        
    def save_rst(self, filename):
        """Saves the coordinates and box information (if present) to a restart file."""
        
        assert len(self.x) == len(self.y) == len(self.z)

        fp = open(filename, "w")
        # Write header
        fp.write("Processed by softcore_setup.py\n")
        fp.write("%6d\n" % len(self.x))

        # Write coordinates
        for i in range(len(self.x)):
            fp.write("%12.7f" % self.x[i])
            fp.write("%12.7f" % self.y[i])
            fp.write("%12.7f" % self.z[i])
            if i % 2 == 1: fp.write("\n")
                
        if len(self.x) % 2 == 1: fp.write("\n")
            
        if self.box is not None:
            for n in self.box:
                fp.write("%12.7f" % n)
            fp.write("\n")
        
        fp.close()
    
    def save_pdb(self, sclist=None, filename=None):
        """Makes a PDB that LEaP will read back in without complaint. This
        differs from ambpdb/ambmask in that it writes beta and occupancy fields
        that are set to 1.0 if that atom is in sclist, and 0.0 if not.
        
        This means:
        - Don't mangle names (e.g. don't do 1H5' instead of H5'1)
        - Use single quotes instead of asterisks
          (which can be accomplished simply by not mangling the prmtop names)
        - Insert a TER record after every MOLECULE
        - Have an END record at the end, I guess
        """
        
        if filename is not None:
            fp = open(filename, "w")
        else:
            fp = sys.stdout
        
        # Iterate over atoms and determine which residue
        # it is in. Save this information in arrays.
        residues = []
        residue_counter = self.num_residues() - 1
        for atom in range(self.num_atoms(), 0, -1):
            residues.append(residue_counter + 1)
            if atom == self.blocks['RESIDUE_POINTER'][residue_counter]:
                residue_counter -= 1
        residues.reverse()
        
        # We count down the number of atoms in each molecule as we go so we
        # can print TER records as needed
        if 'ATOMS_PER_MOLECULE' in self.blocks:
            have_molecules = True
            current_molecule = 0
            atoms_left = self.blocks['ATOMS_PER_MOLECULE'][current_molecule]
        else:
            have_molecules = False
        
        # Print the ATOM records
        for atom in range(self.num_atoms()):
            atom_name = self.blocks['ATOM_NAME'][atom]
            res_name = self.blocks['RESIDUE_LABEL'][residues[atom] - 1]
            
            if sclist is not None and (atom + 1) in sclist:
                is_softcore = 1.0
            else:
                is_softcore = 0.0
            
            fp.write("ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % \
                (atom + 1, atom_name, res_name, residues[atom], \
                self.x[atom], self.y[atom], self.z[atom], \
                is_softcore, is_softcore))
            
            # If there's a ATOMS_PER_MOLECULE block we know where to print the TERs    
            if have_molecules:
                atoms_left -= 1
                if atoms_left == 0:
                    fp.write("TER\n")
                    current_molecule += 1
                    if current_molecule < len(self.blocks['ATOMS_PER_MOLECULE']):
                        atoms_left = self.blocks['ATOMS_PER_MOLECULE'][current_molecule]
        fp.write("END\n")
        
        if filename is not None:
            fp.close()
        
    # The hash we are using to store the blocks does not maintain
    # the order of the keys. We keep an ordered list
    # of keys so when we save the prmtop the blocks will be in that order.
    def new_block(self, name, format_string):
        """Creates a new prmtop block, setting up ancillary data structures as needed."""
        
        if name not in self.blocks:
            assert name not in self.block_list, \
                "AmberSystem block_list inconsistent with blocks"
            self.block_list.append(name)
            self.blocks[name] = []
            self.format_strings[name] = format_string
            
    def residue_for_atom(self, atom):
        """Returns the residue number of a given atom index, or None if not 
        found."""
        
        validate('RESIDUE_POINTER' in self.blocks, \
            "No RESIDUE_POINTER block - Corrupt prmtop")
        if atom < 1 or atom > self.num_atoms(): return None
        
        for residue in range(len(self.blocks['RESIDUE_POINTER']), 0, -1):
            if atom >= self.blocks['RESIDUE_POINTER'][residue - 1]:
                return residue
        
        raise Error("Bizarre residue boundaries - corrupt prmtop?")
        
    def residues_for_atoms(self, atoms):
        """Returns a hash of residue numbers to the provided atom numbers."""
        
        tmp = {}
        for atom in atoms:
            residue = self.residue_for_atom(atom)
            if residue not in tmp:
                tmp[residue] = [atom]
            else:
                tmp[residue].append(atom)
        return tmp
        
    def make_ambmask(self, atoms):
        """Converts a list of atoms into a more human-friendly ambmask string."""
        
        atoms_by_residue = self.residues_for_atoms(atoms)
        mask = ""
        first_residue = True
        for residue in atoms_by_residue:
            if not first_residue: mask += " | "
            first_residue = False
            mask += ":%d@" % residue
            atom_names = [self.blocks['ATOM_NAME'][i - 1].strip() for i in \
                atoms_by_residue[residue]]
            first_atom = True
            for atom_name in sorted(atom_names):
                if not first_atom: mask += ","
                first_atom = False
                mask += atom_name
        return mask
    
    def num_atoms(self):
        """Returns the number of atoms. Must have loaded a prmtop (with
        load_prmtop)."""

        validate('POINTERS' in self.blocks, "Corrupt prmtop")
        return int(self.blocks['POINTERS'][AmberSystem.NATOM])
        
    def num_residues(self):
        """Returns the number of residues."""
        
        validate('POINTERS' in self.blocks, "Corrupt prmtop")
        return int(self.blocks['POINTERS'][AmberSystem.NRES])

def ambmask_to_atom_list(prmtop, rst, ambmask):
    """Runs $AMBERHOME/bin/ambmask, extracts atom indices from its output, and
    returns them as a list"""

    executable = "%s/bin/ambmask" % os.environ['AMBERHOME']
    if os.path.exists(executable) is False:
        raise Error("ambmask executable not found. I checked" \
            " %s. Did you set the AMBERHOME environment variable correctly?" \
            % executable)
    p = subprocess.Popen([executable, '-p', prmtop, '-c', rst, '-find', ambmask],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.communicate()[0]
    if "Error" in out:
        raise Error("There's an error in ambmask %s:\n\n%s" % (ambmask, out))
    # Return the atom indices only
    out = out.split('\n')
    return [int(line[6:11]) for line in out if line[0:4] == 'ATOM']

def display_usage():
    print("""Usage:

To prepare a "perturbed" prmtop/rst for softcore TI:
    softcore_setup.py <prmtop-A> <rst-A> <prmtop-B> <rst-B>
   
    A and B refer to "unperturbed" and "perturbed" states.
   
    New versions of <prmtop-B> and <rst-B> are generated, along with potentially
    useful suggested parameters for sander.

To dump a PDB file to stdout (optionally with specified atoms marked):
    softcore_setup.py <prmtop> <rst> ['ambmask']

    You can edit the resulting PDB to create your "perturbed" structure,
    and use LEaP's loadpdb and saveamberparm commands to create prmtop and
    coordinate files to be further processed as below.
""", file=sys.stderr)   

def make_atom_table(s):
    """Returns a dict with a key for each atom, to make it easy to check whether
    a structure contains an atom of a given name at a given position.
    
    The key contains the atom name and xyz coordinates rounded to the nearest
    tenth of an Angstrom."""

    table = {}
    for i in range(s.num_atoms()):
        key = "%4s%.1f%.1f%.1f" % (s.blocks['ATOM_NAME'][i], s.x[i], s.y[i], s.z[i])
        table[key] = i + 1
    return table

def sclist_too_big_error(sclist, prmtop):
    """Convenience method to print an error when one of the softcore regions is very big, which
    more likely than not indicates that the structures are totally wrong or messed up"""

    raise Error("Predicted a suspiciously large softcore region of %d atoms in %s.\n" \
    "Please check that your structures have the same coordinates except for the\n" \
    "softcore regions." % (len(sclist), prmtop))

def check_block_exists(s, prmtop, name):
    """Check if the named prmtop block exists, but it's not fatal if it doesn't."""

    if name not in s.blocks:
        print("WARNING: %s is missing the %s block. This may cause an error when you\n" \
        "run the simulation. A LEaP with a hacked setBox command can generate it." \
        % (prmtop, name))

def main():
    print("AMBER softcore TI setup utility\n", file=sys.stderr)
    
    # Parse command line options
    if len(sys.argv) == 3:
        print("Dumping a LEaP-friendly PDB to stdout.", file=sys.stderr)
        A = AmberSystem(sys.argv[1], sys.argv[2])
        A.save_pdb()
        sys.exit()
    if len(sys.argv) == 4:
        print("Dumping a LEaP-friendly PDB to stdout. Atoms specified by " \
            "'%s' have beta and\noccupancy values of 1.00," \
            " and all others have beta and occupancy\nvalues of 0.00." % sys.argv[3], file=sys.stderr)
        A = AmberSystem(sys.argv[1], sys.argv[2])
        atom_list = ambmask_to_atom_list(sys.argv[1], sys.argv[2], sys.argv[3])
        A.save_pdb(atom_list)
        sys.exit()
    elif len(sys.argv) == 5:
        (me, prmtop_A, rst_A, prmtop_B, rst_B) = sys.argv
    else:
        display_usage()
        sys.exit()
    
    print("Loading %s and %s..." % (prmtop_A, rst_A), file=sys.stderr)
    A = AmberSystem(prmtop_A, rst_A)
    print("Loading %s and %s..." % (prmtop_B, rst_B), file=sys.stderr)
    B = AmberSystem(prmtop_B, rst_B)
    
    print("Finding softcore regions...", file=sys.stderr)
    
    # Guess the softcore regions
    # sclist_A is atoms not in B; if the softcore regions are really big
    # then the structures are probably messed up
    atoms_A = make_atom_table(A)
    atoms_B = make_atom_table(B)
    sclist_A = [atoms_A[key] for key in atoms_A if key not in atoms_B]
    if len(sclist_A) > 200: sclist_too_big_error(sclist_A, prmtop_A)
    sclist_B = [atoms_B[key] for key in atoms_B if key not in atoms_A]
    if len(sclist_B) > 200: sclist_too_big_error(sclist_B, prmtop_B)
    not_sclist_A = [atoms_A[key] for key in atoms_A if key in atoms_B]
    not_sclist_B = [atoms_B[key] for key in atoms_B if key in atoms_A]

    # Create scmasks that we will report to the user
    scmask_A = A.make_ambmask(sclist_A)
    scmask_B = B.make_ambmask(sclist_B)

    print("""
%s has %d atoms in %d residues: %d softcore atoms, %d other atoms.
%s has %d atoms in %d residues: %d softcore atoms, %d other atoms.""" \
    % (prmtop_A, A.num_atoms(), A.num_residues(), \
    len(sclist_A), len(not_sclist_A), \
    prmtop_B, B.num_atoms(), B.num_residues(), \
    len(sclist_B), len(not_sclist_B)), file=sys.stderr)
        
    assert len(not_sclist_A) == len(not_sclist_B), \
        "Non-softcore regions have different lengths - this is a bug"
    
    # Copy the periodic box
    if 'BOX_DIMENSIONS' in A.blocks:
        B.box = A.box # So the box will be in the rst too
        B.new_block('BOX_DIMENSIONS', A.format_strings['BOX_DIMENSIONS'])
        B.blocks['BOX_DIMENSIONS'] = A.blocks['BOX_DIMENSIONS']
        B.blocks['POINTERS'][AmberSystem.IFBOX] \
            = A.blocks['POINTERS'][AmberSystem.IFBOX]
    else:
        print("No box information was found in %s." % prmtop_A, file=sys.stderr)
    
    # These blocks are important to have...sander seems to want them        
    check_block_exists(A, prmtop_A, 'ATOMS_PER_MOLECULE')
    check_block_exists(B, prmtop_B, 'ATOMS_PER_MOLECULE')
    check_block_exists(A, prmtop_A, 'SOLVENT_POINTERS')
    check_block_exists(B, prmtop_B, 'SOLVENT_POINTERS')
    
    # Come up with sane-looking output filenames
    (prmtop_base_A, prmtop_ext_A) = os.path.splitext(prmtop_A)
    (rst_base_A, rst_ext_A) = os.path.splitext(rst_A)
    (prmtop_base_B, prmtop_ext_B) = os.path.splitext(prmtop_B)
    (rst_base_B, rst_ext_B) = os.path.splitext(rst_B)
    prmtop_out_B = "%s.SC%s" % (prmtop_base_B, prmtop_ext_B)
    rst_out_B = "%s.SC%s" % (rst_base_B, rst_ext_B)
    pdb_out_A = "%s.SC.pdb" % prmtop_base_A
    pdb_out_B = "%s.SC.pdb" % prmtop_base_B

    # Save the prmtop and rst files
    print("\nGenerating %s and %s..." % (prmtop_out_B, rst_out_B), file=sys.stderr)
    B.save_prmtop(prmtop_out_B)
    B.save_rst(rst_out_B)
    
    # Save PDBs with softcore regions "highlighted"
    print("Generating %s and %s..." % (pdb_out_A, pdb_out_B), file=sys.stderr)
    A.save_pdb(sclist_A, pdb_out_A)
    B.save_pdb(sclist_B, pdb_out_B)
    
    print("""
Suggested sander options:
    
 First stage - Charge removal for %s:
    icfe=1, ifsc=0,
    crgmask='%s',
    """ % (prmtop_A, scmask_A))
    print(""" Second stage - Alchemy step for %s:
    icfe=1, ifsc=1,
    crgmask='%s',
    scmask='%s',
    """ % (prmtop_A, scmask_A, scmask_A))
    print(""" Second stage - Alchemy step for %s:
    icfe=1, ifsc=1,
    crgmask='%s',
    scmask='%s',
    """ % (prmtop_out_B, scmask_B, scmask_B))
    print(""" Third stage - Charge restoration for %s:
    icfe=1, ifsc=0,
    crgmask='%s',

IMPORTANT:    
    Please double-check the generated files and parameters before using them.
    You can visually verify the softcore regions by using your favorite PDB 
    viewer to view the generated PDB files and highlighting by temperature
    factor or occupancy.
""" % (prmtop_out_B, scmask_B))


if __name__ == '__main__': main()
    
