# encoding: utf-8

__doc__ = """ 

Intrface to PYMOL with useful pre-processing.

As well various functions use PYMOL as an image generator of
protein structures. 

A key concept in this module is that a viewing frame can be
defined by two residues:

1. The center residue defines the viewing frame where the center 
   residue is in the middle of the screen directly over the
   center of the system.
2. The top residue defines the z-rotation of the above frame such
   that the top residue is directly above the center residue
   in the viewing frame.

This can be used to define reasonably close viewing frames to
any conceivable viewing frame.Several of the functions thus take
center_res and top_res  parameters to define the viewing frame of
reference in PYMOL. 
"""

import os
from . import pdbatoms
from . import v3
from . import util
from . import data
from . import pdbtext
from . import protein


def get_pymol_id_from_res_tag(tag):
  chain_id, res_num, insert = pdbatoms.split_tag(tag)
  if chain_id == " ":
    return "resi %d" % res_num
  else:
    return "(chain %s and resi %d%s)" % (chain_id, res_num, insert)
  
    
# Functions to transform PDB before rendering in PYMOL


def get_scale_max(max_bfactor, upper_bfactor):
  scale_max = 4.0
  if max_bfactor is not None and upper_bfactor is not None:
    if max_bfactor < upper_bfactor:
      scale_max *= max_bfactor / upper_bfactor
      if scale_max < 1:
        scale_max = 1
  return scale_max

  
def add_fake_water_atom(soup, res_type, bfactor):
  dummy_atom = pdbatoms.Atom()
  dummy_atom.pos = soup.atoms()[0].pos.copy()
  dummy_atom.type = "O"
  dummy_atom.bfactor = bfactor
  dummy_res = pdbatoms.Residue(res_type, '', 9999)
  dummy_res.insert_atom(dummy_atom)
  soup.append_residue(dummy_res)


def rescale_positive_bfactors_pdb(pdb, lower_bfactor, upper_bfactor):
  """
  Returns max_bfactor after rescale (needed for worm
  calculation)
  """
  soup = pdbatoms.Soup(pdb)
  bfactors = [a.bfactor for a in soup.atoms()]
  # cut-off max_values
  if upper_bfactor:
    bfactors = [upper_bfactor if b > upper_bfactor else b 
                for b in bfactors]
    # will delete later within pymol script
    add_fake_water_atom(soup, 'XXX', upper_bfactor)
  # cut-off below min_val to zero
  if lower_bfactor:
    for j in range(len(bfactors)):
      bfactors = [0 if b < lower_bfactor else b for b in bfactors]
  for a, bfactor in zip(soup.atoms(), bfactors):
    a.bfactor = bfactor
  new_pdb = util.fname_variant(pdb)
  soup.write_pdb(new_pdb)
  return new_pdb, max(bfactors)


def rescale_positive_bfactors_of_pdbs(pdbs, lower_bfactor, upper_bfactor):
  "Returns list of new_pdbs, and max_bfactor"
  max_bfactor = 0
  new_pdbs = []
  for pdb in pdbs:
    new_pdb, this_max_bfactor = rescale_positive_bfactors_pdb(
      pdb, lower_bfactor, upper_bfactor)
    if this_max_bfactor > max_bfactor:
      max_bfactor = this_max_bfactor
    new_pdbs.append(new_pdb)
  return new_pdbs, max_bfactor


def rescale_positive_negative_bfactors_pdb(
    pdb, lower_bfactor, upper_bfactor):
  """
  Returns max_bfactor after rescale 
  """
  soup = pdbatoms.Soup(pdb)
  bfactors = [a.bfactor for a in soup.atoms()]
  if upper_bfactor is None:
    upper_bfactor = max(bfactors)
  # cut-off max_values
  if upper_bfactor:
    for j in range(len(bfactors)):
      if bfactors[j] > upper_bfactor:
        bfactors[j] = upper_bfactor
      if bfactors[j] < -upper_bfactor:
        bfactors[j] = -upper_bfactor
    # will delete later within pymol script
    add_fake_water_atom(soup, 'XXX', upper_bfactor)
    add_fake_water_atom(soup, 'XXX', -upper_bfactor)
  # cut-off below min_val to zero
  if lower_bfactor:
    for j in range(len(bfactors)):
      if -lower_bfactor < bfactors[j] < lower_bfactor:
        bfactors[j] = 0.0
  for a, bfactor in zip(soup.atoms(), bfactors):
    a.bfactor = bfactor
  new_pdb = util.fname_variant(pdb)
  soup.write_pdb(new_pdb)
  return new_pdb, max(bfactors)


def rescale_positive_negative_bfactor_pdbs(
    pdbs, lower_bfactor, upper_bfactor):
  "Returns list of new_pdbs, and max_bfactor"
  max_bfactor = 0
  new_pdbs = []
  for pdb in pdbs:
    new_pdb, this_max_bfactor = \
        rescale_positive_negative_bfactors_pdb(
            pdb, lower_bfactor, upper_bfactor)
    if this_max_bfactor > max_bfactor:
      max_bfactor = this_max_bfactor
    new_pdbs.append(new_pdb)
  return new_pdbs, max_bfactor


# PYMOL script functions and snippets

def make_bgcolor_script(bg_color):
  return  "cmd.bg_color('%s');\n" % bg_color

  
def make_load_pdbs_script(pdbs):
  "Returns pymol script, name of pdbs"
  script = ""
  for pdb in pdbs:
    script += "load %s\n" % pdb
  script += "hide everything\n"
  return script
  

def make_separate_chain_colors_script(pdbs):
  names = [os.path.basename(p).replace('.pdb', '') 
           for p in pdbs]
  colors = ['util.color_chains("(%s and elem c)")\n' % n 
            for n in names]
  return ''.join(colors)


red_white_gradient_script = """\
run %s
color_b all, gradient=wr
""" % os.path.join(data.data_dir, "color_b.py")


blue_white_red_gradient_script = """\
cmd.spectrum("b", "blue_white_red", selection="all");
"""


cartoon_script = """
set cartoon_flat_sheets, 0
set cartoon_rect_width, 0.2
set cartoon_oval_width, 0.2
set cartoon_loop_radius, 0.2
set cartoon_tube_radius, 0.2
cartoon auto
show cartoon
"""


putty_script = """
set cartoon_flat_sheets, 0
set cartoon_putty_scale_max, %(scale_max)f
set cartoon_putty_radius, 0.4
set cartoon_putty_scale_power, 1
cartoon putty
show cartoon
"""


def make_putty_script(scale_max):
  return putty_script % { 'scale_max': scale_max }


def make_sticks_above_bfactor_script(lower_bfactor):
  script = ""
  script += "select hot, b > %f or b < -%f\n" % \
             (lower_bfactor, lower_bfactor)
  script += "select cold, b < %f and b > -%f\n" % \
             (lower_bfactor, lower_bfactor)
  script += "show stick, hot\n"
  script += "hide stick, cold\n"
  script += "deselect\n"
  return script
    
    
def make_ligands_as_sticks_script(pdbs, color=""):
  script = ""
  for pdb in pdbs:
    name = os.path.basename(pdb).replace('.pdb', '')
    soup = pdbatoms.Soup(pdb)
    for res in soup.residues():
      if res.type not in data.res_name_to_char:
        if res.type not in "HOH":
          chain_id_script = ""
          if res.chain_id.strip():
            chain_id_script = "and chain %s" % res.chain_id
          script += \
              "show stick, %s %s and resn %s and resi %d\n" \
                % (name, chain_id_script, res.type, res.num)
          if color:
            script += \
              "color %s, %s %s and resn %s and resi %d\n" \
                % (color, name, chain_id_script, res.type, res.num)
  script += "show nonbonded\n"
  return script


highlight_res_script = """
select highlight, %s
show stick, highlight
color green, highlight
"""


def make_highlight_res_script(pymol_res_id):
  return highlight_res_script % pymol_res_id


peptide_style_script = """
select bb, name ca+n+h+o+c+oxt+h1+h2+h3+ch3+hh31+hh32+hh33+3hh3+2hh3+1hh3
show sphere, bb
hide stick, bb
util.cbaw bb
select sc, not bb and not hydro
show stick, sc
hide sphere, sc
util.cbag sc
set sphere_quality, 2
"""


hide_backbone_sticks_script = """
hide stick, hydro
select bb, name c+o+n+h+oxt
hide stick, bb
select nuc, resn A+U+T+C+G+A3+U3+T3+C3+G3+A5+U5+T5+C5+G5+DA+DT+DC+DG
select nuc_bb, name P+O1P+O2P+OP1+Op2+O3'+C3'+C2'+C1'+O4'+C4'+C5'+O5'
hide cartoon, nuc and not nuc_bb
deselect
"""


def bfactor_script(pdb, lower_bfactor=None, upper_bfactor=None, 
                   max_bfactor=None, is_putty=False):
  "Returns script that displays bfactors of pdb" 
  script = make_load_pdbs_script([pdb])
  script += red_white_gradient_script
  if is_putty:
    script += make_putty_script(get_scale_max(max_bfactor, upper_bfactor))
  else:
    script += cartoon_script
    script += "cartoon tube\n"
  if not is_putty:
    if lower_bfactor is not None:
      script += make_sticks_above_bfactor_script(lower_bfactor)
    else:
      script += "show stick\n"
  script += make_ligands_as_sticks_script([pdb])
  script += hide_backbone_sticks_script
  return script



def split_resname(resname):
  "Returns (chain_id, res_num)"
  words = resname.split(":")
  if len(words) == 2:
    return (words[0], int(words[1]))
  else:
    return (' ', int(words[0]))


def find_ca_of_resname(atoms, resname):
  chain_id, res_num = split_resname(resname)
  for atom in atoms:
    if chain_id == atom.chain_id and res_num == atom.res_num:
      if "CA" == atom.type:
        return atom
  return None


def get_pdb_transform(pdb, center_res, top_res):
  """
  Returns a transformation matrix that centers pdb to 
  center_res on the z-axis and moves top_res above center_res
  on the y-axis
  """
  soup = pdbatoms.Soup(pdb)
  atoms = soup.atoms()
  soup_center = pdbatoms.get_center(atoms)
  translation = v3.translation(-soup_center)
  soup.transform(translation)
  result = translation

  center_atom = find_ca_of_resname(soup.atoms(), center_res)
  view = v3.vector(0, 0, 1)
  axis = v3.cross(view, center_atom.pos)
  angle = v3.vec_dihedral(view, axis, center_atom.pos)
  rotation = v3.rotation(axis, angle)
  soup.transform(rotation)
  result = v3.combine(rotation, result)

  top_atom = find_ca_of_resname(soup.atoms(), top_res)
  top_dir = v3.vector(0, 1, 0)
  axis = view.copy()
  angle = v3.vec_dihedral(top_dir, axis, top_atom.pos)
  rotation2 = v3.rotation(axis, angle)
  result = v3.combine(rotation2, result)
  
  del soup
  return result



def run_pymol_script(pml, width=500, height=500):
  is_quit = 'quit' in util.words_in_file(pml)
  if is_quit:
    pymol_batch = data.binary("pymol_batch")
    cmd = pymol_batch + ' -c '
  else:
    pymol = data.binary("pymol")
    cmd = pymol + " -q " # no splash screen
  cmd += " -W %d -H %d " % (width, height)
  cmd += pml
  util.run_with_output(cmd)


# Functions to generate PNG's with PYMOL


def pdb_to_bfactor_png(
    bfactor_pdb, png, lower_bfactor=None, upper_bfactor=None, 
    highlight_res=None, is_putty=False, is_sticks=True,
    center_res=None, top_res=None, height=480, width=480):
  """
  Generates a bfactor-colored .png using the standard white-to-red
  color scheme with a useful set of options. highlight_res, 
  center_res & top_res follow the 'A:10' style of residue naming.
  """

  pdb, max_bfactor = rescale_positive_bfactors_pdb(
      bfactor_pdb, lower_bfactor, upper_bfactor)
  temp_fnames = [pdb]

  if center_res or top_res:
    pdb = get_pdb_transform(pdb, center_res, top_res)
    temp_fnames.append(pdb)

  script = ""

  script += bfactor_script(
       pdb, lower_bfactor, upper_bfactor,
       max_bfactor, is_putty)

  if highlight_res is not None:
    pymol_res_id = get_pymol_id_from_res_tag(highlight_res)
    script += make_highlight_res_script(pymol_res_id)
    script += hide_backbone_sticks_script

  script += "clip far, -20\n"
  script += "save %s\n" % png
  script += "quit"

  pml = util.fname_variant('temp.pml')
  open(pml, 'w').write(script)
  run_pymol_script(pml, width, height)
  temp_fnames.append(pml)

  util.clean_fname(*temp_fnames)


def soup_to_bfactor_png(
    soup, png, bfactors, lower_bfactor=None, upper_bfactor=None,
    highlight_res=None, is_putty=False, is_sticks=True,
    center_res=None, top_res=None, height=480, width=480):
  """
  Wrapper around pdb_to_bfactor that loads in external values 
  for the residue bfactors with a reusable soup object.
  """
  pdb = util.temp_fname('.pdb')
  temp_fnames = [pdb]
  soup.load_residue_bfactors(bfactors)
  soup.write_pdb(pdb)
  pdb_to_bfactor_png(
    pdb, png, lower_bfactor, upper_bfactor,
    highlight_res, is_putty, is_sticks,
    center_res, top_res, height, width)
  util.clean_fname(*temp_fnames)


def make_pdbs_png(
    png, pdbs, bgcolor="white", center_res=None, top_res=None,
    highlight_res=None, is_sticks=True, is_putty=False,
    width=480, height=480):

  if 'transparent' in bgcolor:
    script = 'set opaque_background, off\n'
  else: 
    script = make_bgcolor_script(bgcolor)

  temp_fnames = []
  if center_res and top_res:
    transform = get_pdb_transform(pdbs[0], center_res, top_res)
    for i in range(len(pdbs)):
      soup = pdbatoms.Soup(pdbs[i])
      soup.transform(transform)
      new_pdb = util.fname_variant(pdbs[i])
      soup.write_pdb(new_pdb)
      temp_fnames.append(new_pdb)
      pdbs[i] = new_pdb
      del soup

  script += make_load_pdbs_script(pdbs)
  script += make_separate_chain_colors_script(pdbs)

  if is_putty:
    script += make_putty_script(get_scale_max(
        max_bfactor, upper_bfactor))
  else:
    script += cartoon_script

  if not is_sticks:
    script += "hide stick\n"
  else:
    script += "show stick\n"

  script += make_ligands_as_sticks_script(pdbs)

  if highlight_res:
    script += make_highlight_res_script(highlight_res)
  script += hide_backbone_sticks_script

  # script += "clip far, 5\n"
  script += "save %s\n" % png
  script += "quit"

  pml = util.fname_variant('temp.pml')
  open(pml, 'w').write(script)
  run_pymol_script(pml, width, height)
  temp_fnames.append(pml)

  util.clean_fname(*temp_fnames)




