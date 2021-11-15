import subprocess
import parmed
from .template import default_force_field, leap_template
from .utils import easy_call


def _make_leap_template(parm,
                        ns_names,
                        gaplist,
                        sslist,
                        input_pdb,
                        prmtop='prmtop',
                        rst7='rst7'):
    # box
    box = parm.box
    if box is not None:
        box_info = 'set x box { %s  %s  %s }' % (box[0], box[1], box[2])
    else:
        box_info = ''

    # Now we can assume that we are dealing with AmberTools16:
    more_force_fields = ''

    for res in ns_names:
        more_force_fields += '%s = loadmol2 %s.mol2\n' % (res, res)
        more_force_fields += 'loadAmberParams %s.frcmod\n' % res

    #  more_leap_cmds 
    more_leap_cmds = ''
    if gaplist:
        for d, res1, resid1, res2, resid2 in gaplist:
            more_leap_cmds += 'deleteBond x.%d.C x.%d.N\n' % (resid1, resid2)

    #  process sslist
    if sslist:
        for resid1, resid2 in sslist:
            more_leap_cmds += 'bond x.%d.SG x.%d.SG\n' % (resid1+1, resid2+1)

    leap_string = leap_template.format(
        force_fields=default_force_field,
        more_force_fields=more_force_fields,
        box_info=box_info,
        input_pdb=input_pdb,
        prmtop=prmtop,
        rst7=rst7,
        more_leap_cmds=more_leap_cmds)
    return leap_string
