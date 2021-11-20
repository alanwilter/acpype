#
# Module with global parameters for mm_pbsa
#
# Holger Gohlke: 17.04.2002
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_global;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(
                  @GEN_PAR
                  @DC_PAR
                  @GC_PAR @GC_FIL 
                  @AS_PAR @AS_FIL
                  @MM_PAR @MM_FIL
                  @GB_PAR @GB_FIL
                  @PB_PAR 
                  @DE_PAR @DE_FIL 
                  @RL_PAR @RL_FIL
                  @MS_PAR @MS_FIL 
                  @NM_PAR @NM_FIL
                  $HTMLPATH
                 );
@EXPORT_OK   = qw();
$VERSION     = 1.00;

########################################################################

use strict;

########################################################################

# Declaration of global fields
##############################
use vars qw(
            @GEN_PAR
            @DC_PAR
            @GC_PAR @GC_FIL 
            @AS_PAR @AS_FIL
            @MM_PAR @MM_FIL
            @GB_PAR @GB_FIL
            @PB_PAR 
	    @DE_PAR @DE_FIL 
	    @RL_PAR @RL_FIL
            @MS_PAR @MS_FIL 
            @NM_PAR @NM_FIL
            $HTMLPATH
           );

# Definition of global fields that contain info
#   which parameters and files are needed on input
##################################################
@GEN_PAR = ("PREFIX", "PATH", "COMPLEX", "RECEPTOR", "LIGAND",
            "GC", "AS", "DC", "MM", "GB", "PB", "MS", "NM"); # COMPT, RECPT, LIGPT,
                                                             # START, STOP, OFFSET
                                                             #   are not included here.

@DC_PAR  = ("DCTYPE", "COMREC", "COMLIG", "COMPRI", 
                      "RECRES", "RECPRI", "RECMAP", 
                      "LIGRES", "LIGPRI", "LIGMAP");

@GC_PAR  = ("BOX", "NTOTAL", "NSTART", "NSTOP", "NFREQ",
            "NUMBER_LIG_GROUPS", "NUMBER_REC_GROUPS");
@GC_FIL  = ("MAKE_CRD");

@AS_PAR  = ("BOX", "NTOTAL", "NSTART", "NSTOP", "NFREQ",
            "NUMBER_LIG_GROUPS", "NUMBER_REC_GROUPS",
            "NUMBER_MUTANT_GROUPS");
@AS_FIL  = ("MAKE_CRD");

@MM_PAR  = ("DIELC");
@MM_FIL  = ("SANDER");

@NM_PAR  = ("PROC", "MAXCYC", "DRMS", "IGB");
@NM_FIL  = ("NMODE", "SANDER", "NABNMODE");

@PB_PAR  = ("PROC", "REFE", "INDI", "EXDI", "SCALE", "LINIT");
@DE_PAR  = ("FOCUS", "PERFIL", "SURFTEN", "SURFOFF");
@DE_FIL  = ("CHARGE", "SIZE", "DELPHI", "AMBPDB");
@RL_PAR  = ("RADIOPT", "ARCRES");
@RL_FIL  = ("PBSA");

@GB_PAR  = ("IGB", "GBSA", "SALTCON", "EXTDIEL", "INTDIEL", 
            "SURFTEN", "SURFOFF");
@GB_FIL  = ("SANDER");

@MS_PAR  = ("PROBE");
@MS_FIL  = ("MS", "AMBPDB");

$HTMLPATH = "http://ambermd.org/Questions/mm_pbsa.html";

########################################################################

1; # Necessary for package function
