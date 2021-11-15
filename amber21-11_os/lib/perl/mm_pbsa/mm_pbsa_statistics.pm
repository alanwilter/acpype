#
# Module to calc statistics on mm_pbsa output
#
# Holger Gohlke: 18.10.2001
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_statistics;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(calc_stat mm_pbsa_stat);
@EXPORT_OK   = qw($TEMP $R $gammaP $betaP $gammaG $betaG);
$VERSION     = 1.00;

########################################################################

use strict;
use mm_pbsa_global qw($HTMLPATH);

########################################################################

# Declaration of global variables
#################################
use vars qw(
            $TEMP
            $R
            $gammaP
            $betaP
            $gammaG
            $betaG

            $CALC
            $DATA
            $AVG
            $STD
            $DEC
            $VAR
            $VARREG
            $VARDEC
            $VARPAI
           );

# Definitions of physical parameters
####################################
$TEMP   = 300.0;
$R      = 8.314;

# Definition of a reference to an array that
#   contains information about calculations
#   done by mm_pbsa
############################################
$CALC = [];

# Definitions of references to 
#   DATA, AVeraGe and STandarDeviation arrays.
# Each entry is a reference to a hash, i.e.
#   $DATA->[file order no.]->{calc type}->{variable}->[snapshot]     = value
#   $AVG ->[file order no.]->{calc type}->{variable}->[total/residue]= value
#   $STD ->[file order no.]->{calc type}->{variable}->[total/residue]= value
# Entries with indices > no. of files are used
#   to store info like TOTAL = COM - REC - LIG
############################################################################
$DATA = [];
$AVG  = [];
$STD  = [];

# Definition of a reference to DEC array.
# Each entry is a reference to an array that contains
#   mappings of start and stop residues for energy decomposition, i.e.
#   $DEC->[file order no.]->[2j]   = start value for first group 
#   $DEC->[file order no.]->[2j+1] = end   value for first group 
#   and so forth.
# "Mapping" is understood as follows:
#   Residues with indices 0 .. indmax in $DATA->[$i+1 (or 2)] 
#     of receptor (or ligand) are mapped to
#     $DEC->[$i+1 (or 2)]->[2j] .. $DEC->[$i+1 (or 2)]->[2j+1].
#   The same is done for the complex:
#     $DATA->[$i+0] is mapped to 
#     $DEC->[$i+0]->[2j] .. $DEC->[$i+0]->[2j+1].
#   Hence, set union of $DEC->[$i+0] equals 
#          set union of ($DEC->[$i+1] v $DEC->[$i+2]).
# The user has to ensure, that the mapping between residues in different
#   files is correct. No automatic correction of this can be performed.
############################################################################
$DEC = [];

# Definitions of references to global hashs that 
#   contain information about data parameters.
# For each key, the value is a reference to a hash.
#
# This hash contains -
#   as keys:
#     variable names
#   as values:
#     references to arrays with
#       index = 0: value = 0 -> skip parameter
#                  value = 1 -> read parameter
#                  value = 2 -> calc parameter
#                  value = 3 -> special treatment
#       index = 1: value = 0 -> no output
#                  value > 0 -> output (number determines order)
#       index = 2: if value[0] = 1 -> regular expression for reading
#                  if value[0] = 2 -> expression for calculation 
#                                     (a "?" indicates that this data is not mandatory)
#
# New calculation modes (i.e. MM, PB, ...) as well as new variables
#   (i.e. VDWNB, ELENB, ...) may be added w/o further needs to
#   any changes in the code.
####################################################################
$VARREG = {
           "MM" => { 
                     "VDWNB" => [1,0,'VDWAALS = +(-?\d+\.\d+)'],
                     "ELENB" => [1,0,'EEL     = +(-?\d+\.\d+)'],
                     "BOND"  => [1,0,'BOND    = +(-?\d+\.\d+)'],
                     "ANGLE" => [1,0,'ANGLE   = +(-?\d+\.\d+)'],
                     "DIHED" => [1,0,'DIHED      = +(-?\d+\.\d+)'],
                     "VDW14" => [1,0,'1-4 VDW = +(-?\d+\.\d+)'],
                     "ELE14" => [1,0,'1-4 EEL = +(-?\d+\.\d+)'],

                     "ELE"   => [2,1,'+MM_ELENB+MM_ELE14'],
                     "VDW"   => [2,2,'+MM_VDWNB+MM_VDW14'],
                     "INT"   => [2,3,'+MM_BOND+MM_ANGLE+MM_DIHED'],
                     "GAS"   => [2,4,'+MM_ELE+MM_VDW+MM_INT']
                   },

           "PB" => {
                     "PB"       => [1,0,'corrected reaction field energy: +(-?\d+\.\d+)'],
                     "PBDIS"    => [1,0,'EDISPER = +(-?\d+\.\d+)'],
                     "PBCAV"    => [1,0,'ECAVITY = +(\d+\.\d+)'],
                     "PBSUR"    => [3,5],
                     "PBCAL"    => [3,6],
                     "PBSOL"    => [2,7,'+PB_PBCAL+PB_PBSUR+PB_PBDIS'],
                     "PBELE"    => [2,8,'+PB_PBCAL+?MM_ELE'],
                     "PBTOT"    => [2,9,'+PB_PBSOL+?MM_GAS']
                   },

           "PBHyb" => {
		    "ELRAELE"  => [1,0,'ELRAELE = +(-?\d+\.\d+)'],
		    "EPB"      => [1,0,'EPB = +(-?\d+\.\d+)'],
                    "PB"       => [1,0,'corrected reaction field energy: +(-?\d+\.\d+)'],
                    "PBDIS"    => [1,0,'EDISPER = +(-?\d+\.\d+)'],
                    "PBCAV"    => [1,0,'ECAVITY = +(\d+\.\d+)'],
                    "PBSUR"    => [3,5],
                    "PBCAL"    => [3,6],
                    "PBSOL"    => [2,7,'+PBHyb_PBCAL+PBHyb_PBSUR+PBHyb_PBDIS'],
                    "PBELE"    => [2,8,'+PBHyb_PBCAL'],
                    "PBTOT"    => [2,9,'+PBHyb_PBSOL']
                  },

           "GB" => { 
                     "GB"    => [1,11,'EGB        = +(-?\d+\.\d+)'],
                     "GBSUR" => [3,10],
                     "GBSOL" => [2,12,'+GB_GB+GB_GBSUR'],
                     "GBELE" => [2,13,'+GB_GB+MM_ELE'],
                     "GBTOT" => [2,14,'+GB_GBSOL+MM_GAS']
                   },

           "MS" => { 
                     "SURF"  => [1,0,'surface area = +(\d+\.\d+)']
                   },

           "NM" => { 
                     "STRA"  => [1,0,'translational +(\d+\.\d+)'],
                     "SROT"  => [1,0,'rotational +(\d+\.\d+)'],
                     "SVIB"  => [1,0,'vibrational +(\d+\.\d+)'],
                     "STOT"  => [1,0,'Total +(\d+\.\d+)'],
                     "TSTRA" => [3,15],
                     "TSROT" => [3,16],
                     "TSVIB" => [3,17],
                     "TSTOT" => [3,18],
                   },
          };

$VARDEC = {
           "MM" => { 
                                     # TDC   1   1.0
                     "TINT"  => [1,3,'^TDC +\d+ +(-?\d+\.\d+)'],
                     "SINT"  => [1,1,'^SDC +\d+ +(-?\d+\.\d+)'],
                     "BINT"  => [1,2,'^BDC +\d+ +(-?\d+\.\d+)'],

                                     # TDC   1   1.0          2.0
                     "TVDW"  => [1,6,'^TDC +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SVDW"  => [1,4,'^SDC +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BVDW"  => [1,5,'^BDC +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                                     # TDC   1   1.0         2.0          3.0
                     "TELE"  => [1,9,'^TDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SELE"  => [1,7,'^SDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BELE"  => [1,8,'^BDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TGAS"  => [2,12,'+MM_TINT+MM_TVDW+MM_TELE'],
                     "SGAS"  => [2,10,'+MM_SINT+MM_SVDW+MM_SELE'],
                     "BGAS"  => [2,11,'+MM_BINT+MM_BVDW+MM_BELE']
                   },

           "GB" => { 
                                     # TDC   1   1.0         2.0         3.0          4.0
                     "TGB"   => [1,15,'^TDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SGB"   => [1,13,'^SDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BGB"   => [1,14,'^BDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TGBSUR"=> [3,18],
                     "SGBSUR"=> [3,16],
                     "BGBSUR"=> [3,17],

                     "TGBSOL"=> [2,21,'+GB_TGB+GB_TGBSUR'],
                     "SGBSOL"=> [2,19,'+GB_SGB+GB_SGBSUR'],
                     "BGBSOL"=> [2,20,'+GB_BGB+GB_BGBSUR'],

                     "TGBELE"=> [2,0,'+GB_TGB+MM_TELE'],
                     "SGBELE"=> [2,0,'+GB_SGB+MM_SELE'],
                     "BGBELE"=> [2,0,'+GB_BGB+MM_BELE'],

                     "TGBTOT"=> [2,24,'+GB_TGBSOL+MM_TGAS'],
                     "SGBTOT"=> [2,22,'+GB_SGBSOL+MM_SGAS'],
                     "BGBTOT"=> [2,23,'+GB_BGBSOL+MM_BGAS']
                   },                 

           "PB" => {
                                      # PB_TDC   1   1.0         2.0         3.0          4.0
                     "TPB"   => [1,27,'^PB_TDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SPB"   => [1,25,'^PB_SDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BPB"   => [1,26,'^PB_BDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TPBSOL"=> [2,30,'+PB_TPB+GB_TGBSUR'],
                     "SPBSOL"=> [2,28,'+PB_SPB+GB_SGBSUR'],
                     "BPBSOL"=> [2,29,'+PB_BPB+GB_BGBSUR'],

                     "TPBELE"=> [2,0,'+PB_TPB+MM_TELE'],
                     "SPBELE"=> [2,0,'+PB_SPB+MM_SELE'],
                     "BPBELE"=> [2,0,'+PB_BPB+MM_BELE'],

                     "TPBTOT"=> [2,33,'+PB_TPBSOL+MM_TGAS'],
                     "SPBTOT"=> [2,31,'+PB_SPBSOL+MM_SGAS'],
                     "BPBTOT"=> [2,32,'+PB_BPBSOL+MM_BGAS']
                   },

           "MS" => { 
                                     # TDC   1   1.0         2.0         3.0         4.0          5.0
                     "TSURF" => [1,0,'^TDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SSURF" => [1,0,'^SDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BSURF" => [1,0,'^BDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)']
                   },
           "NM" => { 
                     "T_SVIB"  => [1,0,'^TDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "S_SVIB"  => [1,0,'^SDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "B_SVIB"  => [1,0,'^BDC +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "T_TSVIB"  => [3,3,],
                     "S_TSVIB"  => [3,1,],
                     "B_TSVIB"  => [3,2,]
                   },

          };

$VARPAI = {
           "MM" => { 
                                     # TDC   1         1.0
                     "TINT"  => [1,3,'^TDC +\d+-> +\d+ +(-?\d+\.\d+)'],
                     "SINT"  => [1,1,'^SDC +\d+-> +\d+ +(-?\d+\.\d+)'],
                     "BINT"  => [1,2,'^BDC +\d+-> +\d+ +(-?\d+\.\d+)'],

                                     # TDC   1        1.0          2.0
                     "TVDW"  => [1,6,'^TDC +\d+-> +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SVDW"  => [1,4,'^SDC +\d+-> +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BVDW"  => [1,5,'^BDC +\d+-> +\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                                     # TDC   1        1.0         2.0          3.0
                     "TELE"  => [1,9,'^TDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SELE"  => [1,7,'^SDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BELE"  => [1,8,'^BDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TGAS"  => [2,12,'+MM_TINT+MM_TVDW+MM_TELE'],
                     "SGAS"  => [2,10,'+MM_SINT+MM_SVDW+MM_SELE'],
                     "BGAS"  => [2,11,'+MM_BINT+MM_BVDW+MM_BELE']
                   },

           "GB" => { 
                     
                                     # TDC   1        1.0         2.0         3.0          4.0
                     "TGB"   => [1,15,'^TDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SGB"   => [1,13,'^SDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BGB"   => [1,14,'^BDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TGBSUR"=> [3,18],
                     "SGBSUR"=> [3,16],
                     "BGBSUR"=> [3,17],

                     "TGBSOL"=> [2,21,'+GB_TGB+GB_TGBSUR'],
                     "SGBSOL"=> [2,19,'+GB_SGB+GB_SGBSUR'],
                     "BGBSOL"=> [2,20,'+GB_BGB+GB_BGBSUR'],

                     "TGBELE"=> [2,0,'+GB_TGB+MM_TELE'],
                     "SGBELE"=> [2,0,'+GB_SGB+MM_SELE'],
                     "BGBELE"=> [2,0,'+GB_BGB+MM_BELE'],

                     "TGBTOT"=> [2,24,'+GB_TGBSOL+MM_TGAS'],
                     "SGBTOT"=> [2,22,'+GB_SGBSOL+MM_SGAS'],
                     "BGBTOT"=> [2,23,'+GB_BGBSOL+MM_BGAS']
                   },                 

           "PB" => {
                                      # PB_TDC   1       1.0         2.0         3.0          4.0
                     "TPB"   => [1,27,'^PB_TDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SPB"   => [1,25,'^PB_SDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BPB"   => [1,26,'^PB_BDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],

                     "TPBSOL"=> [2,30,'+PB_TPB+GB_TGBSUR'],
                     "SPBSOL"=> [2,28,'+PB_SPB+GB_SGBSUR'],
                     "BPBSOL"=> [2,29,'+PB_BPB+GB_BGBSUR'],

                     "TPBELE"=> [2,0,'+PB_TPB+MM_TELE'],
                     "SPBELE"=> [2,0,'+PB_SPB+MM_SELE'],
                     "BPBELE"=> [2,0,'+PB_BPB+MM_BELE'],

                     "TPBTOT"=> [2,33,'+PB_TPBSOL+MM_TGAS'],
                     "SPBTOT"=> [2,31,'+PB_SPBSOL+MM_SGAS'],
                     "BPBTOT"=> [2,32,'+PB_BPBSOL+MM_BGAS']
                   },

           "MS" => { 
                                     # TDC   1        1.0         2.0         3.0         4.0          5.0
                     "TSURF" => [1,0,'^TDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "SSURF" => [1,0,'^SDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)'],
                     "BSURF" => [1,0,'^BDC +\d+-> +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)']
                   },

          };
########################################################################################
########################################################################################

sub calc_stat(){
################

  # Parameters: \%GEN,\%DEC,\%DEL,\%GBO,\%PRO
  my $r_gen = shift;
  my $r_dec = shift;
  my $r_del = shift;
  my $r_gbo = shift;
  my $r_pro = shift;

  print "\n";
  print "=>> Doing statistics\n";

  my $statin = $r_gen->{"PREFIX"} . "_statistics.in";
  open(OUT,">$statin") || die("$statin not opened\nFor details see: $HTMLPATH#calc_stat_in_not_opened\n");
  if($r_gen->{"COMPLEX"}){
    print OUT $r_gen->{"PREFIX"} . "_com.all.out";
    if($r_gen->{"DC"}){
      print OUT " ",$r_dec->{"COMPRI"};
    }
    print OUT "\n";
  }
  if($r_gen->{"RECEPTOR"}){
    print OUT $r_gen->{"PREFIX"} . "_rec.all.out";
    if($r_gen->{"DC"}){
      if($r_gen->{"COMPLEX"}){
        print OUT " ",$r_dec->{"RECMAP"};
      }
      else{
        print OUT " ",$r_dec->{"RECPRI"};
      }
    }
    print OUT "\n";
  }
  if($r_gen->{"LIGAND"}){
    print OUT $r_gen->{"PREFIX"} . "_lig.all.out";
    if($r_gen->{"DC"}){
      if($r_gen->{"COMPLEX"}){
        print OUT " ",$r_dec->{"LIGMAP"};
      }
      else{
        print OUT " ",$r_dec->{"LIGPRI"};
      }
    }
    print OUT "\n";
  }
  close(OUT);

  my $dctype = 0;
  if($r_gen->{"DC"}){
    if($r_dec->{"DCTYPE"} == 1 || $r_dec->{"DCTYPE"} == 2){
      $dctype = 1;
    }
    elsif($r_dec->{"DCTYPE"} == 3 || $r_dec->{"DCTYPE"} == 4){
      $dctype = 2;
    }
  }

  # Calc statistics
  my $statout = $r_gen->{"PREFIX"} . "_statistics.out";
  my $command;
  if($r_gen->{"COMPLEX"} && $r_gen->{"RECEPTOR"} && $r_gen->{"LIGAND"}){
    &mm_pbsa_stat(1,$dctype,$statin,$statout);
  }
  else{
    &mm_pbsa_stat(0,$dctype,$statin,$statout);
  }
}  
  
sub mm_pbsa_stat(){
###################

  $ | = 1; # Autoflush

  if(scalar(@_) < 4){
    print "\n";
    print "USAGE: &mm_pbsa_stat(<calc_delta? 0..2> <calc_dec? 0..2>\n";
    print "                     <input_file> <output_file>\n";
    print "                     [<snap_min> <snap_max>])";
    die("\nFor details see: $HTMLPATH#mm_pbsa_stat_usage\n");
  }

  # Parameters: $calc_delta,$calc_dec,$input,$output,$snap_min,$snap_max
  my $calc_delta = shift;
  my $calc_dec = shift;
  my $input = shift;
  my $output = shift;
  my ($snap_min, $snap_max);
  if(scalar(@_) == 2){
    $snap_min = shift;
    $snap_max = shift;
  }
  else{
    $snap_min = 1;
    $snap_max = 10e10;
  }

  if($calc_dec == 1){
    $VAR = $VARDEC;
  }
  elsif($calc_dec == 2){
    $VAR = $VARPAI;
  }
  else{
    $VAR = $VARREG;
  }

  # Read files to work with
  my (@files, @decomp);
  &read_input($input,\@files);

  # Reorder files
  &reorder_files(\@files);

  # Read data
  &read_data(\@files,$calc_dec,$snap_min,$snap_max);

  # Print values of parameters
  &print_parameters();

  # Treat special parameters
  &treat_special();

  # Calc missing parameters
  &calc_missing();

  # Calc delta 
  # [i.e. $DATA[3] = $DATA[0] - ($DATA[1]+$DATA[2])
  # and so forth, then(!) averaging]
  if($calc_delta == 1){
    &calc_delta($calc_delta,$calc_dec);
  }

  # Calc average and stddev
  &calc_ave_stdev($calc_delta,$calc_dec);

  # Calc delta 
  # [i.e. $AVG[3] = $AVG[0] - ($AVG[1]+$AVG[2])
  # after(!) averaging]
  if($calc_delta == 2){
    &calc_delta($calc_delta,$calc_dec);
  }

  # Output
  if($calc_dec == 1 || $calc_dec == 2){
    &output_decomp($output,$calc_delta,$calc_dec);
  }
  else{
    &output($output,$calc_delta);
    &output_snap($output.".snap",$calc_delta,$snap_min);
  }
}

sub print_parameters(){
###################

  print "=>> Values of parameters\n";
  print "    TEMP   = $TEMP\n";
  print "    R      = $R\n";
  if(defined $gammaP && defined $betaP){
    print "    gammaP = $gammaP\n";
    print "    betaP  = $betaP\n";
  }
  if(defined $gammaG && defined $betaG){
    print "    gammaG = $gammaG\n";
    print "    betaG  = $betaG\n";
  }
}

sub output(){
#############

  # Parameters: $output,$calc_delta
  my $output = shift;
  my $calc_delta = shift;

  print "=>> Print output to $output\n";
  
  open(OUT,">$output") || die("$output not opened\nFor details see: $HTMLPATH#output_out_not_opened\n");

  # If ! calc_delta, print $AVG/$STD subsequently for every
  #   molecular species,
  # else print for three species in a row and
  #   add delta values below 
  my $no;
  if($calc_delta){
    # scalar(@$AVG) = no. of molecular species + no. of delta data
    #               = no. of molecular species + (no. of molecular species)/3
    $no = scalar(@$AVG)*(1.0 - 1.0/4.0)/3.0;
  }
  else{
    $no = scalar(@$AVG);
  }  

  my $i;
  for($i = 0; $i < $no; $i++){
    if($calc_delta && $i%$no == 0){
      my $j;
      my $header;
      printf OUT "%-10s","#";
      for($j = 0; $j < 3; $j++){
        $header = "COMPLEX" if($j%3 == 0);
        $header = "RECEPTOR" if($j%3 == 1);
        $header = "LIGAND" if($j%3 == 2);
        printf OUT " %15s %7s",$header,"";
      }
      print OUT "\n";
      printf OUT "%-10s","#";
      for($j = 0; $j < 3; $j++){
        printf OUT " %23s","-----------------------";
      }
      print OUT "\n";
      printf OUT "%-10s","#";
      for($j = 0; $j < 3; $j++){
        printf OUT " %12s %10s","MEAN","STD";
      }
      print OUT "\n";
      printf OUT "%-10s","#";
      for($j = 0; $j < 3; $j++){
        printf OUT " %23s","=======================";
      }
      print OUT "\n";
    }
    else{
      my $header;
      printf OUT "%-10s %12s %10s\n","#","MEAN","STD";
      printf OUT "%-10s %23s\n","#","=======================";
    }

    my @delta = ();   
    my $maxout = 0;
    my $isout = 0;
    do{
      my $calc;
      foreach $calc (@$CALC){
        my $var;
        foreach $var (keys %{$VAR->{$calc}}){

          if($VAR->{$calc}->{$var}->[1] > 0){
            # Var. should be printed
            if($VAR->{$calc}->{$var}->[0] > 1){
              die("    Missing values for $calc $var\n    For details see: $HTMLPATH#output_values_missing\n");
            }
            elsif($VAR->{$calc}->{$var}->[0] == 1){
              # No skipping of variable
              $maxout = $VAR->{$calc}->{$var}->[1] > $maxout ? $VAR->{$calc}->{$var}->[1] : $maxout;

              if($isout == $VAR->{$calc}->{$var}->[1]){
                # Print variable now
                printf OUT "%-10s",$var;
                if($calc_delta){
                  # Print data for single species (3 in one row)
                  my $ind = $i * 3;
                  my $j;
                  for($j = 0; $j < 3; $j++){
                    printf OUT " %12.2f %10.2f",$AVG->[$ind + $j]->{$calc}->{$var}->[0],
                                                $STD->[$ind + $j]->{$calc}->{$var}->[0];
                  }
                  
                  # Save data for delta-output
                  $ind = $no * 3 + $i;
                  my $str = sprintf ("%-10s",$var);
                  $str .= sprintf " %12.2f %10.2f\n",$AVG->[$ind]->{$calc}->{$var}->[0],
                                                     $STD->[$ind]->{$calc}->{$var}->[0];
                  push @delta,$str;
                }
                else{
                  # Print data for single species
                  my $ind = $i;
                  printf OUT " %12.2f %10.2f",$AVG->[$ind]->{$calc}->{$var}->[0],
                                              $STD->[$ind]->{$calc}->{$var}->[0];
                }
                print OUT "\n";
              }
            }
          }
        } # end foreach
      } # end foreach
      $isout++;
    }while($isout <= $maxout);
    print OUT "\n";

    if($calc_delta){
      # Print missing delta data
      my $header;
      printf OUT "%-10s %15s %7s\n","#","DELTA","";
      printf OUT "%-10s %23s\n","#","-----------------------";
      printf OUT "%-10s %12s %10s\n","#","MEAN","STD";
      printf OUT "%-10s %23s\n","#","=======================";

      my $str;
      foreach $str (@delta){
        print OUT $str;
      }
    }
  } # end for

  close(OUT);
}

sub output_snap(){
##################

  # Parameters: $output,$calc_delta,$snap_min
  my $output = shift;
  my $calc_delta = shift;
  my $snap_min = shift;

  print "=>> Print output (snap) to $output\n";
  
  open(OUT,">$output") || die("$output not opened\nFor details see: $HTMLPATH#output_snap_not_opened\n");

  my $no;
  if($calc_delta == 1){
    # scalar(@$DATA) = no. of molecular species + no. of delta data
    #                = no. of molecular species + (no. of molecular species)/3
    $no = scalar(@$DATA)*(1.0 - 1.0/4.0)/3.0;
  }
  elsif($calc_delta == 2){
    print "    That operation doesn't make sense with calc_delta == 2\n";
    die("\n    For details see: $HTMLPATH#output_snap_nonsense_operation\n");
  }
  else{
    $no = scalar(@$DATA);
  }  

  my $snap = -1;
  my $i;
  for($i = 0; $i < scalar(@$DATA); $i++){

    # Print header line
    if($i < $no * 3){
      my $header;
      if($calc_delta){
        $header = "COMPLEX" if($i%($no * 3) == 0);
        $header = "RECEPTOR" if($i%($no * 3) == 1);
        $header = "LIGAND" if($i%($no * 3) == 2);
      }
      else{
        $header = "RECEPTOR";
      }
      print OUT "\n=>> $header\n\n";
    }
    else{
      print OUT "\n=>> DELTA\n\n";
    }
    printf OUT "%10s","Snapshot";
    my $maxout = 0;
    my $isout = 0;
    do{
      my $calc;
      foreach $calc (@$CALC){
        my $var;
        foreach $var (keys %{$VAR->{$calc}}){

          if($VAR->{$calc}->{$var}->[1] > 0){
            # Var. should be printed
            if($VAR->{$calc}->{$var}->[0] > 1){
              die("    Missing values for $calc $var\n    For details see: $HTMLPATH#output_snap_missing_value\n");
            }
            elsif($VAR->{$calc}->{$var}->[0] == 1){
              # No skipping of variable
              $maxout = $VAR->{$calc}->{$var}->[1] > $maxout ? $VAR->{$calc}->{$var}->[1] : $maxout;

              if($isout == $VAR->{$calc}->{$var}->[1]){
                # Print variable now
                printf OUT " %10s",$var;
                # Get no. of snapshots
                if($snap < 0){
                  $snap = scalar(@{$DATA->[$i]->{$calc}->{$var}});
                }
                elsif($snap != scalar(@{$DATA->[$i]->{$calc}->{$var}})){
                  die("    Non-equal no. of snapshots: $snap != scalar(@{$DATA->[$i]->{$calc}->{$var}})\n    For details see: $HTMLPATH#output_snap_nonequalnr_snaps\n");
                }
              }
            }
          }
        } # end foreach
      } # end foreach
      $isout++;
    }while($isout <= $maxout);
    print OUT "\n";

    my $j;
    for($j = 0; $j < $snap; $j++){
      printf OUT "%10i",$j+$snap_min; # snap_min is supposed to be >= 1      

      # Print "raw" values   
      my $maxout = 0;
      my $isout = 0;
      do{
        my $calc;
        foreach $calc (@$CALC){
          my $var;
          foreach $var (keys %{$VAR->{$calc}}){
            if($VAR->{$calc}->{$var}->[1] > 0){
              # Var. should be printed
              if($VAR->{$calc}->{$var}->[0] > 1){
                die("    Missing values for $calc $var\n    For details see: $HTMLPATH#output_snap_missing_value\n");
              }
              elsif($VAR->{$calc}->{$var}->[0] == 1){
                # No skipping of variable
                $maxout = $VAR->{$calc}->{$var}->[1] > $maxout ? $VAR->{$calc}->{$var}->[1] : $maxout;

                if($isout == $VAR->{$calc}->{$var}->[1]){
                  # Print variable now
                  printf OUT " %10.2f",$DATA->[$i]->{$calc}->{$var}->[$j];
                }
              }
            }
          } # end foreach
        } # end foreach
        $isout++;
      }while($isout <= $maxout);
      print OUT "\n";
    } # end for
  } # end for

  close(OUT);
}

sub output_decomp(){
####################

  # Parameters: $output,$calc_delta,$calc_dec
  my $output = shift;
  my $calc_delta = shift;
  my $calc_dec = shift;

  print "=>> Print output (decomp) to $output\n";
  
  open(OUT,">$output") || die("$output not opened\nFor details see: $HTMLPATH#output_dec_not_opened\n");

  my $no = scalar(@$AVG);
  my $decno;
  if($calc_delta){
    # scalar(@$AVG) = no. of molecular species + no. of delta data
    #               = no. of molecular species + (no. of molecular species)/3
    $decno = scalar(@$AVG)*(1.0 - 1.0/4.0);
  }
  else{
    $decno = scalar(@$AVG);
  }  

  my $i;
  for($i = 0; $i < $no; $i++){

    # Print header line
    if($i < $decno){
      my $header;
      if($calc_delta){
        $header = "COMPLEX" if($i%$no == 0);
        $header = "RECEPTOR" if($i%$no == 1);
        $header = "LIGAND" if($i%$no == 2);
      }
      else{
        $header = "RECEPTOR";
      }
      print OUT "\n=>> $header\n\n";
    }
    else{
      print OUT "\n=>> DELTA\n\n";
    }
    printf OUT "%10s %10s","Number","Residue";
    if($calc_dec == 2){
      printf OUT " -> %10s","Residue";
    }
    my $maxout = 0;
    my $isout = 0;
    do{
      my $calc;
      foreach $calc (@$CALC){
        my $var;
        foreach $var (keys %{$VAR->{$calc}}){

          if($VAR->{$calc}->{$var}->[1] > 0){
            # Var. should be printed
            if($VAR->{$calc}->{$var}->[0] > 1){
              die("    Missing values for $calc $var\n    For details see: $HTMLPATH#output_dec_missing_value\n");
            }
            elsif($VAR->{$calc}->{$var}->[0] == 1){
              # No skipping of variable
              $maxout = $VAR->{$calc}->{$var}->[1] > $maxout ? $VAR->{$calc}->{$var}->[1] : $maxout;

              if($isout == $VAR->{$calc}->{$var}->[1]){
                # Print variable now
                printf OUT " %10s %10s",$var,"";
              }
            }
          }
        } # end foreach
      } # end foreach
      $isout++;
    }while($isout <= $maxout);
    print OUT "\n";

    # Print "raw" values   
    my $decind;
    if($i >= $decno){
      $decind = ($i - $decno) * 3;
    }
    else{
      $decind = $i;
    }
    my $decsum = 0;
    my $j;
    for($j = 0; $j < scalar(@{$DEC->[$decind]}); $j+=2){
      $decsum += $DEC->[$decind]->[$j+1] - $DEC->[$decind]->[$j] + 1;
    }
    if($calc_dec == 2){
      $decsum *= $decsum;
    }

    for($j = 0; $j < $decsum; $j++){
      my ($mapres1, $mapres2);
      if($calc_dec == 1){
        $mapres1 = &map_res($j,$decind);
      }
      elsif($calc_dec == 2){
        &map_res_pair($j,$decind,\$mapres1,\$mapres2);
      }

      printf OUT "%10i %10i",$j,$mapres1;      
      if($calc_dec == 2){
        printf OUT " -> %10i",$mapres2;
      }
      my $maxout = 0;
      my $isout = 0;
      do{
        my $calc;
        foreach $calc (@$CALC){
          my $var;
          foreach $var (keys %{$VAR->{$calc}}){

            if($VAR->{$calc}->{$var}->[1] > 0){
              # Var. should be printed
              if($VAR->{$calc}->{$var}->[0] > 1){
                die("    Missing values for $calc $var\n    For details see: $HTMLPATH#output_dec_missing_value\n");
              }
              elsif($VAR->{$calc}->{$var}->[0] == 1){
                # No skipping of variable
                $maxout = $VAR->{$calc}->{$var}->[1] > $maxout ? $VAR->{$calc}->{$var}->[1] : $maxout;

                if($isout == $VAR->{$calc}->{$var}->[1]){
                  # Print variable now
                  printf OUT " %10.2f %10.2f",$AVG->[$i]->{$calc}->{$var}->[$j],
					      $STD->[$i]->{$calc}->{$var}->[$j];
                }
              }
            }
          } # end foreach
        } # end foreach
        $isout++;
      }while($isout <= $maxout);
      print OUT "\n";
    } # end for
  } # end for
  close(OUT);
}

sub calc_delta(){
#################

  # Parameters: $calc_delta,$calc_dec
  my $calc_delta = shift;
  my $calc_dec = shift;

  my $ARR;
  if($calc_delta == 1){
    print "=>> Calc delta from raw data\n";
    $ARR = $DATA;
  }
  elsif($calc_delta == 2){
    print "=>> Calc delta from avg data\n";
    $ARR = $AVG;
  }

  my $no = scalar(@$ARR);
  if($no%3 != 0.0){
    die("    Wrong no of molecular species\n    For details see: $HTMLPATH#calc_delta_wrong_no_molecules\n");
  }

  my $i;
  for($i = 0; $i < $no; $i=+3){
    # Take together 3 species to calc. 
    # delta = species0 - species1 - species2
    # and so forth
    push @$ARR,{};
    if($calc_delta == 2){
      push @$STD,{};
    }
    my $lind = scalar(@$ARR)-1; 
    my $j;
    my $decsumcom = 0;
    for($j = 0; $j < scalar(@{$DEC->[$i+0]}); $j+=2){
      $decsumcom += $DEC->[$i+0]->[$j+1] - $DEC->[$i+0]->[$j] + 1;
    }
    my $decsumrec = 0;
    for($j = 0; $j < scalar(@{$DEC->[$i+1]}); $j+=2){
      $decsumrec += $DEC->[$i+1]->[$j+1] - $DEC->[$i+1]->[$j] + 1;
    }
    my $decsumlig = 0;
    for($j = 0; $j < scalar(@{$DEC->[$i+2]}); $j+=2){
      $decsumlig += $DEC->[$i+2]->[$j+1] - $DEC->[$i+2]->[$j] + 1;
    }
    if($calc_dec == 2){
      $decsumcom *= $decsumcom;
      $decsumrec *= $decsumrec;
      $decsumlig *= $decsumlig;
    }

    my $calc;
    foreach $calc (@$CALC){
      my $var;
      foreach $var (keys %{$VAR->{$calc}}){
        if($VAR->{$calc}->{$var}->[0] > 1){
          die("    Missing values for $calc $var\n    For details see: $HTMLPATH#calc_delta_missing_value\n");
        }
        elsif($VAR->{$calc}->{$var}->[0] == 1){
          if(! exists $ARR->[$lind]->{$calc}->{$var}){
            $ARR->[$lind]->{$calc}->{$var} = [];
          }

          my $ndata = scalar(@{$ARR->[$i+0]->{$calc}->{$var}});
          my $j;
          for($j = 0; $j < $ndata; $j++){
            if(! defined $ARR->[$i+0]->{$calc}->{$var}->[$j]){
              my $mes = "    No data for " . "$i" . "+0 $calc $var $j\n    For details see: $HTMLPATH#calc_delta_no_data_plus_zero\n";
              die($mes);
            }

            if($calc_dec == 0){
              # No decomposition
              if(! defined $ARR->[$i+1]->{$calc}->{$var}->[$j]){
                my $mes = "    No data for " . "$i" . "+1 $calc $var $j\n    For details see: $HTMLPATH#calc_delta_no_data_plus_one\n";
                die($mes);
              }
              elsif(! defined $ARR->[$i+2]->{$calc}->{$var}->[$j]){
                my $mes = "    No data for " . "$i" . "+2 $calc $var $j\n    For details see: $HTMLPATH#calc_delta_no_data_plus_two\n";
                die($mes);
              }
              else{
                push @{$ARR->[$lind]->{$calc}->{$var}},$ARR->[$i+0]->{$calc}->{$var}->[$j] - 
                                                       $ARR->[$i+1]->{$calc}->{$var}->[$j] -
                                                       $ARR->[$i+2]->{$calc}->{$var}->[$j];
                if($calc_delta == 2){
                  # Calc std according to error propagation rule
                  push @{$STD->[$lind]->{$calc}->{$var}},sqrt(
                                                           ($STD->[$i+0]->{$calc}->{$var}->[$j])**2 + 
                                                           ($STD->[$i+1]->{$calc}->{$var}->[$j])**2 +
                                                           ($STD->[$i+2]->{$calc}->{$var}->[$j])**2
                                                         );
		}
              }
            }
            else{
              # Consider decomposition - mapping of residues necessary
              my ($maprec, $maplig);
              if($calc_dec == 1){
                $maprec = &map_ind($j%$decsumcom,$i+0,$i+1);
                $maplig = &map_ind($j%$decsumcom,$i+0,$i+2);
              }
              elsif($calc_dec == 2){
                &map_ind_pair($j%$decsumcom,$i,\$maprec,\$maplig);
              }

              my $data = 0.0;
              my $std  = 0.0;
              if(($maprec >= 0 && $maplig >= 0) || 
                 ($maprec <  0 && $maplig <  0)){
                if($calc_dec == 1){
                  print $i,"  ",$calc,"  ",$var,"  ",$ARR->[$i+0]->{$calc}->{$var}->[$j],"\n";
                  die("    Error while mapping for $j: $maprec $maplig\n    For details see: $HTMLPATH#calc_delta_mapping_error\n");
                }
                elsif($calc_dec == 2){
                  #print "    Skipping ",$j%$decsumcom,"\n";
                }
              }
              elsif($maprec >= 0){
                $maprec = int($j/$decsumcom) * $decsumrec + $maprec;
                if(! defined $ARR->[$i+1]->{$calc}->{$var}->[$maprec]){
                  my $mes = "    No data for " . "$i" . "+1 $calc $var $maprec\n    For details see: $HTMLPATH#calc_delta_no_data_plus_one_maprec\n";
                  die($mes);
                }
                else{
                  $data = $ARR->[$i+1]->{$calc}->{$var}->[$maprec];
                  if($calc_delta == 2){
                    $std = $STD->[$i+1]->{$calc}->{$var}->[$maprec];
                  }
                }
              }
              elsif($maplig >= 0){
                $maplig = int($j/$decsumcom) * $decsumlig + $maplig;
                if(! defined $ARR->[$i+2]->{$calc}->{$var}->[$maplig]){
                  my $mes = "    No data for " . "$i" . "+2 $calc $var $maplig\n    For details see: $HTMLPATH#calc_delta_no_data_plus_two_maplig\n";
                  die($mes);
                }
                else{
                  $data = $ARR->[$i+2]->{$calc}->{$var}->[$maplig];
                  if($calc_delta == 2){
                    $std = $STD->[$i+2]->{$calc}->{$var}->[$maplig];
  		  }
                }
              }
              else{
                print "    Actually can't be ... ???\nFor details see: $HTMLPATH#calc_delta_actually_cant_be\n";
              }
              push @{$ARR->[$lind]->{$calc}->{$var}},$ARR->[$i+0]->{$calc}->{$var}->[$j] -
                                                     $data;
              if($calc_delta == 2){
                push @{$STD->[$lind]->{$calc}->{$var}},sqrt(
                                                         ($STD->[$i+0]->{$calc}->{$var}->[$j])**2 +
                                                         $std**2
                                                       );
              }
            } # end else
          } # end for
        }
      } # end foreach
    } # end foreach
  } # end for
}

sub map_res_pair(){
###################

  # Parameters: $j,$decind,\$mapres1,\$mapres2
  my $calcres = shift;
  my $mol = shift;
  my $r_mapres1 = shift;
  my $r_mapres2 = shift;

  # Find residue no. in mol that fits to calcres 
  $$r_mapres1 = -1;
  $$r_mapres2 = -1;
  my $count = 0;
  my $i;
  for($i = 0; $i < scalar(@{$DEC->[$mol]}); $i+=2){
    my $start1 = $DEC->[$mol]->[$i];
    my $end1   = $DEC->[$mol]->[$i+1];
    my $j;
    for($j = $start1; $j <= $end1; $j++){
      my $k;
      for($k = 0; $k < scalar(@{$DEC->[$mol]}); $k+=2){
        my $start2 = $DEC->[$mol]->[$k];
        my $end2   = $DEC->[$mol]->[$k+1];
        my $l;
        for($l = $start2; $l <= $end2; $l++){
          if($count == $calcres){
            $$r_mapres1 = $j;
            $$r_mapres2 = $l;
            last;
          }
          else{
            $count++;
          }
        }
        last if($$r_mapres1 >= 0 && $$r_mapres2 >= 0);
      }
      last if($$r_mapres1 >= 0 && $$r_mapres2 >= 0);
    }
    last if($$r_mapres1 >= 0 && $$r_mapres2 >= 0);
  }
}

sub map_res(){
##############

  # Parameters: $calcres,$mol
  my $calcres = shift;
  my $mol = shift;

  # Find residue no. in mol that fits to calcres 
  my $molres = -1;
  my $count = 0;
  my $i;
  for($i = 0; $i < scalar(@{$DEC->[$mol]}); $i+=2){
    my $start = $DEC->[$mol]->[$i];
    my $end   = $DEC->[$mol]->[$i+1];
    my $j;
    for($j = $start; $j <= $end; $j++){
      if($count == $calcres){
        $molres = $j;
        last;
      }
      else{
        $count++;
      }
    }
    last if($molres >= 0);
  }
  return $molres;
}

sub map_ind_pair(){
###################

  # Parameters: $j%$decsumcom,$mol,\$maprec,\$maplig
  my $calcres = shift;
  my $com = shift;
  my $r_maprec = shift;
  my $r_maplig = shift;

  # Find residue no. in com that fits to calcres 
  my $comres1 = -1;
  my $comres2 = -1;
  my $count = 0;
  my $i;
  for($i = 0; $i < scalar(@{$DEC->[$com]}); $i+=2){
    my $start1 = $DEC->[$com]->[$i];
    my $end1   = $DEC->[$com]->[$i+1];
    my $j;
    for($j = $start1; $j <= $end1; $j++){
      my $k;
      for($k = 0; $k < scalar(@{$DEC->[$com]}); $k+=2){
        my $start2 = $DEC->[$com]->[$k];
        my $end2   = $DEC->[$com]->[$k+1];
        my $l;
        for($l = $start2; $l <= $end2; $l++){
          if($count == $calcres){
            $comres1 = $j;
            $comres2 = $l;
            last;
          }
          else{
            $count++;
          }
        }
        last if($comres1 >= 0 && $comres2 >= 0);
      }
      last if($comres1 >= 0 && $comres2 >= 0);
    }
    last if($comres1 >= 0 && $comres2 >= 0);
  }

  # Find residue no. in mol that fits to com
  my $m;
  for($m = 1; $m <= 2; $m++){
    my $mol = $com + $m;
    my $molind = -1;
    $count = 0;
    for($i = 0; $i < scalar(@{$DEC->[$mol]}); $i+=2){
      my $start1 = $DEC->[$mol]->[$i];
      my $end1   = $DEC->[$mol]->[$i+1];
      my $j;
      for($j = $start1; $j <= $end1; $j++){
        my $k;
        for($k = 0; $k < scalar(@{$DEC->[$mol]}); $k+=2){
          my $start2 = $DEC->[$mol]->[$k];
          my $end2   = $DEC->[$mol]->[$k+1];
          my $l;
          for($l = $start2; $l <= $end2; $l++){
            if($comres1 == $j && $comres2 == $l){
              $molind = $count;
              last;
            }
            else{
              $count++;
            }
          }
          last if($molind >= 0);
        }
        last if($molind >= 0);
      }
      last if($molind >= 0);
    }

    if($m == 1){
      $$r_maprec = $molind;
    }
    elsif($m == 2){
      $$r_maplig = $molind;
    }

  }
}

sub map_ind(){
##############

  # Parameters: $calcres,$com,$mol
  my $calcres = shift;
  my $com = shift;
  my $mol = shift;

  # Find residue no. in com that fits to calcres 
  my $comres = -1;
  my $count = 0;
  my $i;
  for($i = 0; $i < scalar(@{$DEC->[$com]}); $i+=2){
    my $start = $DEC->[$com]->[$i];
    my $end   = $DEC->[$com]->[$i+1];
    my $j;
    for($j = $start; $j <= $end; $j++){
      if($count == $calcres){
        $comres = $j;
        last;
      }
      else{
        $count++;
      }
    }
    last if($comres >= 0);
  }

  # Find residue no. in mol that fits to com
  my $molind = -1;
  $count = 0;
  for($i = 0; $i < scalar(@{$DEC->[$mol]}); $i+=2){
    my $start = $DEC->[$mol]->[$i];
    my $end   = $DEC->[$mol]->[$i+1];
    my $j;
    for($j = $start; $j <= $end; $j++){
      if($comres == $j){
        $molind = $count;
        last;
      }
      else{
        $count++;
      }
    }
    last if($molind >= 0);
  }

  return $molind;
}

sub calc_ave_stdev(){
#####################

  # Parameters: $calc_delta,$calc_dec
  my $calc_delta = shift;
  my $calc_dec = shift;

  print "=>> Calc average and stddev\n";

  my $no = scalar(@$DATA);
  my $decno;
  if($calc_delta == 1){
    # scalar(@$DATA) = no. of molecular species + no. of delta data
    #                = no. of molecular species + (no. of molecular species)/3
    $decno = scalar(@$DATA)*(1.0 - 1.0/4.0);
  }
  else{
    $decno = $no;
  }  

  my $calc;
  foreach $calc (@$CALC){
    my $var;
    foreach $var (keys %{$VAR->{$calc}}){
      if($VAR->{$calc}->{$var}->[0] > 1){
        die("    Missing values for $calc $var\n    For details see: $HTMLPATH#calc_ave_stdev_missing_value\n");
      }
      elsif($VAR->{$calc}->{$var}->[0] == 1){
        my $i;
        for($i = 0; $i < $no; $i++){
           if(! defined $AVG->[$i]){
             push @$AVG,{};
           }
           $AVG->[$#{$AVG}]->{$calc}->{$var} = [];
           if(! defined $STD->[$i]){
             push @$STD,{};
           }
           $STD->[$#{$STD}]->{$calc}->{$var} = [];

           my ($ave, $adev, $sdev, $vari, $skew, $curt);
           my $decind;
           if($i >= $decno){
             $decind = ($i-$decno) * 3;
           }
           else{
             $decind = $i;
           }
           if($calc_dec > 0 && scalar(@{$DEC->[$decind]})){
             # Energy decomposition - collect data from corresponding residues
             # before averaging
             my $decsum = 0;
             my $j;
             for($j = 0; $j < scalar(@{$DEC->[$decind]}); $j+=2){
               $decsum += $DEC->[$decind]->[$j+1] - $DEC->[$decind]->[$j] + 1;
             }
             if($calc_dec == 2){
               $decsum *= $decsum;
             }

             my $snap = scalar(@{$DATA->[$i]->{$calc}->{$var}}) / $decsum;
             for($j = 0; $j < $decsum; $j++){
               my @data = ();
               my $k;
               for($k = 0; $k < $snap; $k++){
                 push @data, $DATA->[$i]->{$calc}->{$var}->[$k * $decsum + $j];
               }
               if(scalar(@data) > 1){
                 &calc_moments_of_distr(\@data,
                                        \$ave,\$adev,\$sdev,\$vari,\$skew,\$curt);
                 push @{$AVG->[$i]->{$calc}->{$var}}, $ave;
                 push @{$STD->[$i]->{$calc}->{$var}}, $sdev;
               }
               else{
                 push @{$AVG->[$i]->{$calc}->{$var}}, $data[0];
                 push @{$STD->[$i]->{$calc}->{$var}}, 0.0;
               }
             }
           }
           elsif($calc_dec == 0){
             if(scalar(@{$DATA->[$i]->{$calc}->{$var}}) > 1){
               &calc_moments_of_distr($DATA->[$i]->{$calc}->{$var},
                                      \$ave,\$adev,\$sdev,\$vari,\$skew,\$curt);
               push @{$AVG->[$i]->{$calc}->{$var}}, $ave;
               push @{$STD->[$i]->{$calc}->{$var}}, $sdev;
             }
             else{
               push @{$AVG->[$i]->{$calc}->{$var}}, $DATA->[$i]->{$calc}->{$var}->[0];
               push @{$STD->[$i]->{$calc}->{$var}}, 0.0;
             } 
             #print $i," ",$calc," ",$var," ",$AVG->[$i]->{$calc}->{$var}," ",
             #                                $STD->[$i]->{$calc}->{$var},"\n";
          }
        }
      }
    } # end foreach
  } # end foreach  
}

sub treat_special(){
####################

  print "=>> Treat special parameters\n";

  my $calc;
  foreach $calc (@$CALC){
    my $var;
    foreach $var (keys %{$VAR->{$calc}}){
      if($VAR->{$calc}->{$var}->[0] == 3){
          
        my $no = scalar(@$DATA);
        my $i;
        for($i = 0; $i < $no; $i++){

          if(scalar(%{$DATA->[$i]})){

	    if(($calc eq "PB" || $calc eq "PBHyb") && $var eq "PBCAL"){
	      # Convert PB from kT to kcal/mol
	      my $kt2kcal = $R * $TEMP / 4.184 / 1000.0;
	      $DATA->[$i]->{$calc}->{$var} = [];
	      my $j;
	      for($j = 0; $j < scalar(@{$DATA->[$i]->{$calc}->{"PB"}}); $j++){
		push @{$DATA->[$i]->{$calc}->{$var}}, $kt2kcal * $DATA->[$i]->{$calc}->{"PB"}->[$j];
	      }
	    }
	    elsif(($calc eq "PB" || $calc eq "PBHyb") && $var eq "PBSUR"){
	      # Calc nonpolar contribution for PB
	      $DATA->[$i]->{$calc}->{$var} = [];
	      my $j;
	      for($j = 0; $j < scalar(@{$DATA->[$i]->{$calc}->{"PBCAV"}}); $j++){
		push @{$DATA->[$i]->{$calc}->{$var}}, $gammaP * $DATA->[$i]->{$calc}->{"PBCAV"}->[$j] + $betaP;
		
		# PBCAV, PBDIS, ELRAELE and EPB are missing when doing iAPBS calculation
		# zero them out so the rest of the script can proceed
		if ((! defined ($DATA->[$i]->{"PB"}->{"PBDIS"}->[$j])) &&
		    (! defined ($DATA->[$i]->{"PB"}->{"ELRAELE"}->[$j])) &&
		    (! defined ($DATA->[$i]->{"PB"}->{"EPB"}->[$j]))){
		  print "treat_special(): No PBDIS, ELRAELE and EPB found, assuming iAPBS calculation.\nFor details see: $HTMLPATH#treat_special_assuming_iapbs\n" if ($j == 0);
		  $DATA->[$i]->{"PB"}->{"PBDIS"}->[$j] = 0.;
		  $DATA->[$i]->{"PB"}->{"ELRAELE"}->[$j] = 0.;
		  $DATA->[$i]->{"PB"}->{"EPB"}->[$j] = 0.;
		}
	      }
	    }
	    elsif($calc eq "GB" && $var =~ /GBSUR/){
	      my $surf;
	      if($var eq "GBSUR"){
		$surf = "SURF";
	      }
	      elsif($var eq "TGBSUR"){
		$surf = "TSURF";
	      }
	      elsif($var eq "SGBSUR"){
		$surf = "SSURF";
	      }
	      elsif($var eq "BGBSUR"){
		$surf = "BSURF";
	      }
	      else{
		die("    No surf specification given in treat_special\n    For details see: $HTMLPATH#treat_special_no_surf_spec\n");
	      }
	      # Calc nonpolar contribution for GB
	      $DATA->[$i]->{$calc}->{$var} = [];
	      my $j;
	      for($j = 0; $j < scalar(@{$DATA->[$i]->{"MS"}->{"$surf"}}); $j++){
		push @{$DATA->[$i]->{$calc}->{$var}}, $gammaG * $DATA->[$i]->{"MS"}->{"$surf"}->[$j] + $betaG;
	      }
	    }
	    elsif($calc eq "NM" && $var =~ /^TS.../ || $var =~ /^[TSB]_TSVIB/){
	      my $entr;
	      if($var eq "TSTRA"){
		$entr = "STRA";
	      }
	      elsif($var eq "TSROT"){
		$entr = "SROT";
	      }
	      elsif($var eq "TSVIB"){
		$entr = "SVIB";
	      }
	      elsif($var eq "TSTOT"){
		$entr = "STOT";
	      }
	      elsif($var eq "T_TSVIB"){
		$entr = "T_SVIB";
	      }
	      elsif($var eq "S_TSVIB"){
		$entr = "S_SVIB";
	      }
	      elsif($var eq "B_TSVIB"){
		$entr = "B_SVIB";
	      }
	      else{
		die("    No entr specification given in treat_special\n    For details see: $HTMLPATH#treat_special_no_entr_spec");
	      }
	      # Convert cal/molK in kcal/mol
	      $DATA->[$i]->{$calc}->{$var} = [];
	      my $j;
	      for($j = 0; $j < scalar(@{$DATA->[$i]->{"NM"}->{"$entr"}}); $j++){
		push @{$DATA->[$i]->{$calc}->{$var}}, 
		  $TEMP * $DATA->[$i]->{"NM"}->{"$entr"}->[$j] / 1000.0;
	      }
	    }
	    else{
	      die("    No rule found for $calc $var\n    For details see: $HTMLPATH#treat_special_no_rules_found\n");
	    }
          } # end if
        } # end for

        $VAR->{$calc}->{$var}->[0] = 1;
      } # end if
    } # end foreach
  } # end foreach
}

sub calc_missing(){
###################

  print "=>> Calc missing parameters\n";

  my $trial = 0;
  my $todo;
  do{ # do loop necessary since hash does not 
      # return keys in sequential order
    $todo = 0;
    my $calc;
    foreach $calc (@$CALC){
      my $var;
      foreach $var (keys %{$VAR->{$calc}}){

        if($VAR->{$calc}->{$var}->[0] == 2){
          print "    Processing $calc $var\n";

          my $no = scalar(@$DATA);
          my $i;
          for($i = 0; $i < $no; $i++){

            if(scalar(%{$DATA->[$i]})){

            $DATA->[$i]->{$calc}->{$var} = [];
            my $expr = $VAR->{$calc}->{$var}->[2];

            my $skip = 0;
            my $count = 0;
            while(length($expr)){
              # Extract calculation information from expression
              $expr =~ s/^(\+)(.+)/$2/; # Extract + or -
              my $sign = $1 eq "+" ? 1.0 : -1.0;
              $expr =~ s/^(\??)(.+)/$2/; # Extract possible ?
              my $quest = $1 eq "?" ? 1 : 0;
              $expr =~ s/^(.+?)_(.+)/$2/; # Extract letters in front of "_" (+? operator is non-greedy)
              my $tocalc = $1;
              $expr =~ s/^(.+?)(\+.+|$)/$2/; # Extract letters behind "_" (|$ takes care of last operand)
              my $tovar = $1;

              print "        Doing $sign $tocalc $tovar\n";

              if(! exists $VAR->{$tocalc}->{$tovar} || $VAR->{$tocalc}->{$tovar}->[0] == 0){
                die("    No entries for $tocalc" . "_" . "$tovar existing\n    For details see: $HTMLPATH#calc_missing_no_entries\n");
              }
              elsif($VAR->{$tocalc}->{$tovar}->[0] > 1){
                # Values not yet available -> consider eventually later again
		if($trial > 0 && $quest){
		  print("    Value for $tocalc" . "_" . "$tovar not mandatory -> Will be omitted\n    For details see: $HTMLPATH#calc_missing_value_not_mand_omitting\n")
		}
		else{
		  print("    No values for $tocalc" . "_" . "$tovar existing -> Skipping\n    For details see: $HTMLPATH#calc_missing_no_value_skipping\n");
		  $skip = 1;
		  $todo = 1;
		  last;
		}
              }

              # Calc missing data
	      my $last;
	      if(defined scalar(@{$DATA->[$i]->{$tocalc}->{$tovar}})){
		$last = scalar(@{$DATA->[$i]->{$tocalc}->{$tovar}});
	      }
	      elsif(defined scalar(@{$DATA->[$i]->{$calc}->{$var}})){
		$last = scalar(@{$DATA->[$i]->{$calc}->{$var}});
	      }
	      else{
		die("    No value for \"last\" available\n    For details see: $HTMLPATH#calc_missing_no_value\n");
	      }

              my $j;
              for($j = 0; $j < $last; $j++){
		my $val;
		if($trial > 0 && $quest && $VAR->{$tocalc}->{$tovar}->[0] > 1){
		  $val = 0.0;
		}
		else{
		  $val = $sign * $DATA->[$i]->{$tocalc}->{$tovar}->[$j];
		}

                if($count == 0){
                  push @{$DATA->[$i]->{$calc}->{$var}}, $val;
                }
                else{
                  $DATA->[$i]->{$calc}->{$var}->[$j] += $val;
                }
              }

              $count++;
            } # end while

            if($skip){
              # Skip also for other species
              last;
            }
            else{
              # Calculation done, data finished
              $VAR->{$calc}->{$var}->[0] = 1;
            }

            } # end if

          } # end for
        } # end if
      } # end foreach
    } # end foreach
    $trial++;
  } while($todo)
}

sub read_data(){
################

  # Parameters: \@files,$calc_dec,$snap_min,$snap_max
  my $r_fil = shift;
  my $calc_dec = shift;
  my $snap_min = shift;
  my $snap_max = shift;

  print "=>> Reading files\n";
  
  my $no = scalar(@$r_fil);
  my $i;
  for($i = 0; $i < $no; $i++){
    push @$DATA,[];
    $DATA->[$i] = {};

    my $decsum = 0;
    my $j;
    for($j = 0; $j < scalar(@{$DEC->[$i]}); $j+=2){
      $decsum += ($DEC->[$i]->[$j+1] - $DEC->[$i]->[$j] + 1);
    }
    if($calc_dec == 2){
      $decsum *= $decsum;
    } 

    my $snap = 0;
    my $check_calc = 0;
    my $check_MS = 0;
    my @tmp_calc = ();
    my @tmp_data = ();
    my $file = $r_fil->[$i];

    if($calc_dec == 0 ||
       ($calc_dec > 0 && scalar(@{$DEC->[$i]}) > 0)){
      print "    Reading $file\n";
      open(IN,"$file") || die("    $file not opened\n    For details see: $HTMLPATH#read_data_not_opened\n");
      my $line;
      my $index = -1;
      my $readflg = 0;
      while(defined($line = <IN>)){
        chomp($line);
        next if($line =~ /^\#/ || $line =~ /^$/);

        if($line =~ /^(\d+)/){
          # Snapshot number at beginning of each entry
          if($index >= 0 && $readflg){
            # Check comleteness of data
            &check_completeness($i,$index,$decsum,$snap_min);
          }
          # Get snapshot number
          $readflg = 0;
          $snap = $1;
          if($snap >= $snap_min && $snap <= $snap_max){
            $readflg = 1;
            $index++;
          }
        }
        elsif($snap == 0 && $i == 0){
          # Get modes of calculation at beginning of first file
	  if($line =~ /^\D+$/){
	    push @$CALC,$line;
	    $check_MS = 1 if $line eq "MS" || $line eq "PB"; # Ensure that non-polar energy is also calculated if MS=0 but PB=1
	  }
	  elsif ($line =~ /^(\D+) +(-?\d+\.\d+)$/){
	    if($check_MS){
	      $gammaP = $2 if($1 eq "PB_SURFTEN");
	      $betaP  = $2 if($1 eq "PB_SURFOFF");
	      $gammaG = $2 if($1 eq "GB_SURFTEN");
	      $betaG  = $2 if($1 eq "GB_SURFOFF");
	    }
	    else{
	      $gammaP = $gammaG = 1.0;
	      $betaP  = $betaG  = 0.0;
	    }
	  }
        }
        elsif($snap == 0 && $i > 0){
          # Get modes of calc and data at beginning of other files
          # for check of consistency
	  if($line =~ /^\D+$/){
	    push @tmp_calc, $line;
	  }
	  elsif ($line =~ /^(\D+) +(\d+\.\d+)$/){
	    push @tmp_data, [ $1, $2 ];
	  }
        }
        elsif($readflg){
          if(! $check_calc && $i > 0){
            # Check consistency of modes of calculation
            &check_calc(\@tmp_calc);
            &check_data(\@tmp_data);
            $check_calc = 1;
          }

          # Get data
          my $calc;
          foreach $calc (@$CALC){
            my $var;
            foreach $var (keys %{$VAR->{$calc}}){
              if($VAR->{$calc}->{$var}->[0] == 1){
                if($line =~ /$VAR->{$calc}->{$var}->[2]/){
                  # Line matches expression -> store data
                  my $value = $1;
                  if(! exists $DATA->[$i]->{$calc}->{$var}){
                    $DATA->[$i]->{$calc}->{$var} = [];
                  }
                  if(defined($DATA->[$i]->{$calc}->{$var}->[$index]) && ($decsum == 0)){
                    ###if($calc eq "PB" && $var eq "PB"){
                    ###  print "    Found 2nd entry for $i PB PB $index -> subtracting\n";
                    ###  $DATA->[$i]->{$calc}->{$var}->[$index] -= $value;
                    ###}
                    ###else{
                      print "    Entry $i $calc $var ",$index," exists -> overwriting\nFor details see: $HTMLPATH#read_data_exists_overwriting\n";
                      $DATA->[$i]->{$calc}->{$var}->[$index] = $value;
                    ###}
                  }
                  else{
                    push @{$DATA->[$i]->{$calc}->{$var}}, $value;
                  }
                }
              }
            }
          }
        } # end else
      } # end while
      close(IN);       

      # Check completeness of data for last snapshot
      &check_completeness($i,$index,$decsum,$snap_min) if($readflg);
    } # end if
  }# end for
}

sub check_completeness(){
#########################

  # Parameters: $fileno,$snap,$decsum,$snap_min
  my $fileno = shift;
  my $index = shift;
  my $decsum = shift;
  my $snap_min = shift;

  #print "    Checking DATA for $snap\n";

  my $calc;
  foreach $calc (@$CALC){
    my $var;
    foreach $var (keys %{$VAR->{$calc}}){
      if($VAR->{$calc}->{$var}->[0] == 1){
        if($decsum == 0){
          # No decomposition
          if(! defined($DATA->[$fileno]->{$calc}->{$var}->[$index])){
            print "    WARNING: Missing $var for $calc in $index -> Taken from ",$index-1,"\n    For details see: $HTMLPATH#check_completeness_missing_var_taken_from\n";
            push @{$DATA->[$fileno]->{$calc}->{$var}}, $DATA->[$fileno]->{$calc}->{$var}->[$index-1];
            #die("    Missing $var for $calc in $index+$snap_min\n");
          }
        }
        else{
          # Decomposition - consider values for single residues
          my $j;
          for($j = $index*$decsum; $j < ($index+1)*$decsum; $j++){
            if(! defined($DATA->[$fileno]->{$calc}->{$var}->[$j])){
              die("    Missing $var for $calc in " . ($index+$snap_min) . " (residue " . ($j+1) .")\n    For details see: $HTMLPATH#check_completeness_missing_var\n");
            }
          }
        }
      }
    }
  }
}

sub check_calc(){
#################

  # Parameters: \@tmp_calc
  my $r_calc = shift;

  print "    Checking consistency of CALC\n";

  my $must;
  foreach $must (sort @$CALC){
    my $found = 0;
    my $is;
    foreach $is (sort @$r_calc){
      if($must eq $is){
        $found = 1;
        last;
      }
    }

    if(!$found){
      die("Calc mode $must not found\nFor details see: $HTMLPATH#check_calc_mode_not_found\n");
    }
  }
}
      
sub check_data(){
#################

  # Parameters: \@tmp_data
  my $r_data = shift;

  print "    Checking consistency of DATA\n";

  my $consistent = 1;
  foreach my $entry (@$r_data){
    my $var = $entry->[0];
    my $val = $entry->[1];

    if($var eq "PB_SURFTEN" && defined $gammaP && $gammaP != $val ||
       $var eq "PB_SURFOFF" && defined $betaP  && $betaP  != $val ||
       $var eq "GB_SURFTEN" && defined $gammaG && $gammaG != $val ||
       $var eq "GB_SURFOFF" && defined $betaG  && $betaG  != $val){
      $consistent = 0;
    }
  }

  if(! $consistent){
    die("SURFTEN / SURFOFF values are not consistent across *.all.out files\nFor details see: $HTMLPATH#check_data_inconsistent_surften_surfoff\n");
  }
}

sub reorder_files(){
####################

  # Parameters: \@files
  my $r_fil = shift;

  print "=>> Reordering files\n";

  # Something to implement for mutational analysis

  print "    Final order:\n";
  my $no = scalar(@$r_fil);
  my $i;
  for($i = 0; $i < $no; $i++){
    print "    ", $i+1 . ". " . $r_fil->[$i] . ":";
    if(scalar(@{$DEC->[$i]})){
      my $j;
      for($j = 0; $j < scalar(@{$DEC->[$i]}); $j+=2){
        print "   ",$DEC->[$i]->[$j]," - ",$DEC->[$i]->[$j+1];
      }
    }
    else{
      print " - ";
    }
    print "\n";
  }
}

sub read_input(){
#################

  # Parameters: $ARGV[0],\@files
  my $input = shift;
  my $r_fil = shift;

  print "=>> Reading input\n";

  open(IN,"$input") || die("    $input not opened\n    For details see: $HTMLPATH#read_input_not_opened\n");
  my $line;
  while(defined($line = <IN>)){
    chomp($line);
    next if($line =~ /^\#/ || $line =~ /^$/);
    my @tmp = split/ +/,$line;
    push @$r_fil,$tmp[0];
    push @$DEC,[];
    my $lind = scalar(@$DEC) - 1;
    my $i;
    for($i = 1; $i < scalar(@tmp); $i++){
      if($tmp[$i] =~ /0-0/ || $tmp[$i] =~ /^0$/){
        # Do nothing
      }
      elsif($tmp[$i] =~ /(\d+)-(\d+)/){
        push @{$DEC->[$lind]},$1;
        push @{$DEC->[$lind]},$2;
      }
      else{
        die("    $tmp[$i] is not conform with \\d+-\\d+ or \"0\"\n    For details see: $HTMLPATH#read_input_not_conform\n");
      }
    }
  }
  close(IN);
}

sub calc_moments_of_distr(){
############################

  # Parameters: \@data,\$ave,\$adev,\$sdev,\$var,\$skew,\$curt
  my $r_data = shift;
  my $r_ave = shift;
  my $r_adev = shift;
  my $r_sdev = shift;
  my $r_var = shift;
  my $r_skew = shift;
  my $r_curt = shift;

  my $n = scalar(@$r_data);

  if($n <= 1){
    print STDERR "\n";
    print STDERR "Error: n must be at least two in calc_moments_of_distr\nFor details see: $HTMLPATH#calc_moments_n_lt_two\n";
    die("\n");
  }

  # Calc average
  my $i;
  my $sum = 0.;
  for($i = 0;$i < $n;$i++){
    $sum += $r_data->[$i];
  }

  $$r_ave = $sum / $n;

  $$r_adev = 0.;
  $$r_var = 0.;
  $$r_skew = 0.;
  $$r_curt = 0.;
  my $ep = 0.;

  for($i = 0;$i < $n;$i++){
    $sum = $r_data->[$i] - $$r_ave;
    $ep += $sum;
    $$r_adev += abs($sum);
    my $p = $sum * $sum;
    $$r_var += $p;
    $p *= $sum;
    $$r_skew += $p;
    $p *= $sum;
    $$r_curt += $p;
  }
  
  $$r_adev /= $n;
  $$r_var = ($$r_var - ($ep * $ep)/$n) / ($n-1);
  $$r_sdev = sqrt($$r_var);

  if($$r_var != 0.){
    $$r_skew /= ($n * $$r_sdev * $$r_sdev * $$r_sdev);
    $$r_curt /= ($n * $$r_var * $$r_var) - 3.;
  }
  else{
    print STDERR "No skew or curtosis when zero variance in moment\nFor details see: $HTMLPATH#calc_moments_skew\n";
    $$r_skew = 0.;
    $$r_curt = 0.;
  }
  
  return;
}

1; # Necessary for package function
