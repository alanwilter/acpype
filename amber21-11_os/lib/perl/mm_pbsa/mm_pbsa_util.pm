#
# Module with utility functions for mm_pbsa
#
# Holger Gohlke: 17.04.2002
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_util;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(init_data final_clean_up);
@EXPORT_OK   = qw();
$VERSION     = 1.00;

########################################################################

use strict;
use mm_pbsa_global qw($HTMLPATH);

########################################################################

sub init_data(){
################

  # Parameters: \%PRO
  my $r_pro = shift;

  print "\n";
  print "=>> Init data\n";

  my $amberexe;
  if(! exists $ENV{"AMBERHOME"} || ! defined $ENV{"AMBERHOME"}){
    print "    Can't find environment variable AMBERHOME\n";
    die("For details see: $HTMLPATH#init_data_cant_find_amberhome\n");
  }
  else{
    $amberexe = $ENV{"AMBERHOME"};
    $amberexe .= "/bin";
    print "    Presuming executables of amber suite to be in $amberexe\n";
  }

  $r_pro->{"SANDER"}   = "$amberexe" . "/" . "sander";
  $r_pro->{"PBSA"}     = "$amberexe" . "/" . "sander"; # Was "pbsa" up to Amber9.
  $r_pro->{"IAPBS"}    = "$amberexe" . "/" . "sander.APBS";
  $r_pro->{"GBNSR6"}   = "$amberexe" . "/" . "gbnsr6";
  $r_pro->{"NMODE"}    = "$amberexe" . "/" . "nmode";
  $r_pro->{"AMBPDB"}   = "$amberexe" . "/" . "ambpdb";
  $r_pro->{"MS"}       = "$amberexe" . "/" . "molsurf";
  $r_pro->{"MAKE_CRD"} = "$amberexe" . "/" . "make_crd_hg";
  $r_pro->{"NABNMODE"} = "$amberexe" . "/" . "mm_pbsa_nabnmode";
}

sub final_clean_up(){
#####################

  # Parameters: \%GEN,\%DEL
  my $r_gen = shift;
  my $r_del = shift;

  if($r_gen->{"PB"}){
    if($r_del->{"PROC"} == 1){
      my $runs = $r_del->{"FOCUS"} + 1;
      my $i;
      for($i = 0; $i < $runs; $i++){
	my $name = "delphi.prm" . "." . $i;
	unlink $name if(-e $name && ! $r_gen->{"VERBOSE"});
      }
    }
    elsif($r_del->{"PROC"} == 2 && ! $r_gen->{"VERBOSE"}){
      unlink "pbsa_com.in" if(-e "pbsa_com.in");
      unlink "pbsa_rec.in" if(-e "pbsa_rec.in");
      unlink "pbsa_lig.in" if(-e "pbsa_lig.in");
    }
    elsif($r_del->{"PROC"} == 3){
      unlink "iapbs.in" if(-e "iapbs.in");
    }
  }

  if( ! $r_gen->{"VERBOSE"}){
    unlink "restrt"        if(-e "restrt");
    unlink "mdinfo"        if(-e "mdinfo");
    unlink "make_crd.in"   if(-e "make_crd.in");
    unlink "sander_com.in" if(-e "sander_com.in");
    unlink "sander_rec.in" if(-e "sander_rec.in");
    unlink "sander_lig.in" if(-e "sander_lig.in");
    unlink "gbnsr6_com.in" if(-e "gbnsr6_com.in");
    unlink "gbnsr6_rec.in" if(-e "gbnsr6_rec.in");
    unlink "gbnsr6_lig.in" if(-e "gbnsr6_lig.in");
    unlink "sanmin_com.in" if(-e "sanmin_com.in");
    unlink "sanmin_rec.in" if(-e "sanmin_rec.in");
    unlink "sanmin_lig.in" if(-e "sanmin_lig.in");
    unlink "nmode_com.in"  if(-e "nmode_com.in");
    unlink "nmode_rec.in"  if(-e "nmode_rec.in");
    unlink "nmode_lig.in"  if(-e "nmode_lig.in");
    unlink "pbsa_com.in"   if(-e "pbsa_com.in");
    unlink "pbsa_rec.in"   if(-e "pbsa_rec.in");
    unlink "pbsa_lig.in"   if(-e "pbsa_lig.in");
  }
}

1; # Necessary for package function
