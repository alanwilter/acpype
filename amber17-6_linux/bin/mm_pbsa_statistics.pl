#!/usr/bin/perl -w
#
# Script to calc statistics on mm_pbsa output
#
# Holger Gohlke: 18.10.2001
#
# Last update: 14.02.2014
#

use strict;
use lib ("$ENV{AMBERHOME}/AmberTools/src/mm_pbsa");
use mm_pbsa_statistics;
use mm_pbsa_global qw($HTMLPATH);

########################################################################################
########################################################################################

MAIN:{
######

  $ | = 1; # Autoflush

  if(scalar(@ARGV) < 4){
    print "\n";
    print "USAGE: mm_pbsa_statistics.pl <calc delta ? 0..2> <calc decomp ? 0..2>\n";
    print "                             <input file> <output file>\n";
    print "                             [<snap_min> <snap_max>]\n";
    die("\nFor details see: $HTMLPATH#mm_pbsa_stat_usage\n");
  }
  my $calc_delta = $ARGV[0];
  my $calc_dec = $ARGV[1];
  my $input = $ARGV[2];
  my $output = $ARGV[3];
  my ($snap_min, $snap_max);
  if(scalar(@ARGV) > 4){
    $snap_min = $ARGV[4];
    $snap_max = $ARGV[5];
  }
  else{
    $snap_min = 1;
    $snap_max = 10e10;
  }

  &mm_pbsa_stat($calc_delta, $calc_dec, $input, $output, $snap_min, $snap_max);
}

