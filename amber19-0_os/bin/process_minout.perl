#!/usr/bin/perl -w

use strict;
my $ARGV;
my $filein;
my $i;
my $checkfile;
my @sortedkeys;
my %NSTEP;
my %outarray;
my $j;
my $outarray;
my @f;
my $status;
my $debug;
my $string;
my $maxcyc;
my $maxcyc10;
my $nstep;
my $energy;
my $rms;
my $gmax;
my $name;
my $number;
my $bond;
my $angle;
my $dihedral;
my $vdwaals;
my $eel;
my $hbond;
my $vdw14;
my $eel14;
my $restraint;
my $NSTEP;
my $ENERGY;
my $RMS;
my $GMAX;
my $NAME;
my $NUMBER;
my $BOND;
my $ANGLE;
my $DIHEDRAL;
my $VDWAALS;
my $EEL;
my $HBOND;
my $VDW14;
my $EEL14;
my $RESTRAINT;
my %ENERGY;
my %RMS;
my %GMAX;
my %NAME;
my %NUMBER;
my %BOND;
my %ANGLE;
my %DIHEDRAL;
my %VDWAALS;
my %EEL;
my %HBOND;
my %VDW14;
my %EEL14;
my %RESTRAINT;

$#ARGV >= 0 || die "usage: $0 mdout_filename [more_mdout_filenames...] \n";

# Flush the buffer
$| = 1;

foreach $i ( 0..$#ARGV ) {
    $filein = $ARGV[$i];
    $checkfile = $filein;
    $checkfile =~ s/\.Z//;
    if ( $filein ne $checkfile ) {
	open(INPUT, "zcat $filein |") || 
	    die "Cannot open compressed $filein -- $!\n";
    } else {
	open(INPUT, $filein) || die "Cannot open $filein -- $!\n";
    }
    print "Processing sander output file ($filein)...\n";
    &process_input;
    close(INPUT);
}

print "Starting output...\n";
@sortedkeys = sort {$a <=> $b} keys(%NSTEP);


foreach $i ( "NSTEP", "ENERGY", "RMS", "GMAX", "NAME", "NUMBER", "BOND", "ANGLE", "DIHEDRAL", "VDWAALS", "EEL", "HBOND", "VDW14", "EEL14", "RESTRAINT" ) {
    print "Outputing summary.$i\n";
    open(OUTPUT, "> summary.$i");
    %outarray = eval "\%$i";
    foreach $j ( @sortedkeys ) {
	print OUTPUT "$j  ", $outarray{$j}, "\n";
    }
    close (OUTPUT);
}

sub process_input {
	
	my @f;
	
    $status = 0;
    $debug = 0;
    while ( <INPUT> ) {
	$string = $_;
	
	next if $string =~ /^\s*$/;
	
	# Initialise the output variables
	undef $bond;
	undef $angle;
	undef $dihedral;
	undef $vdwaals;
	undef $eel;
	undef $hbond;
	undef $vdw14;
	undef $eel14;
	undef $restraint;
	
	if (/maxcyc.*ncyc/) {
		@f = split /\s+/, $_;
		$maxcyc = $f[3];
		$maxcyc =~ s/,$//;
		$maxcyc10 = $maxcyc/10;
	}
	
	if (/NSTEP/) {
		$_ = <INPUT>;
	   	@f = split /\s+/, $_;
		$nstep = $f[1];
		$energy = $f[2];
		$rms = $f[3];
		$gmax = $f[4];
		$name = $f[5];
		$number = $f[6];
		if ($debug) {
			foreach (@f){
				print "\"$_\" ";
			}
			print "\n";
		} else {
			# Only print every 10% steps
			if ($nstep % $maxcyc10 == 0) {
				print "Processing step $nstep of a possible $maxcyc...\n";
			}
		}
	    if ( $debug ) {
		print $_;
		print "nstep is $nstep, energy is $energy, rms is $rms, gmax is $gmax, name is $name, number is $number\n";
	    }
	    $_ = <INPUT>;
	    $_ = <INPUT>;

	    if (/BOND.*ANGLE.*DIHED/) {
		($bond, $angle, $dihedral) =
		    /BOND.*=(.*\d*\.?\d*).*ANGLE.*=(.*\d*\.?\d*).*DIHED.*=(.*\d*\.?\d*)/;
		$bond =~ s/^\s*\*+\s*$/99999999.9999/;
		$angle =~ s/^\s*\*+\s*$/99999999.9999/;
		$dihedral =~ s/^\s*\*+\s*$/99999999.9999/;
		if ( $debug ) {
		    print $_;
		    print "bond is $bond, angle is $angle, dihedral is $dihedral\n";
		}
		$_ = <INPUT>;
	    }
	    if (/VDWAALS/) {
		($vdwaals, $eel, $hbond) =
		    /VDWAALS.*=(.*\d*\.?\d*).*EEL.*=(.*\d*\.?\d*).*HBOND.*=(.*\d*\.?\d*)/;
		$vdwaals =~ s/^\s*\*+\s*$/99999999.9999/;
		$eel =~ s/^\s*\*+\s*$/99999999.9999/;
		$hbond =~ s/^\s*\*+\s*$/99999999.9999/;
		if ( $debug ) {
		    print $_;
		    print "vdwaals is $vdwaals, eel is $eel, hbond is $hbond\n";
		}
		$_ = <INPUT>;
	    }
	    if (/1-4 VDW/) {
		($vdw14, $eel14, $restraint) =
		    /1-4 VDW.*=(.*\d*\.?\d*).*1-4 EEL.*=(.*\d*\.?\d*).*RESTRAINT.*=(.*\d*\.?\d*)/;
		$vdw14 =~ s/^\s*\*+\s*$/99999999.9999/;
		$eel14 =~ s/^\s*\*+\s*$/99999999.9999/;
		$restraint =~ s/^\s*\*+\s*$/99999999.9999/;
		if ( $debug ) {
		    print $_;
		    print "vdw14 is $vdw14, eel14 is $eel14, restraint is $restraint\n";
		}
		$_ = <INPUT>;
#
#               check to see if EAMBER is in the mdout file (present when
#               NTR=1) - removed from process_minout.perl (in process_mdout)
#
#       update arrays

		$NSTEP{$nstep}     = $nstep;
		$ENERGY{$nstep}    = $energy;
		$RMS{$nstep}       = $rms;
		$GMAX{$nstep}      = $gmax;
		$NAME{$nstep}      = $name;
		$NUMBER{$nstep}    = $number;
		$BOND{$nstep}      = $bond;
		$ANGLE{$nstep}     = $angle;
		$DIHEDRAL{$nstep}  = $dihedral;
		$VDWAALS{$nstep}   = $vdwaals;
		$EEL{$nstep}       = $eel;
		$HBOND{$nstep}     = $hbond;
		$VDW14{$nstep}     = $vdw14;
		$EEL14{$nstep}     = $eel14;
		$RESTRAINT{$nstep} = $restraint;
	    }

	}
    }
}



