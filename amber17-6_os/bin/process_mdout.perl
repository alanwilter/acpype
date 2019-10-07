#!/usr/bin/perl


$#ARGV >= 0 || die "usage: $0 mdout_filename [more_mdout_filenames...] \n";


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
@sortedkeys = sort by_number keys(%TIME);
@sortedavgkeys = sort by_number keys(%AVG_TIME);

foreach $i ( TEMP, TSOLUTE, TSOLVENT, PRES, EKCMT, ETOT, EKTOT, EPTOT, DENSITY, VOLUME, ESCF ) {
    print "Outputing summary.$i\n";
    open(OUTPUT, "> summary.$i");
    %outarray = eval "\%$i";
    foreach $j ( @sortedkeys ) {
	print OUTPUT "$j  ", $outarray{$j}, "\n";
    }
    close (OUTPUT);

    print "Outputing summary_avg.$i\n";
    open(OUTPUT, "> summary_avg.$i");
    %outarray = eval "\%AVG_$i";
    foreach $j ( @sortedavgkeys ) {
	print OUTPUT "$j  ", $outarray{$j}, "\n";
    }
    close (OUTPUT);

    print "Outputing summary_rms.$i\n";
    open(OUTPUT, "> summary_rms.$i");
    %outarray = eval "\%RMS_$i";
    foreach $j ( @sortedavgkeys ) {
	print OUTPUT "$j  ", $outarray{$j}, "\n";
    }
    close (OUTPUT);


}


sub by_number {
    if ($a < $b) {
	-1;
    } elsif ($a == $b) {
	0;
    } elsif ($a > $b) {
	1;
    }
}

sub process_input {

    $status = 0;
    $debug = 0;
    while ( <INPUT> ) {
	$string = $_;
	
	print $_ if ( ! /NB-upda/ && $debug );

	if (/A V E R A G E S/) {
	    $averages = 1;
	    ($averages_over) = /.*O V E R.*(\d*).*S T E P S/;
	}

	$rms = 1 if (/R M S/);

	if (/NSTEP/) {
	    ($time, $temp, $pres) =
		/NSTEP =.*TIME.* =(.*\d*\.\d*).*TEMP.* =(.*\d*\.\d*).*PRESS = (.*\d*\.\d*)/;
	    if ( $debug ) {
		print $_;
		print "time is $time, temp is $temp, pres is $pres\n";
	    }
	    $_ = <INPUT>;

	    if (/Etot/) {
		($etot, $ektot, $eptot) =
		    /Etot.*=(.*\d*\.\d*).*EKtot.*=(.*\d*\.\d*).*EPtot.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "Etot is $etot, ektot is $ektot, eptot is $eptot\n";
		}
		$_ = <INPUT>;
	    }
	    if (/BOND.*ANGLE.*DIHED/) {
		($bond, $angle, $dihedral) =
		    /BOND.*=(.*\d*\.\d*).*ANGLE.*=(.*\d*\.\d*).*DIHED.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "bond is $bond, angle is $angle, dihedral is $dihedral\n";
		}
		$_ = <INPUT>;
	    }
	    if (/1-4 NB/) {
		($nb14, $eel14, $nb) =
		    /1-4 NB.*=(.*\d*\.\d*).*1-4 EEL.*=(.*\d*\.\d*).*VDWAALS.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "nb14 is $nb14, eel14 is $eel14, vdwaals is $nb\n";
		}
		$_ = <INPUT>;
	    }
	    if (/EELEC/) {
		($eel, $ehbond, $constraint) =
		    /EELEC.*=(.*\d*\.\d*).*EHBOND.*=(.*\d*\.\d*).*CONSTRAINT.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "eel is $eel, ehbond is $ehbond, constraint is $constraint\n";
		}
		$_ = <INPUT>;
#
#               check to see if EAMBER is in the mdout file (present when
#               NTR=1)
#
	        if ( /EAMBER/ ) {
		    $_ = <INPUT>;
		}
	    }
	    if (/EKCMT/) {
		($ekcmt, $virial, $volume) =
		    /EKCMT.*=(.*\d*\.\d*).*VIRIAL.*=(.*\d*\.\d*).*VOLUME.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "Ekcmt is $ekcmt, virial is $virial, volume is $volume\n";
		}
		$_ = <INPUT>;
	    }
	    if (/T_SOLUTE/) {
		($tsolute, $tsolvent) =
		    /T_SOLUTE =(.*\d*\.\d*).*T_SOLVENT =(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "Temp solute is $tsolute, temp solvent is $tsolvent\n";
		}
		$_ = <INPUT>;
	    }

	    if (/Density/) {
		($density) = /.*Density.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "Density is $density\n";
		}
		$_ = <INPUT>;
	    }

	    if (/Etot/) {
		($etot, $ektot, $eptot) =
		    /Etot.*=(.*\d*\.\d*).*EKtot.*=(.*\d*\.\d*).*EPtot.*=(.*\d*\.\d*)/;
		if ( $debug ) {
		    print $_;
		    print "Etot is $etot, ektot is $ektot, eptot is $eptot\n";
		}
		$_ = <INPUT>;
	    }
            if (/ESCF/) {
                ($escf) = 
                    /.*ESCF.*=(.*\d*\.\d*)/;
                if ( $debug ) {
                    print $_;
                    print "ESCF is $escf\n";
                }
                $_ = <INPUT>;
            }

#       update arrays

	    if ( $averages == 1 ) {
		$AVG_TIME{$time}       = $time;
		$AVG_TEMP{$time}       = $temp;
		$AVG_PRES{$time}       = $pres;
		$AVG_ETOT{$time}       = $etot;
		$AVG_EKTOT{$time}      = $ektot;
		$AVG_EPTOT{$time}      = $eptot;
		$AVG_BOND{$time}       = $bond;
		$AVG_ANGLE{$time}      = $angle;
		$AVG_DIHEDRAL{$time}   = $dihedral;
		$AVG_NB14{$time}       = $nb14;
		$AVG_EEL14{$time}      = $eel14;
		$AVG_NB{$time}         = $nb;
		$AVG_EEL{$time}        = $eel;
		$AVG_EHBOND{$time}     = $ehbond;
		$AVG_CONSTRAINT{$time} = $constraint;
		$AVG_EKCMT{$time}      = $ekcmt;
		$AVG_VIRIAL{$time}     = $virial;
		$AVG_VOLUME{$time}     = $volume;
		$AVG_TSOLUTE{$time}    = $tsolute;
		$AVG_TSOLVENT{$time}   = $tsolvent;
		$AVG_DENSITY{$time}    = $density;
		$AVG_ESCF{$time}       = $escf;
		$averages = 0;
	    } elsif ( $rms == 1 ) {
		$RMS_TIME{$time}       = $time;
		$RMS_TEMP{$time}       = $temp;
		$RMS_PRES{$time}       = $pres;
		$RMS_ETOT{$time}       = $etot;
		$RMS_EKTOT{$time}      = $ektot;
		$RMS_EPTOT{$time}      = $eptot;
		$RMS_BOND{$time}       = $bond;
		$RMS_ANGLE{$time}      = $angle;
		$RMS_DIHEDRAL{$time}   = $dihedral;
		$RMS_NB14{$time}       = $nb14;
		$RMS_EEL14{$time}      = $eel14;
		$RMS_NB{$time}         = $nb;
		$RMS_EEL{$time}        = $eel;
		$RMS_EHBOND{$time}     = $ehbond;
		$RMS_CONSTRAINT{$time} = $constraint;
		$RMS_EKCMT{$time}      = $ekcmt;
		$RMS_VIRIAL{$time}     = $virial;
		$RMS_VOLUME{$time}     = $volume;
		$RMS_TSOLUTE{$time}    = $tsolute;
		$RMS_TSOLVENT{$time}   = $tsolvent;
		$RMS_DENSITY{$time}    = $density;
		$RMS_ESCF{$time}       = $escf;
		
		$rms = 0;
	    } else {
		$TIME{$time}       = $time;
		$TEMP{$time}       = $temp;
		$PRES{$time}       = $pres;
		$ETOT{$time}       = $etot;
		$EKTOT{$time}      = $ektot;
		$EPTOT{$time}      = $eptot;
		$BOND{$time}       = $bond;
		$ANGLE{$time}      = $angle;
		$DIHEDRAL{$time}   = $dihedral;
		$NB14{$time}       = $nb14;
		$EEL14{$time}      = $eel14;
		$NB{$time}         = $nb;
		$EEL{$time}        = $eel;
		$EHBOND{$time}     = $ehbond;
		$CONSTRAINT{$time} = $constraint;
		$EKCMT{$time}      = $ekcmt;
		$VIRIAL{$time}     = $virial;
		$VOLUME{$time}     = $volume;	
		$TSOLUTE{$time}    = $tsolute;
		$TSOLVENT{$time}   = $tsolvent;
		$DENSITY{$time}    = $density;
		$ESCF{$time}       = $escf;
	    }

	}
    }
}



