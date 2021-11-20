#
# Module to calc energy and entropy contributions for mm_pbsa
#
# Holger Gohlke: 18.10.2001
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_calceneent;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(calc_energy_entropy);
@EXPORT_OK   = qw();
$VERSION     = 1.00;

########################################################################

### HG
use strict;
use lib ("$ENV{AMBERHOME}/AmberTools/src/mm_pbsa");
use ForkManager;
use vars qw( $procscnt @procs %procs_data );
use mm_pbsa_global qw($HTMLPATH);

########################################################################

sub calc_energy_entropy(){
##########################

  # Parameters: \%GEN,\%DEC,\%DEL,\%GBO,\%MOL,\%NMO,\%PRO
  my $r_gen = shift;
  my $r_dec = shift;
  my $r_del = shift;
  my $r_gbo = shift;
  my $r_mol = shift;
  my $r_nmo = shift;
  my $r_pro = shift;

  print "\n";
  print "=>> Calculating energy / entropy contributions\n";

  my $suf_list = "";
  $suf_list .= "_com " if($r_gen->{"COMPLEX"});
  $suf_list .= "_rec " if($r_gen->{"RECEPTOR"});
  $suf_list .= "_lig " if($r_gen->{"LIGAND"});

  # Loop over COM/REC/LIG
  #######################
  my $suf;
  foreach $suf (split/ +/,$suf_list){
    # Generic names
    my $ncrd  = $r_gen->{"PATH"} . $r_gen->{"PREFIX"} . $suf . ".crd.";
    my $nbase = $r_gen->{"PREFIX"} . $suf;
    my $npdb  = $nbase . ".pdb.";
    my $npqr  = $nbase . ".pqr.";
    my $nlog  = $nbase . ".delphilog.";
    my $nxyz  = $nbase . ".xyz.";
    my $nmol  = $nbase . ".mslog.";
    my $nout  = $nbase . ".all.out";

    # Result file
    if(! -e $nout){
      # Open new output file
      open(OUT,">$nout") || die("$nout not opened\nFor details see: $HTMLPATH#all_out_not_opened\n");
      print OUT "MM\n"    if($r_gen->{"MM"});
      print OUT "GB\n"    if($r_gen->{"GB"});
      print OUT "PB\n"    if($r_gen->{"PB"} && (! exists $r_del->{"IVCAP"} || exists $r_del->{"IVCAP"} && $r_del->{"IVCAP"} == 0));
      print OUT "PBHyb\n" if($r_gen->{"PB"} && exists $r_del->{"IVCAP"} && $r_del->{"IVCAP"} != 0);
      print OUT "MS\n"    if($r_gen->{"GB"} # MS printed to *all.out file -> 'surface area' output in *all.out file considered for statistical analysis if GB!=0
                          ||($r_gen->{"MS"})); # MS printed to *all.out file -> 'surface area' output in *all.out file considered for statistical analysis if MS!=0 
      print OUT "NM\n"    if($r_gen->{"NM"});

      if($r_gen->{"PB"} && ($r_del->{"PROC"} == 3) && (! $r_gen->{"MS"})){
        print OUT "PB_SURFTEN ", "1.0", "\n";
        print OUT "PB_SURFOFF ", "0.0", "\n";
        print OUT "GB_SURFTEN ", "1.0", "\n";
        print OUT "GB_SURFOFF ", "0.0", "\n";
      }
      else{
        print OUT "PB_SURFTEN ", $r_del->{"SURFTEN"}, "\n" if(exists($r_del->{"SURFTEN"}));
        print OUT "PB_SURFOFF ", $r_del->{"SURFOFF"}, "\n" if(exists($r_del->{"SURFOFF"}));
        print OUT "GB_SURFTEN ", $r_gbo->{"SURFTEN"}, "\n" if(exists($r_gbo->{"SURFTEN"}));
        print OUT "GB_SURFOFF ", $r_gbo->{"SURFOFF"}, "\n" if(exists($r_gbo->{"SURFOFF"}));
      }
    }
    else{
      # Re-open existing file
      open(OUT,">>$nout") || die("$nout not opened\nFor details see: $HTMLPATH#all_out_not_opened\n");
    }

    # Parmtop file
    my $npar;
    $npar = $r_gen->{"COMPT"} if($suf eq "_com");
    $npar = $r_gen->{"RECPT"} if($suf eq "_rec");
    $npar = $r_gen->{"LIGPT"} if($suf eq "_lig");

    # Generic input files
    my $nsan  = "sander" . $suf . ".in";
    my $npbs  = "pbsa" . $suf . ".in";
    my $nscl  = "sanmin" . $suf . ".in";
    my $nnmo  = "nmode" . $suf . ".in";
    my $gbnsr6in = "gbnsr" . $suf . ".in";

    if($r_gen->{"NM"} && $r_nmo->{"PROC"} == 2){
      if(! -e $nscl){
        die("    $nscl not found\n    For details see: $HTMLPATH#sanmin_in_not_found\n");
      }
      if(! -e $nnmo){
        die("    $nnmo not found\n    For details see: $HTMLPATH#nmode_in_not_found\n");
      }
    }
    if(! -e $npar){
      die("$npar not found\nFor details see: $HTMLPATH#prmtop_not_found\n");
    }

    # Prepare parallel execution
    ############################
    my $pnumber = -1;
    $procscnt = 0;
    @procs = ();
    %procs_data = ();
    $procs_data{"START"} = $r_gen->{"START"};
    $procs_data{"OFFSET"} = $r_gen->{"OFFSET"};
    my $pm=new Parallel::ForkManager($r_gen->{"PARALLEL"});
    $pm->run_on_start( 
		      sub { my($pid, $pnumber) = @_;
			    &mark_started_proc($pid, $pnumber);
			  } 
		     );
    $pm->run_on_finish( 
		      sub { my($pid, $exit_code, $pnumber) = @_;
			    &analyze_finished_proc($pid, $exit_code, $pnumber);
			  } 
		     );

    # Loop over all coordinates
    ###########################
    for(my $number =  $r_gen->{"START"};
	   $number <= $r_gen->{"STOP"};
	   $number += $r_gen->{"OFFSET"}){

      # Prepare fork data
      ###################
      $procs_data{$number}->{OUT} = *OUT;
      $procs_data{$number}->{r_gen} = $r_gen;
      $procs_data{$number}->{r_del} = $r_del;
      $procs_data{$number}->{r_dec} = $r_dec;
      $procs_data{$number}->{r_mol} = $r_mol;
      $procs_data{$number}->{r_gbo} = $r_gbo;
      $procs_data{$number}->{suf} = $suf;
      my $sanout1 = "sander" . $suf . "." . $number . ".out";
      $procs_data{$number}->{sanout1} = $sanout1; 
      my $sanout2 = "sanmin" . $suf . "." . $number . ".out";
      $procs_data{$number}->{sanout2} = $sanout2; 
      my $sanres = "sanmin" . $suf . "." . $number . ".restrt";
      $procs_data{$number}->{sanres} = $sanres; 
      my $nmodeout = "nmode" . $suf . "." . $number . ".out";
      $procs_data{$number}->{nmodeout} = $nmodeout; 
      my $outbase = $nlog . $number;
      $procs_data{$number}->{outbase} = $outbase; 
      my $pbsaout = "pbsa" . $suf . "." . $number . ".out";
      $procs_data{$number}->{pbsaout} = $pbsaout; 
      my $iapbsout = "iapbs" . $suf . "." . $number . ".out";
      $procs_data{$number}->{iapbsout} = $iapbsout; 
      my $gbnsr6out = "gbnsr6" . $suf . "." . $number . ".out";
      $procs_data{$number}->{gbnsr6out} = $gbnsr6out;
      my $multruns = 1;
      if((exists($r_del->{"ISTRNG"}) && ($r_del->{"ISTRNG"} > 0.0)) ||
	 (exists($r_del->{"REFE"})   && ($r_del->{"REFE"} > 0))){
	$multruns = 2;
      }
      $procs_data{$number}->{multruns} = $multruns; 
      my $pqr = $npqr . $number;
      $procs_data{$number}->{pqr} = $pqr; 
      my $mol = $nmol . $number;
      $procs_data{$number}->{mol} = $mol; 

      my $ncrdfile = $ncrd . "$number";
      last if(! -e $ncrdfile); # exit the loop (before forking a new thread) if there is no new input file

      # Start fork
      ############
      $pnumber++;
      $pm->start($pnumber) and next;

      print "    Calc contrib for $ncrdfile\n";
      if($number == $r_gen->{"START"} && 
	 ($r_gen->{"MM"} || $r_gen->{"GB"} || $r_gen->{"PB"} || $r_gen->{"MS"} || $r_gen->{"NM"})){
	# Check if atom numbers in prmtop and crdfile agree
	&check_atom_numbers($ncrdfile, $npar);
      }

      # Calc MM and GB and SAS (if implicit)
      ######################################
      if($r_gen->{"MM"} || $r_gen->{"GB"}){
        &calc_MM_GB_SAS(*OUT,$nsan,$sanout1,$ncrdfile,$npar,$r_pro);
      }

      # Calc GBNSR6 (if IGB == 6)
      #####################################
      if(($r_gen->{"GB"}) && ($r_gbo->{"IGB"} == 66)){
        &calc_gbnsr6(*OUT,$gbnsr6in,$gbnsr6out,$ncrdfile,$npar,$r_pro);
      }

      # Generate pdb file for PB, NABnmode, or MS calc
      ################################################
      my $pdb = "";
      if(($r_gen->{"PB"} > 0 && $r_del->{"PROC"} == 1) || 
	 ($r_gen->{"NM"} > 0 && $r_nmo->{"PROC"} == 1) ||
	 ($r_gen->{"MS"} > 0)){
        $pdb = $npdb . $number;
	my $tmp_pdb = "tmp.pdb." . $number;
        &generate_pdb($pdb,$tmp_pdb,$ncrdfile,$npar,$r_pro);
        &center_pdb($pdb,$tmp_pdb);
      }

      # Calc entropy
      ##############
      if($r_gen->{"NM"}){
	if($r_nmo->{"PROC"} == 1){
	  &calc_NABnmode(*OUT,$npar,$pdb,$nmodeout,$r_nmo,$r_pro);
	}
	elsif($r_nmo->{"PROC"} == 2){
	  &calc_NM(*OUT,$nscl,$nnmo,$ncrdfile,$npar,$sanout2,$sanres,$nmodeout,$r_pro);
	}
      }

      # Calc PB
      #########
      if($r_gen->{"PB"}){
	if($r_del->{"PROC"} == 1){
	  my $pdbbase = $npdb . $number;
	  &calc_delphi(*OUT,$pdbbase,$outbase,$multruns,$r_del,$r_pro);
	}
	elsif($r_del->{"PROC"} == 2){
	  &calc_pbsa(*OUT,$npbs,$pbsaout,$ncrdfile,$npar,$r_pro);
	}
	elsif($r_del->{"PROC"} == 3){
	  &calc_iapbs(*OUT,"iapbs.in",$iapbsout,$ncrdfile,$npar,$r_gen,$r_pro,$r_del,$suf);
	}
      }

      # Calc MS
      ###########
      if($r_gen->{"MS"}){
        &generate_pqr($pdb,$pqr,$r_mol,$r_del);
        &calc_ms(*OUT,$pqr,$mol,$r_mol,$r_pro);
      }

      # Final clean up's
      if(-e "${npdb}${number}" && ! $r_gen->{"VERBOSE"}){
        unlink $npdb . $number;
      }

      # Finish fork (ATTENTION: valid exit codes are only 0-255)
      #############
      $pm->finish($pnumber%256);

      # Delete procs_data
      ###################
      delete $procs_data{$number};
    } # end for

    # Wait for all forks to finish and check for remaining analyses
    ###############################################################
    $pm->wait_all_children;
    &check_analyses();

    close(OUT);
  } # end foreach
}

sub mark_started_proc(){
########################

  #Parameters: $pid, $pnumber;

  my $pid = shift;
  my $pnumber = shift;

  #print "    Process $pnumber started with PID $pid\n";
  if(! exists $procs[$pnumber]){
    for(my $i = $#procs + 1; $i < $pnumber; $i++){
      push @procs, 0; # Mark procs as "nothing has been done yet"
    }
  }

  $procs[$pnumber] = 1; # Mark procs as "being started"
}

sub analyze_finished_proc(){
############################

  #Parameters: $pid, $exit_code, $pnumber;

  my $pid = shift;
  my $exit_code = shift;
  my $pnumber = shift;

  # ATTENTION: exit codes can only be in the range of 0-255, thus, do a modulo on the snapshot number
  if($pnumber%256 != $exit_code){
    die("Finished process $pnumber with PID $pid has wrong exit code $exit_code\nFor details see: $HTMLPATH#ana_finished_proc\n");
  }

  #print "    Process $pnumber finished with PID $pid\n";
  $procs[$pnumber] = 2; # Mark procs as "being finished"

  while(defined ($procs[$procscnt]) && $procs[$procscnt] == 2){
    my $cnt = $procscnt * $procs_data{"OFFSET"} + $procs_data{"START"};
    print "    Analyzing $cnt\n";

    # Do analysis here
    my $OUT = $procs_data{$cnt}->{OUT};
    print OUT $cnt,"\n";
    my $r_gen = $procs_data{$cnt}->{r_gen};
    my $r_dec = $procs_data{$cnt}->{r_dec};
    my $r_del = $procs_data{$cnt}->{r_del};
    my $r_mol = $procs_data{$cnt}->{r_mol};
    my $r_gbo = $procs_data{$cnt}->{r_gbo};

    # Ana MM/GB
    my @MMenergies_out;
    if($r_gen->{"MM"} || $r_gen->{"GB"}){
      my $MMenergies_out = &ana_MM_GB_SAS(*OUT,$procs_data{$cnt}->{sanout1},$r_gen,$r_gbo);
      @MMenergies_out = @{$MMenergies_out};
    }

    # Ana GBNSR6
    if(($r_gen->{"MM"} || $r_gen->{"GB"}) && ($r_gbo->{IGB} == 66)){
      &ana_GBNSR6(*OUT, $procs_data{$cnt}->{gbnsr6out},$r_gen,$r_gbo,\@MMenergies_out);
    }

    # Ana NM
    if($r_gen->{"NM"}){
      &ana_NM(*OUT,$procs_data{$cnt}->{suf},$procs_data{$cnt}->{sanres},
	      $procs_data{$cnt}->{sanout2},$procs_data{$cnt}->{nmodeout},$r_gen,$r_dec);
    }

    # Ana PB
    if($r_gen->{"PB"}){
      if($r_del->{"PROC"} == 1){
	&ana_delphi(*OUT,$procs_data{$cnt}->{outbase},$procs_data{$cnt}->{multruns},$r_gen,$r_del);
      }
      elsif($r_del->{"PROC"} == 2){
	&ana_pbsa(*OUT,$procs_data{$cnt}->{pbsaout},$r_gen,$r_del);
      }
      elsif($r_del->{"PROC"} == 3){
	&ana_iapbs(*OUT,$procs_data{$cnt}->{iapbsout},$r_gen);
      }
    }

    # Ana MS
    if($r_gen->{"MS"}){
      &ana_ms(*OUT,$procs_data{$cnt}->{pqr},$procs_data{$cnt}->{mol},$r_mol,$r_gen,$r_del);
    }

    $procs[$procscnt] = 3; # Mark procs as "being analyzed"
    $procscnt++;
  }
}

sub check_analyses(){
#####################

  for(my $i = 0; $i <= $#procs; $i++){
    if($procs[$i] < 3){
      my $pnumber = $i;
      die("Process $pnumber was not successfully analyzed\nFor details see: $HTMLPATH#proc_not_succ_analyzed\n");
    }
  }
}

sub calc_MM_GB_SAS(){
#####################

  # Parameters: *OUT,$nsan,$sanout,$ncrdfile,$npar,$r_pro
  my $OUT = shift;
  my $nsan = shift;
  my $sanout = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $r_pro = shift;

  print "        Calc MM/GB/SAS\n";

  my $command = $r_pro->{"SANDER"} . " -O"   .
                                     " -i "  . $nsan .
                                     " -o "  . $sanout .
                                     " -c "  . $ncrdfile .
                                     " -p "  . $npar;

  system($command) && die("        $command not successful\n        For details see: $HTMLPATH#mmgbsas_command_not_successful\n");

}

sub ana_MM_GB_SAS(){
####################

  # Parameters: *OUT,$sanout,$r_gen,$r_gbo
  my $OUT = shift;
  my $sanout = shift;
  my $r_gen = shift;
  my $r_gbo = shift;
  my @MMenergies_out;

  print "        Ana MM/GB/SAS\n";

  # Write MM results
  open(IN,"$sanout") || die("        $sanout not opened\n        For details see: $HTMLPATH#mmgbsas_sanout_not_opened\n");
  my $finalflg = 0;
  my $line;
  while(defined($line = <IN>)){
    if($r_gen->{"DC"}){
      if($line =~ /^TDC/ || $line =~ /^SDC/ || $line =~ /^BDC/){
        print OUT $line;
      }
    }
    else{
      if($line =~ /FINAL RESULTS/){
        $finalflg = 1;
      }
      elsif($finalflg > 0 &&
            ($line =~ /BOND +=/    ||
             $line =~ /ANGLE +=/   ||
             $line =~ /DIHED +=/   ||
             $line =~ /1-4 VDW +=/ ||
             $line =~ /1-4 EEL +=/ ||
             $line =~ /VDWAALS +=/ ||
             $line =~ /EEL +=/     ||
             ($line =~ /EGB +=/ && $r_gen->{"GB"} > 0) ||
             ($line =~ /ESURF +=/ && $r_gen->{"GB"} > 0 && $r_gen->{"MS"} == 0))
           ){
        $line =~ s/ ESURF +=/surface area =/; # Equal output as for MS
        if(! ($r_gbo->{"IGB"} == 66)){
          print OUT $line;
        }
        elsif($r_gbo->{"IGB"} == 66){
          push(@MMenergies_out,$line);
        }
      }
    }
  }
  close(IN);

  # Clean up MM
  unlink $sanout if( ! $r_gen->{"VERBOSE"});

  # Return array
  return \@MMenergies_out;
}

sub calc_gbnsr6(){
#####################

  # Parameters: *OUT,$gbnsr6in,$gbnsr6out,$npar,$ncrdfile,$r_pro
  my $OUT = shift;
  my $gbnsr6in = shift;
  my $gbnsr6out = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $r_pro = shift;

  print "        Calc GBNSR6\n";

  my $command = $r_pro->{"GBNSR6"} . " -O"   .
                                     " -i "  . $gbnsr6in .
                                     " -o "  . $gbnsr6out .
                                     " -c "  . $ncrdfile .
                                     " -p "  . $npar;

  system($command) && die("        $command not successful\n        For details see: $HTMLPATH#mmgbsas_command_not_successful\n");
}

sub ana_GBNSR6(){
####################

  # Parameters: *OUT,$gbnsr6out,$r_gen
  my $OUT = shift;
  my $gbnsr6out = shift;
  my $r_gen = shift;
  my $r_gbo = shift;
  my $MMenergies_out = shift;
  my @MMenergies_out = @{$MMenergies_out};

  print "        Ana GBNSR6\n";

  # Write GBNSR6 results
  open(IN,"$gbnsr6out") || die("        $gbnsr6out not opened\n");
  my $finalflg = 0;
  my $line;
  my $gb_ele;
  my $gb_ele14;
  my $egb;
  my $sa;

  while(defined($line = <IN>)){
    if($line =~ /FINAL RESULTS/){
      $finalflg = 1;
    }
    elsif($finalflg > 0 &&
      ($line =~ /1-4 EEL +=/    ||
       $line =~ /EELEC +=/     ||
      ($line =~ /EGB +=/ && $r_gen->{"GB"} > 0) ||
      ($line =~ /ESURF +=/ && $r_gen->{"GB"} > 0 && $r_gen->{"MS"} == 0))
      ){
        if($line =~ /EELEC += +([-]?\d+)\.(\d+)/){
           $gb_ele = $1.".".$2;
        }
        if($line =~ /1-4 EEL += +([-]?\d+)\.(\d+)/){
           $gb_ele14 = $1.".".$2;
        }
        if($line =~ /EGB += +([-]?\d+)\.(\d+)/){
           $egb = $1.".".$2;
        }
        if($line =~ /ESURF += +([-]?\d+)\.(\d+)/){
           $sa = $1.".".$2;
        }
     }
  }

  for my $MMenerg (@MMenergies_out){
     if(($MMenerg =~ /VDWAALS +=/)&&($MMenerg =~ /EEL +=/)&&($MMenerg =~ /EGB +=/)){
        my $vdw;
        if($MMenerg =~ /VDWAALS += +([-]?\d+)\.(\d+)/){
           $vdw = $1.".".$2;
        }
        printf OUT ("%1s%9s%14s%2s%9s%14s%2s%12s%14s\n" ," ", "VDWAALS =", $vdw, "  ", "EEL     =",$gb_ele, "  ", "EGB        =", $egb);
     }
     elsif(($MMenerg =~ /1\-4 VDW +=/)&&($MMenerg =~ /1\-4 EEL +=/)&&($MMenerg =~ /RESTRAINT +=/)){
        my $vdw14;
        my $restraint;
        if($MMenerg =~ /1\-4 VDW += +([-]?\d+)\.(\d+)/){
           $vdw14 = $1.".".$2;
        }
        if($MMenerg =~ /RESTRAINT += +([-]?\d+)\.(\d+)/){
           $restraint = $1.".".$2;
        }
        printf OUT ("%1s%9s%14s%2s%9s%14s%2s%12s%14s\n" ," ", "1-4 VDW =", $vdw14, "  ", "1-4 EEL =",$gb_ele14, "  ", "RESTRAINT  =", $restraint);
     }
     elsif($MMenerg =~ /surface area +=/){
        if($r_gbo->{"GBSA"} == 6){
           printf OUT ("%1s%14s%14s\n", " ", "surface area =", $sa);
        }
        else{
           print OUT $MMenerg;
        }
     }
     else{
        print OUT $MMenerg;
     }
  }
  close(IN);

  # Clean up MM
  unlink $gbnsr6out if( ! $r_gen->{"VERBOSE"});
}

sub calc_pbsa(){
################

  # Parameters: $OUT, $pbsain, $pbsaout, $ncrdfile, $npar, $r_pro

  my $OUT = shift;
  my $pbsain = shift;
  my $pbsaout = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $r_pro = shift;

  print "        Calc PBSA\n";

  my $command = $r_pro->{"PBSA"} . " -O"   .
                                   " -i "  . $pbsain .
                                   " -o "  . $pbsaout .
                                   " -c "  . $ncrdfile .
                                   " -p "  . $npar;

  system($command) && die("        $command not successful\n        For details see: $HTMLPATH#pbsa_command_not_successful\n");
}

sub ana_pbsa(){
################

  # Parameters: $OUT, $pbsaout, $r_gen, $r_del

  my $OUT = shift;
  my $pbsaout = shift;
  my $r_gen = shift;
  my $r_del = shift;

  print "        Ana PBSA\n";

  # Write MM results
  open(IN,"$pbsaout") || die("        $pbsaout not opened\n        For details see: $HTMLPATH#pbsaout_not_opened\n");
  my $finalflg = 0;
  my ($lraele, $lravdw);
  my $line;
  while(defined($line = <IN>)){
    if($r_gen->{"DC"}){
      if($line =~ /^TDC/ || $line =~ /^SDC/ || $line =~ /^BDC/){
        print OUT "PB_" . $line;
      }
    }
    elsif($line =~ /FINAL RESULTS/){
      $finalflg = 1;
    }
    elsif($line =~ /Protein-solvent interactions: +(-?\d+\.\d+) +(-?\d+\.\d+)/){
      $lraele = $1;
      $lravdw = $2;
    }
    elsif($finalflg > 0){
      # Convert reaction field energy to kT, because back-conversion to kcal/mol is done in
      #   mm_pbsa_statistics, subroutine treat_special (to be compatible with delphi calculation)
      my $kcal2kt = 1000.0 * 4.184 / (300.0 * 8.314);

      if ($line =~ /VDWAALS.+EEL.+EPB += +(-?\d+\.\d+)/){
	# print OUT $line;
	my $rfe;
	if(exists $r_del->{"IVCAP"} && ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5)){
	  # Use hybrid solvent approach to estimate solvation free energy
	  $rfe = (0.5 * $lraele + $1) * $kcal2kt;
	  print OUT "EDISPER = ${lravdw}\n";
	  print OUT "ELRAELE = ", 0.5 * ${lraele} * $kcal2kt, "\n";
	  print OUT "EPB = ", $1 * $kcal2kt, "\n";
	}
	else{
	  # Use pure implicit solvent to estimate solvation free energy
	  $rfe = $1 * $kcal2kt;
	}
	printf OUT "%s %15.6f\n", "corrected reaction field energy:", $rfe;
      }
      elsif($line =~ /ECAVITY +=/ && $r_gen->{"PB"} > 0 && $r_gen->{"MS"} == 0){
	print OUT $line;
      }
    }
  }
  close(IN);

  # Clean up MM
  unlink $pbsaout if( ! $r_gen->{"VERBOSE"});
}

sub calc_iapbs(){
################

  # Parameters: $OUT, $pbsain, $pbsaout, $ncrdfile, $npar, $r_gen, $r_pro, $r_del, $suf

  my $OUT = shift;
  my $pbsain = shift;
  my $pbsaout = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $r_gen = shift;
  my $r_pro = shift;
  my $r_del = shift;
  my $suf = shift;

  print "        Calc Amber/iAPBS\n";

  my $command = $r_pro->{"IAPBS"} . " -O"   .
                                    " -i "  . $pbsain .
                                    " -o "  . $pbsaout .
                                    " -c "  . $ncrdfile .
                                    " -p "  . $npar;


  my $pqrfname = "";
  if($r_del->{"RADIOPT"} == 2){
      $pqrfname  = "com.pqr" if ($suf eq "_com");
      $pqrfname  = "rec.pqr" if ($suf eq "_rec");
      $pqrfname  = "lig.pqr" if ($suf eq "_lig");
      symlink $pqrfname, " ./pqr" || die ("Linking pqr file $pqrfname failed!\nFor details see: $HTMLPATH#iapbs_linking_pqr_failed\n");
  }
  system($command) && die("        $command not successful\n        For details see: $HTMLPATH#iapbs_command_not_successful\n");
  if($r_del->{"RADIOPT"} == 2){
      unlink "./pqr";
  }
}

sub ana_iapbs(){
################

  # Parameters: $OUT, $pbsaout, $r_gen

  my $OUT = shift;
  my $pbsaout = shift;
  my $r_gen = shift;

  print "        Ana Amber/iAPBS\n";

  # Write iapbs results
  open(IN,"$pbsaout") || die("        $pbsaout not opened\n        For details see: $HTMLPATH#ana_iapbs_not_opened");
  my $finalflg = 0;
  my $line;
  while(defined($line = <IN>)){
    if($line =~ /FINAL RESULTS/){
      $finalflg = 1;
    }
    elsif($finalflg > 0){
      if ($line =~ /VDWAALS.+EEL.+EPB. += +(-?\d+\.\d+)/){
        # Convert to kT, because back-conversion to kcal/mol is done in
        #   mm_pbsa_statistics, subroutine treat_special
	#  (to be compatible with delphi calculation)
        my $kcal2kt = 1000.0 * 4.184 / (300.0 * 8.314);
        my $rfe = $1 * $kcal2kt;
        printf OUT "%s %15.6f\n", "corrected reaction field energy:", $rfe;
      }
      elsif($line =~ /ENPOLAR +=/ && $r_gen->{"PB"} > 0 && $r_gen->{"MS"} == 0){
        $line =~ s/ENPOLAR +=/ECAVITY =/; # Equal output as for PB
        print OUT $line;
      }
    }
  }
  close(IN);

  # Clean up MM
  unlink $pbsaout if( ! $r_gen->{"VERBOSE"});
}

sub calc_delphi(){
##################

  # Parameters: $OUT,$pdbbase,$outbase,$multruns,\%DEL,\%PRO
  my $OUT = shift;
  my $pdbbase = shift;
  my $outbase = shift;
  my $multruns = shift;
  my $r_del = shift;
  my $r_pro = shift;

  print "        Calc delphi\n";

  # Delphi calculation

  # Input PDB file
  unlink "fort.13" if(-l "fort.13");
  unlink "fort.13" if(-e "fort.13");
  symlink $pdbbase, "fort.13"; # Because of this no parallel comp. is possible with delphi

  # Do delphi (with focussing)
  my $focruns = $r_del->{"FOCUS"} + 1;

  my $j;
  for($j = 0; $j < $multruns; $j++){
    my $i;
    my $delphiout = "";
    for($i = 0; $i < $focruns; $i++){
      print "            Run $i $j\n";
      unlink "ARCDAT" if(-e "ARCDAT");
      unlink $delphiout if(-e $delphiout); # Remove delphiout from last run
      $delphiout = $outbase . "." . $i . "." . $j;
      my $delphiprm = "delphi.prm." . $i . "." . $j;
      my $command = $r_pro->{"DELPHI"} . " " . $delphiprm . " > " . $delphiout;
      system($command) && die("        $command not running properly\n        For details see: $HTMLPATH#delphi_not_running_properly\n");
    }
  }

  return;
}

sub ana_delphi(){
#################

  # Parameters: $OUT,$outbase,$multruns,\%GEN,\%DEL
  my $OUT = shift;
  my $outbase = shift;
  my $multruns = shift;
  my $r_gen = shift;
  my $r_del = shift;

  print "        Ana delphi\n";

  # Write PB results
  my $solv = 0.0;
  my $j;
  for($j = 0; $j < $multruns; $j++){
    my $delphiout = $outbase . "." . $r_del->{"FOCUS"} . "." . $j;
    open(IN,"$delphiout") || die("        $delphiout not opened\n        For details see: $HTMLPATH#ana_delphi_not_opened\n");
    my $line;
    while(defined($line = <IN>)){
      if($r_del->{"REFE"} == 0){
	if($j == 0 && $line =~ /corrected reaction field energy: +(-?\d+\.\d+)/){
	  $solv += $1;
	}
	if($multruns == 2 && $line =~ /total energy \(including grid energy\): +(-?\d+\.\d+)/){
	  # solv = solv(ISTRNG = 0.0) + total(ISTRNG > 0.0) - total(ISTRNG = 0.0)
	  #        delphi run 1         delphi run 2          delphi run 1
	  if($j == 0){
	    $solv -= $1;
	  }
	  else{
	    $solv += $1;
	  }
	}
      }
      else{ 
	# REFE > 0 does not work with ISTRNG > 0.0 yet
	if($j == 0 && $line =~ /corrected reaction field energy: +(-?\d+\.\d+)/){
	  $solv -= $1; # 1.0/INDI calc
	}
	elsif($line =~ /corrected reaction field energy: +(-?\d+\.\d+)/){
	  $solv += $1; # EXDI/INDI calc
	}
      }
    }
    close(IN);
    unlink $delphiout if( ! $r_gen->{"VERBOSE"});
  }
  printf OUT "%s %15.6f\n", "corrected reaction field energy:", $solv; 

  # Clean up PB
  unlink "fort.13";
  unlink "fort.14";
  unlink "ARCDAT";

  return;
}

sub calc_NABnmode(){
####################

  # Parameters: *OUT,$npar,$pdb,$nmodeout,\%NMO,\%PRO

  my $OUT = shift;
  my $npar = shift;
  my $npdb = shift;
  my $nmodeout = shift;
  my $r_nmo = shift;
  my $r_pro = shift;

  print "        Minimize structure and calc entropy\n";
  my $mmopt = $r_nmo->{"MMOPT"};
  my $command = $r_pro->{"NABNMODE"} . " $npdb $npar" 
                                     . ' "' . $mmopt . '" ' 
                                     . $r_nmo->{"MAXCYC"} . " " 
                                     . $r_nmo->{"DRMS"} 
                                     . " > " . $nmodeout . " 2>&1";

  system($command) && die("        $command not running properly\n        For details see: $HTMLPATH#calc_NABnmode_not_running_properly\n");
}

sub calc_NM(){
##############

  # Parameters: *OUT,$nscl,$nnmo,$ncrdfile,$npar,$sanout,$sanres,$nmodeout,\%PRO
  my $OUT = shift;
  my $nscl = shift;
  my $nnmo = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $sanout = shift;
  my $sanres = shift;
  my $nmodeout = shift;
  my $r_pro = shift;

  print "        Minimize structure\n";
  my $command = $r_pro->{"SANDER"} . " -O"  .
                                     " -i " . $nscl .
                                     " -o " . $sanout .
                                     " -c " . $ncrdfile .
                                     " -p " . $npar .
                                     " -r " . $sanres; 
  system($command) && die("        $command not running properly\n        For details see: $HTMLPATH#calc_NM_min_not_running_properly\n");

  print "        Calc entropy\n";
  $command = $r_pro->{"NMODE"} . " -O"  .
                                 " -i " . $nnmo .
                                 " -o " . $nmodeout .
                                 " -c " . $sanres .
                                 " -p " . $npar;
  system($command) && die("        $command not running properly\n        For details see: $HTMLPATH#calc_NM_entropy_not_running_properly\n");

}

sub ana_NM(){
##############

  # Parameters: *OUT,$suf,$sanres,$sanout,$nmodeout,\%GEN,\%DEC
  my $OUT = shift;
  my $suf = shift;
  my $sanres = shift;
  my $sanout = shift;
  my $nmodeout = shift;
  my $r_gen = shift;
  my $r_dec = shift;

  print "        Ana entropy\n";

  open(IN,"$nmodeout") || die("         $nmodeout not opened\n         For details see: $HTMLPATH#ana_NM_entropy_not_opened\n");
  my $line;
  while(defined($line = <IN>)){
    if($r_gen->{"DC"}){
      if($line =~ /^TDC/ || $line =~ /^SDC/ || $line =~ /^BDC/){
        my $res = substr($line,3,7);
        my $tmp;
        if($suf eq "_com"){$tmp = $r_dec->{"COMPRI"}}
        if($suf eq "_rec"){$tmp = $r_dec->{"RECPRI"}}
        if($suf eq "_lig"){$tmp = $r_dec->{"LIGPRI"}}
        my @ranges = split(' ', $tmp);
        foreach my $range (@ranges){
          my @limit = split('-', $range);
          if($res >= $limit[0] and $res <= $limit[1]){
            print OUT $line
          }
        }
      }
    }
    elsif($line =~ /(Total|translational|rotational|vibrational):? +\S+ +\S+ +(\d+\.\d+)/){
      print OUT "$1 $2\n";
    }
  }
  close(IN);

  # Clean up
  if( ! $r_gen->{"VERBOSE"}){
    unlink $sanout;
    unlink $sanres;
    unlink $nmodeout;
  }
}

sub calc_ms(){
##############

  # Parameters: $out,$pqr,$mol,$r_mol,$r_pro
  my $OUT = shift;
  my $pqr = shift;
  my $mol = shift;
  my $r_mol = shift;
  my $r_pro = shift;

  print "        Calc MS\n";
  my $command = $r_pro->{"MS"} . " $pqr " . $r_mol->{"PROBE"} . " > " . $mol;
  system($command);

}

sub ana_ms(){
#############

  # Parameters: $out,$pqr,$mol,$r_mol,$r_gen,$r_del
  my $OUT = shift;
  my $pqr = shift;
  my $mol = shift;
  my $r_mol = shift;
  my $r_gen = shift;
  my $r_del = shift;

  print "        Ana MS\n";

  open(IN,"$mol") || die("         $mol not opened\n         For details see: $HTMLPATH#ana_ms_not_opened\n");
  my $line;
  while(defined($line = <IN>)){
    if($line =~ /surface area =/){
      print OUT $line;
      if ($r_gen->{"PB"} > 0 ){
        $line =~ s/surface area =/ECAVITY =/; # Equal output as for PB
        print OUT $line;
	if( ! exists $r_del->{"IVCAP"} || ($r_del->{"IVCAP"} != 1 && $r_del->{"IVCAP"} != 5)){
	  print OUT "EDISPER = 0.0000\n";
	}
      }
    }
  }
  close(IN);

  if( ! $r_gen->{"VERBOSE"}){
    unlink $pqr;
    unlink $mol;
  }
}

sub generate_pqr(){
###################

  # This implementation uses bondi radii

  # Parameters: $pdb,$pqr,$r_mol,$r_gen
  my $pdb = shift;
  my $pqr = shift;
  my $r_mol = shift;
  my $r_del = shift;

  # Bondi radii + 1.4A and probe radius of 0.0A yields SAS
  # Bondi radii + 0.0A and probe radius of 1.4A yields molecular surface
  # Bondi radii + 0.0A and probe radius of 0.0A yields vdW surface
  my %exp_rad;

  if(exists $r_del->{"IVCAP"} && ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5)){
    # Prepare for calc of molecular surface
    %exp_rad = (
		   "N"    => 1.550,
		   "H"    => 1.200,
		   "C"    => 1.700,
		   "O"    => 1.500,
		   "P"    => 1.800,
		   "S"    => 1.800,
		   "FE"   => 1.300,
		   "Na+"  => 1.200,
		   "Cl-"  => 1.700,
		   "MG"   => 1.180,
		   "F"    => 1.470, # A. Bondi, J. Phys. Chem. 1964,68, 441-451.
		   "Br"   => 1.850, # A. Bondi, J. Phys. Chem. 1964,68, 441-451.  
		   );
   }
   else{
     # Prepare for calc of molecular surface
     %exp_rad = (
		    "N"    => 1.550 + 1.400,
		    "H"    => 1.200 + 1.400,
		    "C"    => 1.700 + 1.400,
		    "O"    => 1.500 + 1.400,
		    "P"    => 1.800 + 1.400,
		    "S"    => 1.800 + 1.400,
		    "FE"   => 1.300 + 1.400,
		    "Na+"  => 1.200 + 1.400,
		    "Cl-"  => 1.700 + 1.400,
		    "MG"   => 1.180 + 1.400,
		    "F"    => 1.470 + 1.400, # A. Bondi, J. Phys. Chem. 1964,68, 441-451.
		    "Br"   => 1.850 + 1.400, # A. Bondi, J. Phys. Chem. 1964,68, 441-451.
		   );
   }

  print "        Generate PQR\n";
  make_pqr_file($pdb,$pqr,$r_del,\%exp_rad);
}

sub make_pqr_file(){
####################

  # Parameters: $pdb,$pqr,$r_del,\%exp_rad
  my $pdb = shift;
  my $pqr = shift;
  my $r_del = shift;
  my $r_exp_rad = shift;

  open(PDB,"$pdb")  || die("        $pdb not opened\n        For details see: $HTMLPATH#make_pqr_pdb_not_opened\n");
  open(PQR,">$pqr") || die("        $pqr not opened\n        For details see: $HTMLPATH#make_pqr_pqr_not_opened\n");
  my $line;
  while(defined($line = <PDB>)){
    chomp($line);
    if($line =~ /(ATOM|HETATM)/){
      my(        $card, $atnum, $atm1, $atm2, $alt, $resname, $resno, $x, $y, $z) =
	 unpack("a6     a5 x    a      a3     a     a3 x2     a4 x4   a8  a8  a8", $line);

      if(exists $r_del->{"IVCAP"} && ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5)&&
	 # Don't consider solvent molecules or counter ions for surface calculation
	 ($resname eq "WAT" ||
	  $resname eq "DC4" ||
	  $resname eq "PL3" ||
	  $resname eq "SPC" ||
	  $resname eq "SPF" ||
	  $resname eq "T4E" ||
	  $resname eq "TP3" ||
	  $resname eq "TP4" ||
	  $resname eq "TP5" ||
	  $resname eq "CIO" ||
	  $resname eq "Cl-" ||
	  $resname eq "Cs+" ||
	  $resname =~ "IB"  ||
	  $resname =~ "K+"  ||
	  $resname eq "Li+" ||
	  $resname eq "MG2" ||
	  $resname eq "Na+" ||
	  $resname eq "Rb+")){
	next;
      }

      my $atm = $atm1 . $atm2;
      my $atmtmp = $atm;  # Remove leading and trailing white spaces because patterns below
      $atmtmp =~ s/\s//g; #   are of the type /^\S+$/

      my $rad;
      if(exists $r_exp_rad->{"$atmtmp"}){
	# First check full atom name
	$rad = $r_exp_rad->{"$atmtmp"};
      }
      else{
	# Then check character only
	$atmtmp =~ m/([A-Za-z]*)/;
	$atmtmp = $1;
	if(! exists $r_exp_rad->{"$atmtmp"}){
		$atmtmp = substr($atmtmp, 0, 1);
	}
	 
	if(exists $r_exp_rad->{"$atmtmp"}){
	  $rad = $r_exp_rad->{"$atmtmp"};
	}
	else{
	  print "        No radius found for $atm $atnum in residue $resname $resno\n";
	  die("        For details see: $HTMLPATH#make_pqr_no_radius_found\n");
	}
      }

      $line = sprintf("%-6s%5d %1s%-3s%1s%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
		      $card,
		           $atnum,
		               $atm1, $atm2,
                                      $alt,
                                         $resname,
                                               $resno,$x,  $y,  $z,  "00.000", $rad);

    } # endif
    print PQR $line,"\n";
  } # endwhile
  close(PQR);
  close(PDB);
}

sub generate_pdb(){
###################

  # Parameters: $pdb,$ncrdfile,$npar,$r_pro
  my $pdb = shift;
  my $tmp_pdb = shift;
  my $ncrdfile = shift;
  my $npar = shift;
  my $r_pro = shift;

  print "        Generate PDB\n";

  unlink "$tmp_pdb" if(-e "$tmp_pdb");
  my $command = $r_pro->{"AMBPDB"} . " -p " . $npar .
                                     " -c " . $ncrdfile . " > " . $tmp_pdb . " 2> /dev/null";
  system("$command") && die("        $command not successful\n        For details see: $HTMLPATH#gen_pdb_command_not_successful\n");
  
  # ----- The following implements the former "check_names" program
  #       However, ambpdb produces "better" pdb files than anal
  #         -> this part should be redundant.
  open(IN,"$tmp_pdb") || die("        $tmp_pdb not opened\n        For details see: $HTMLPATH#gen_pdb_not_opened\n");
  open(PDB,">$pdb")  || die("        $pdb not opened\n        For details see: $HTMLPATH#gen_pdb_not_opened");
  my $line;
  while(defined($line = <IN>)){
    #              ATOM            1000       1   HB   rest_of_line
    #              $1              $2     $3  $4  $5   $6
    chomp($line);
    if($line =~ /(ATOM +|HETATM +)(\d+)( +)( |\d)(\S+)( +.+$)/){
       $line = $1 . $2 . $3 . $5 . $4 . $6;
    }
    print PDB $line,"\n";
  }
  close(PDB);
  close(IN);
  # -----

  unlink "$tmp_pdb" if(-e "$tmp_pdb");
}

sub center_pdb(){
#################

  # Parameters: $pdb
  my $pdb = shift;
  my $tmp_pdb = shift;

  print "        Center PDB\n";

  unlink "$tmp_pdb" if(-e "$tmp_pdb");

  open(IN,"$pdb") || die("        $pdb not opened\n        For details see: $HTMLPATH#center_pdb_not_opened");
  my @cont = <IN>;
  close(IN);

  my $line;
  my ($x,$y,$z);
  my $xcent = 0.0;
  my $ycent = 0.0;
  my $zcent = 0.0;
  my $count = 0;
  foreach $line (@cont){
    chomp($line);
                #ATOM      1 N    MET     1      55.891  31.265  43.091
    if($line =~ /ATOM(.{26})( *-?\d+\.\d+)( *-?\d+\.\d+)( *-?\d+\.\d+)/){
      $x = substr $line, 30, 8;
      $y = substr $line, 38, 8;
      $z = substr $line, 46, 8;
      $xcent += $x;
      $ycent += $y;
      $zcent += $z;
      $count++;
    }
    elsif($line =~ /ATOM/){
      die("        Unrecognized ATOM format in pdb: $line\n        For details see: $HTMLPATH#center_pdb_unrec_atom_format");
    }
  }

  $xcent /= $count;
  $ycent /= $count;
  $zcent /= $count;

  open(TMP,">$tmp_pdb") || die("        $tmp_pdb not opened\n        For details see: $HTMLPATH#center_pdb_not_opened");
  foreach $line (@cont){
    chomp($line);
                #ATOM      1 N    MET     1      55.891  31.265  43.091
    if($line =~ /ATOM(.{26})( *-?\d+\.\d+)( *-?\d+\.\d+)( *-?\d+\.\d+)/){
      $x = substr $line, 30, 8;
      $y = substr $line, 38, 8;
      $z = substr $line, 46, 8;
      $x -= $xcent;
      $y -= $ycent;
      $z -= $zcent;
      printf TMP "ATOM%s%8.2f%8.2f%8.2f\n", $1, $x, $y, $z;
    }
    else{
      print TMP $line,"\n";
    }
  }
  close(TMP);

  rename "$tmp_pdb", $pdb;
}

sub check_atom_numbers(){
#########################

  # Parameters: $ncrdfile, $npar
  my $ncrdfile = shift;
  my $npar = shift;

  print "        Checking atom numbers\n";

  open(IN,"$ncrdfile") || die("        $ncrdfile not opened\n        For details see: $HTMLPATH#check_atom_nr_ncrd_not_opened");
  my $line;
  my $cnt = 0;
  while($cnt < 2 && defined($line = <IN>)){
    $cnt++;
  }
  chomp($line);
  $line =~ s/^\s+//;
  my $crdno = $line;
  close(IN);

  open(IN,"$npar") || die("        $npar not opened\n        For details see: $HTMLPATH#check_atom_nr_npar_not_opened");
  my $prmtopno = 0;
  $line = <IN>;
  if($line =~ /^%/){
    # New prmtop format detected
    while(defined($line = <IN>)){
      last if($line =~ /%FLAG POINTERS/);
    }
    $line = <IN>; # Read format line
    $line = <IN>;
    chomp($line);
    $line =~ s/^\s+//;
    my @linecont = split(/ +/, $line);
    $prmtopno = $linecont[0];
  }
  else{
    # Old parmtop format detected
    $line = <IN>;
    chomp($line);
    $line =~ s/^\s+//;
    my @linecont = split(/ +/, $line);
    $prmtopno = $linecont[0];
  }
  close(IN);

  if($crdno != $prmtopno){
    die("        Prmtop $npar and coord.file $ncrdfile have different atom numbers\n        For details see: $HTMLPATH#check_atom_nr_prmtop_and_crd_differ");
  }
}

1; # Necessary for package function
