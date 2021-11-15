#
# Module with functions to read mm_pbsa input and check for sanity
#
# Holger Gohlke: 17.04.2002
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_readinput;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(read_input check_sanity);
@EXPORT_OK   = qw();
$VERSION     = 1.00;

########################################################################

use strict;
use lib ("$ENV{AMBERHOME}/AmberTools/src/mm_pbsa");
use mm_pbsa_global qw(
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

########################################################################

sub read_input(){
#################

  # Parameters: $file,\%GEN,\%DEC,\%SAN,\%DEL,\%GBO,\%MOL,\%NMO,\%MAK,\%PRO,\$TRA
  my $file  = shift;
  my $r_gen = shift;
  my $r_dec = shift;
  my $r_san = shift;
  my $r_del = shift;
  my $r_gbo = shift;
  my $r_mol = shift;
  my $r_nmo = shift;
  my $r_mak = shift;
  my $r_pro = shift;
  my $r_tra = shift;

  print "\n";
  print "=>> Reading input parameters\n";

  open(IN,"$file") || die("$file not opened\nFor details see: $HTMLPATH#read_input_not_opened\n");
  my $line;
  my $atflg = "";
  my ($key,$value);

  while(defined($line = <IN>)){
    next if($line =~ /^\#/ || $line =~ /^$/);
    chomp($line);
    $line =~ s/^\s+//; # remove leading whitespace

    # Get @-entry
    if($line !~ /^@/ && $atflg eq ""){
      die("    At-entry must be first\n    For details see: $HTMLPATH#read_input_at_entry_first\n");
    }
    elsif($line =~ /^@(\S+)/){
      $atflg = $1;
    }
    else{
      my @tmp = (split/\s+/,$line);
      $key = $tmp[0];
      $value = join(" ",(@tmp)[1..$#tmp]);    
      print "    Found $key => $value\n";

      # GENERAL parameters
      if($atflg eq "GENERAL"){
        if(exists($r_gen->{$key})){
          die("    $key already exists in GEN\n    For details see: $HTMLPATH#read_input_key_existing_gen\n");
        }
        else{
          $r_gen->{$key} = $value;
        }
      }
      # DECOMP parameters
      elsif($atflg eq "DECOMP"){
        if(exists($r_dec->{$key})){
          die("    $key already exists in DEC\n    For details see: $HTMLPATH#read_input_key_existing_dec\n");
        }
        else{
          $r_dec->{$key} = $value;
        }
      }  
      # PB parameters
      elsif($atflg eq "PB"){
        if($key eq "PERFIL" || $key eq "SCALE"){
	  if(! exists($r_del->{$key})){
	    $r_del->{$key} = [];
	  }
	  push @{$r_del->{$key}},$value;
        }
        elsif(exists($r_del->{$key})){
          die("    $key already exists in DEL\n    For details see: $HTMLPATH#read_input_key_existing_del\n");
        }
        else{
          $r_del->{$key} = $value;
        }
      }
      # MM parameters
      elsif($atflg eq "MM"){
        if(exists($r_san->{$key})){
          die("    $key already exists in SAN\n    For details see: $HTMLPATH#read_input_key_existing_san\n");
        }
        else{
          $r_san->{$key} = $value;
        }
      }
      # GB parameters
      elsif($atflg eq "GB"){
        if(exists($r_gbo->{$key})){
          die("    $key already exists in GBO\n    For details see: $HTMLPATH#read_input_key_existing_gbo\n");
        }
        else{
          $r_gbo->{$key} = $value;
        }
      }
      # MS parameters
      elsif($atflg eq "MS"){
        if(exists($r_mol->{$key})){
          die("    $key already exists in MOL\n    For details see: $HTMLPATH#read_input_key_existing_mol\n");
        }
        else{
          $r_mol->{$key} = $value;
        }
      }
      # NM parameters
      elsif($atflg eq "NM"){
        if(exists($r_nmo->{$key})){
          die("    $key already exists in NMO\n    For details see: $HTMLPATH#read_input_key_existing_nmo\n");
        }
        else{
          $r_nmo->{$key} = $value;
        }
      }
      # TRAJECTORY
      elsif($atflg eq "TRAJECTORY"){
        if($key ne "TRAJECTORY"){
          die("    Unknown key $key where TRAJECTORY expected\n    For details see: $HTMLPATH#read_input_unk_key_trajectory\n");
        }
        else{
          if($$r_tra eq ""){
            $$r_tra .= $value;
          }
          else{
            $$r_tra .= " " . $value;
          }
        }
      }
      # MAKECRD parameters
      elsif($atflg eq "MAKECRD"){
        # LSTART/LSTOP
        if($key eq "LSTART" || $key eq "LSTOP"){
          if(! exists($r_mak->{"NUMBER_LIG_GROUPS"})){
            die("    NUMBER_LIG_GROUPS not yet found\n    For details see: $HTMLPATH#read_input_num_lig_grp_not_found\n");
          }
          else{
            if($key eq "LSTART"){
              if(! exists($r_mak->{"LSTART"})){
                $r_mak->{"LSTART"} = [];
              }
              push @{$r_mak->{"LSTART"}},$value;
            }
            elsif($key eq "LSTOP"){
              if(! exists($r_mak->{"LSTOP"})){
                $r_mak->{"LSTOP"} = [];
              }
              push @{$r_mak->{"LSTOP"}},$value;
            }
          }
        }
        # RSTART/RSTOP
        elsif($key eq "RSTART" || $key eq "RSTOP"){
          if(! exists($r_mak->{"NUMBER_REC_GROUPS"})){
            die("    NUMBER_REC_GROUPS not yet found\n    For details see: $HTMLPATH#read_input_num_rec_grp_not_found\n");
          }
          else{
            if($key eq "RSTART"){
              if(! exists($r_mak->{"RSTART"})){
                $r_mak->{"RSTART"} = [];
              }
              push @{$r_mak->{"RSTART"}},$value;
            }
            elsif($key eq "RSTOP"){
              if(! exists($r_mak->{"RSTOP"})){
                $r_mak->{"RSTOP"} = [];
              }
              push @{$r_mak->{"RSTOP"}},$value;
            }
          }
        }
        # Save rest
        else{
          if(exists($r_mak->{$key})){
           die("    $key already exists in MAK\n    For details see: $HTMLPATH#read_input_key_existing_mak\n");
          }
          else{
            $r_mak->{$key} = $value;
          }
        }
      }
      # ALASCAN parameters -> also stored in $r_mak !
      elsif($atflg eq "ALASCAN"){
        # MUTANTS
        if($key eq "MUTANT_ATOM1" || $key eq "MUTANT_ATOM2" ||
              $key eq "MUTANT_KEEP"  || $key eq "MUTANT_REFERENCE"){
          if(! exists($r_mak->{"NUMBER_MUTANT_GROUPS"})){
            die("    NUMBER_MUTANT_GROUPS not yet found\n    For details see: $HTMLPATH#read_input_num_mut_grp_not_found\n");
          }
          else{
            if($key eq "MUTANT_ATOM1"){
              if(! exists($r_mak->{"MUTANT_ATOM1"})){
                $r_mak->{"MUTANT_ATOM1"} = [];
              }
              push @{$r_mak->{"MUTANT_ATOM1"}},$value;
            }
            elsif($key eq "MUTANT_ATOM2"){
              if(! exists($r_mak->{"MUTANT_ATOM2"})){
                $r_mak->{"MUTANT_ATOM2"} = [];
              }
              push @{$r_mak->{"MUTANT_ATOM2"}},$value;
            }
            elsif($key eq "MUTANT_KEEP"){
              if(! exists($r_mak->{"MUTANT_KEEP"})){
                $r_mak->{"MUTANT_KEEP"} = [];
              }
              push @{$r_mak->{"MUTANT_KEEP"}},$value;
            }
            elsif($key eq "MUTANT_REFERENCE"){
              if(! exists($r_mak->{"MUTANT_REFERENCE"})){
                $r_mak->{"MUTANT_REFERENCE"} = [];
              }
              push @{$r_mak->{"MUTANT_REFERENCE"}},$value;
            }
          }
        }
        # Save rest
        else{
          if(exists($r_mak->{$key})){
           die("    $key already exists in MAK\n    For details see: $HTMLPATH#read_input_key_existing_mak\n");
          }
          else{
            $r_mak->{$key} = $value;
          }
        }
      }
      # PROGRAMS
      elsif($atflg eq "PROGRAMS"){
        $r_pro->{$key} = $value; # Overwriting of entries allowed here since 
                                 #   default is set in &init_data()
      }
      # Wrong input
      else{
        die("    Found unknown atflg $atflg\n   For details see: $HTMLPATH#read_input_found_unk_atflg\n");
      }
    } # end else
  } #end while
  close(IN);
  return;
}

sub check_sanity(){
###################

  # Parameters: \%GEN,\%DEC,\%SAN,\%DEL,\%GBO,\%MOL,\%NMO,\%MAK,\%PRO,\$TRA
  my $r_gen = shift;
  my $r_dec = shift;
  my $r_san = shift;
  my $r_del = shift;
  my $r_gbo = shift;
  my $r_mol = shift;
  my $r_nmo = shift;
  my $r_mak = shift;
  my $r_pro = shift;
  my $r_tra = shift;

  print "\n";
  print "=>> Checking sanity\n";

  # Test GENERAL parameters
  print "    Checking GENERAL\n";
  &test_par(\@GEN_PAR,[$r_gen]);

  if(! ($r_gen->{"COMPLEX"} > 0 && $r_gen->{"RECEPTOR"} > 0 && $r_gen->{"LIGAND"} > 0) &&
     (($r_gen->{"COMPLEX"} > 0 && $r_gen->{"RECEPTOR"} > 0 && $r_gen->{"LIGAND"} ==0) ||
      ($r_gen->{"COMPLEX"} > 0 && $r_gen->{"RECEPTOR"} ==0 && $r_gen->{"LIGAND"} > 0) ||
      ($r_gen->{"COMPLEX"} ==0 && $r_gen->{"RECEPTOR"} > 0 && $r_gen->{"LIGAND"} > 0))){
    die("    Something wrong with COMPLEX/RECEPTOR/LIGAND\n    For details see: $HTMLPATH#check_sanity_something_wrong_comreclig\n");
  }

  if(! exists $r_gen->{"START"} ||
     ! defined $r_gen->{"START"}){
    $r_gen->{"START"} = 1;
    print "    Setting START to default 1\n";
  }
  if(! exists $r_gen->{"STOP"} ||
     ! defined $r_gen->{"STOP"}){
    $r_gen->{"STOP"} = 10e10;
    print "    Setting STOP to default 10e10\n";
  }
  if(! exists $r_gen->{"OFFSET"} ||
     ! defined $r_gen->{"OFFSET"}){
    $r_gen->{"OFFSET"} = 1;
    print "    Setting OFFSET to default 1\n";
  }
  if($r_gen->{"START"} > $r_gen->{"STOP"}){
    die("    START must be smaller than STOP\n    For details see: $HTMLPATH#check_sanity_start_gt_stop\n");
  }
  if($r_gen->{"OFFSET"} < 1){
    die("    OFFSET must be larger than 0\n    For details see: $HTMLPATH#check_sanity_offset_lt_zero\n");
  }

  if(! exists $r_gen->{"VERBOSE"} ||
     ! defined $r_gen->{"VERBOSE"}){
    $r_gen->{"VERBOSE"} = 0;
    print "    Setting VERBOSE to default 0\n";
  }

  if(exists $r_gen->{"PARALLEL"} &&
     defined $r_gen->{"PARALLEL"} &&
    $r_gen->{"PARALLEL"} < 0){
    die "    PARALLEL must be larger than or equal to 0\n    For details see: $HTMLPATH#check_sanity_parallel_gt_zero\n";
  }
  elsif(! exists $r_gen->{"PARALLEL"} ||
	! defined $r_gen->{"PARALLEL"}){
    $r_gen->{"PARALLEL"} = 0;
    print "    Setting PARALLEL to default 0\n";
  }

  if($r_gen->{"GC"} == 0 && 
     $r_gen->{"AS"} == 0 &&
     $r_gen->{"DC"} == 0 &&
     $r_gen->{"MM"} == 0 && 
     $r_gen->{"NM"} == 0 && 
     $r_gen->{"PB"} == 0 && 
     $r_gen->{"GB"} == 0 && 
     $r_gen->{"MS"} == 0){
    die("    No route for calculation\n    For details see: $HTMLPATH#check_sanity_no_route\n");
  }

  if($r_gen->{"GB"} > 0 && $r_gbo->{"GBSA"} == 0){
    print "    If GB calculation is requested surface area must be calulated with GBSA. Please set GBSA != 0;\n";
    die();
  }
  if(($r_gen->{"GB"} > 0) && ($r_gbo->{"GBSA"} == 6) && ($r_gbo->{"IGB"} != 66)){
    print "    GBSA = 6 can currently only be used with IGB = 66.\n";
    die();
  }

  if(($r_gen->{"PB"} > 0) && ($r_gen->{"DC"} > 0) && ($r_del->{"INDI"} != 1.0)){
     print "    PB decomposition can only be performed with INDI = 1.\n";
     die();
  }

  if($r_gen->{"GB"} > 0 && $r_gen->{"MS"} == 0){
    print "    Implicit SAS calc by sander\n";
  }
  #if($r_gen->{"PB"} > 0 && $r_gen->{"GB"} == 0 && $r_gen->{"MS"} == 0){
  #  print "    Without specifying GB, PB calc needs SAS calc by MS\n";
  #  die();
  #}

  if($r_gen->{"DC"} > 0 && ($r_gen->{"GC"} > 0 ||
                            $r_gen->{"AS"} > 0 || 
                            $r_gen->{"MS"} > 0)){
    print "    Energy decomp is not (yet) possible with GC, AS, MS\n";
    die("For details see: $HTMLPATH#check_sanity_no_enedec_gc_as_ms\n");
  }
  if($r_gen->{"DC"} > 0 && $r_gen->{"GB"} && $r_gbo->{"GBSA"} == 1){
    print "    Energy decomp only works with ICOSA algorithm (i.e. GBSA == 2)\n";
    die("For details see: $HTMLPATH#check_sanity_enedec_only_with_icosa\n");
  }
  if($r_gen->{"DC"} > 0 && $r_gen->{"NM"} && $r_nmo->{"PROC"} == 1){
    print "    Energy decomp only works with original nmode (i.e. PROC == 2 in NM section)\n";
    die("For details see: $HTMLPATH#check_sanity_enedec_only_with_ori_nmode\n");
  }
  if($r_gen->{"DC"} > 0 && $r_gen->{"NM"} && $r_dec->{"DCTYPE"} > 2){
    print "    Normal mode decomp only works for DCTYPE == 1 or 2\n";
    die("For details see: $HTMLPATH#check_sanity_enedec_nmode_only_perresidue\n");
  }

  if($r_gen->{"PATH"} !~ /\/$/){
    $r_gen->{"PATH"} .= "/";
  }
  if(! -e $r_gen->{"PATH"} &&
     ! -d $r_gen->{"PATH"}){
    print "    Cannot find directory ",$r_gen->{"PATH"},"\n";
    die("For details see: $HTMLPATH#check_sanity_cannot_find_dir\n");
  }

  if($r_gen->{"MM"} || $r_gen->{"GB"} || $r_gen->{"PB"} || $r_gen->{"MS"} || $r_gen->{"NM"}){
    # Test if approriate parmtop files exist
    if($r_gen->{"COMPLEX"} && 
      (! exists $r_gen->{"COMPT"} || ! defined $r_gen->{"COMPT"} || ! -e $r_gen->{"COMPT"})){
      print "    COMPT must be specified (correctly)\n";
      die("For details see: $HTMLPATH#check_sanity_compt_not_specified\n");
    }
    if($r_gen->{"RECEPTOR"} && 
      (! exists $r_gen->{"RECPT"} || ! defined $r_gen->{"RECPT"} || ! -e $r_gen->{"RECPT"})){
      print "    RECPT must be specified (correctly)\n";
      die("For details see: $HTMLPATH#check_sanity_recpt_not_specified\n");
    }
    if($r_gen->{"LIGAND"} && 
      (! exists $r_gen->{"LIGPT"} || ! defined $r_gen->{"LIGPT"} || ! -e $r_gen->{"LIGPT"})){
      print "    LIGPT must be specified (correctly)\n";
      die("For details see: $HTMLPATH#check_sanity_ligpt_not_specified\n");
    }
  }
 
  # Test DC parameters
  if($r_gen->{"DC"}){
    print "    Checking DC\n";
    &test_par(\@DC_PAR,[$r_dec]);

    if((&get_res_no($r_dec->{"COMREC"}) + &get_res_no($r_dec->{"COMLIG"})) != 
       &get_res_no($r_dec->{"COMPRI"})){
      die("    \#COMREC + \#COMLIG != \#COMPRI\n    For details see: $HTMLPATH#check_dc_comrec_n_comlig_lt_compri\n");
    }

    if(&get_res_no($r_dec->{"RECRES"}) != &get_res_no($r_dec->{"RECPRI"})){
      die("    \#RECRES != \#RECPRI\n    For details see: $HTMLPATH#check_dc_recres_lt_recpri\n");
    }

    if(&get_res_no($r_dec->{"RECMAP"}) > 0 && 
       &get_res_no($r_dec->{"RECMAP"}) != &get_res_no($r_dec->{"RECPRI"})){
      die("    \#RECMAP != \#RECPRI\n    For details see: $HTMLPATH#check_dc_recmap_unequ_recpri\n");
    }

    if(&get_res_no($r_dec->{"LIGRES"}) != &get_res_no($r_dec->{"LIGPRI"})){
      die("    \#LIGRES != \#LIGPRI\n    For details see: $HTMLPATH#check_dc_ligres_lt_ligpri\n");
    }

    if(&get_res_no($r_dec->{"LIGMAP"}) > 0 && 
       &get_res_no($r_dec->{"LIGMAP"}) != &get_res_no($r_dec->{"LIGPRI"})){
      die("    \#LIGMAP != \#LIGPRI\n    For details see: $HTMLPATH#check_dc_ligmap_unequ_ligpri\n");
    }

    if((&get_res_no($r_dec->{"RECMAP"}) + &get_res_no($r_dec->{"LIGMAP"})) !=
       &get_res_no($r_dec->{"COMPRI"})){
      die("    \#RECMAP + \#LIGMAP != \#COMPRI\n    For details see: $HTMLPATH#check_dc_recmap_n_ligmap_unequ_compri\n");
    }
  }

  # Test GC parameters and TRAJECTORY
  if($r_gen->{"GC"} || $r_gen->{"AS"}){
    print "    Checking GC\n";
    &test_par(\@GC_PAR,[$r_mak]);
    &test_fil(\@GC_FIL,[$r_pro]);

    if($r_gen->{"AS"}){
      print "    Checking AS\n";
      &test_par(\@AS_PAR,[$r_mak]);
      &test_fil(\@AS_FIL,[$r_pro]);
    }

    print "    Checking TRAJ\n";
    if(! defined $$r_tra || $$r_tra eq ""){
      die("    TRAJECTORY not defined\n    For details see: $HTMLPATH#check_traj_not_defined\n");
    }
    my $traj;
    foreach $traj (split/ +/,$$r_tra){
      if(! -e $traj){ 
        die("    $traj not found\n    For details see: $HTMLPATH#check_traj_not_found\n");
      }
    }
  }
  
  # Test MM parameters
  if($r_gen->{"MM"}){
    print "    Checking MM\n";
    &test_par(\@MM_PAR,[$r_san]);
    &test_fil(\@MM_FIL,[$r_pro]);
  }
  
  # Test NM parameters
  if($r_gen->{"NM"}){
    print "    Checking NM\n";
    &test_par(\@NM_PAR,[$r_nmo]);
    &test_fil(\@NM_FIL,[$r_pro]);

    if($r_nmo->{"PROC"} == 1){
      # Test NABnmode parameters
      if($r_nmo->{"IGB"} > 0){
	# Test GB parameters for NABnmode
	if($r_nmo->{"IGB"} != 1){
	  die("    If IGB>0, IGB must be 1 in NM section\n    For details see: $HTMLPATH#check_nm_igb_must_be_one\n");
	}
	if(! exists $r_nmo->{"SALTCON"} ||
	   ! exists $r_nmo->{"EXTDIEL"}||
	   ! exists $r_nmo->{"SURFTEN"}){
	  die("    SALTCON and EXTDIEL and SURFTEN must be set for IGB>0 in NM section\n    For details see: $HTMLPATH#check_nm_saltcon_extdiel_surften\n");
	}
      }
      else{
	# Test DIELC parameter for NABnmode
	if(! exists $r_nmo->{"DIELC"}){
	  die("    If IGB=0, DIELC must be set in NM section\n    For details see: $HTMLPATH#check_nm_dielc_not_set\n");
	}
      }
    }
    elsif($r_nmo->{"PROC"} == 2){
      # Test original nmode parameters
      if($r_nmo->{"IGB"} != 0){
	die("    If PROC= 2, IGB must be 0 in NM section\n    For details see: $HTMLPATH#check_nm_igb_must_be_zero\n");
      }
      else{
	if(! exists $r_nmo->{"DIELC"}){
	  die("    If IGB=0, DIELC must be set in NM section\n    For details see: $HTMLPATH#check_nm_dielc_must_be_set\n");
	}
      }
    }
    else{
      die("    PROC parameter in NM section must be 1 or 2\n    For details see: $HTMLPATH#check_nm_proc_must_ne_one_two\n");
    }
  }

  # Test PB parameters
  if($r_gen->{"PB"}){
    print "    Checking PB\n";
    &test_par(\@PB_PAR,[$r_del]);
    if($r_del->{"PROC"} == 1){
      # Test delphi parameters
      &test_par(\@DE_PAR,[$r_del]);
      &test_fil(\@DE_FIL,[$r_del,$r_pro]);

      # Test for enough SCALE and PERFIL parameters
      if($r_del->{"FOCUS"} > 0){
	if(scalar(@{$r_del->{"SCALE"}}) - 1 != $r_del->{"FOCUS"}){
	  die("    Nof SCALE parameters wrong\n    For details see: $HTMLPATH#check_pb_scale_wrong\n");
	}
	if(scalar(@{$r_del->{"PERFIL"}}) - 1 != $r_del->{"FOCUS"}){
	  die("    Nof PERFIL parameters wrong\n    For details see: $HTMLPATH#check_pb_perfil_wrong\n");
	}
      }

      # No focussing calc for REFE > 0
      if($r_del->{"REFE"} > 0 && $r_del->{"FOCUS"} > 0){
	die("    Focussing does not work with REFE > 0\n    For details see: $HTMLPATH#check_pb_focussing_bad_refe\n");
      }

      # Test for IVCAP == 1,5
      if(exists $r_del->{"IVCAP"} && defined $r_del->{"IVCAP"} && 
	 ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5)){
	die("    IVCAP == 1,5 only works with PBSA, not DELPHI.\n    For details see: $HTMLPATH#check_pb_ivcap_not_delphi\n");
      }

      # No energy decomp with delphi
      if($r_gen->{"DC"} > 0){
        die("    PB energy decomp does not work with delphi\n    For details see: $HTMLPATH#check_pb_no_dec_with_delphi\n");
      }

      # No parallel computation with delphi
      if(exists $r_gen->{"PARALLEL"} &&
	 defined $r_gen->{"PARALLEL"} &&
	 $r_gen->{"PARALLEL"} > 0){
        die("    Parallel computation does not work with delphi\n    For details see: $HTMLPATH#check_pb_no_parallel_with_delphi\n");
      }
	  
      if($r_gen->{"GB"} == 0 && $r_gen->{"MS"} == 0){
        print "    Without specifying GB, PB calc with Delphi needs SAS calc by MS\n";
        die("For details see: $HTMLPATH#check_pb_delphi_ms\n");
      } 
    }
    elsif($r_del->{"PROC"} == 2){
      # Backward compatibility, setup default value for newer options
      if(! exists $r_del->{"ARCRES"}){
	$r_del->{"ARCRES"}="0.25";
      }

      # Test pbsa parameters
      &test_par(\@RL_PAR,[$r_del]);
      &test_fil(\@RL_FIL,[$r_del,$r_pro]);

      # No focussing calc for REFE > 0
      if($r_del->{"REFE"} > 0){
	die("    PB calculation with REFE > 0 doesn't work with PROC equals 2\n    For details see: $HTMLPATH#check_pb_refe_proc\n");
      }

      # Tests for IVCAP == 1 or 5
      if(exists $r_del->{"IVCAP"} && defined $r_del->{"IVCAP"} && $r_del->{"IVCAP"} == 1 &&
	 ( ! exists $r_del->{"CUTCAP"} || ! defined $r_del->{"CUTCAP"} ||
	   ! exists $r_del->{"XCAP"}   || ! defined $r_del->{"XCAP"}   ||
	   ! exists $r_del->{"YCAP"}   || ! defined $r_del->{"YCAP"}   ||
	   ! exists $r_del->{"ZCAP"}   || ! defined $r_del->{"ZCAP"})){
	die("    For IVCAP == 1, CUTCAP, XCAP, YCAP, and ZCAP must also be given.\n    For details see: $HTMLPATH#check_pb_ivcap_one_params\n");
      }
      if(exists $r_del->{"IVCAP"} && defined $r_del->{"IVCAP"} && $r_del->{"IVCAP"} == 5 &&
	 ( ! exists $r_del->{"CUTCAP"} || ! defined $r_del->{"CUTCAP"})){
	die("    For IVCAP == 5, CUTCAP must also be given.\n    For details see: $HTMLPATH#check_pb_ivcap_five_cutcap\n");
      }
      if(exists $r_del->{"IVCAP"} && defined $r_del->{"IVCAP"} && 
	 ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5) &&
	 ( ! $r_gen->{"MS"})){
	die("    For IVCAP == 1,5, MS must be used to determine cavity energy.\n    For details see: $HTMLPATH#check_pb_ivcap_one_five_ms\n");
      }
      if(exists $r_del->{"IVCAP"} && defined $r_del->{"IVCAP"} && 
	 ($r_del->{"IVCAP"} == 1 || $r_del->{"IVCAP"} == 5) &&
	 $r_gen->{"MM"}){
	print("    For IVCAP == 1,5, MM cannot be used in parallel.\n");
	print("    Instead, do an IVCAP == 0 calculation with MM == 1 first, using snapshots without explicit waters.\n");
	print("    Then do an IVCAP == 5 calculation with MM == 0, using snapshots with explict waters.\n");
	die  ("    Results of both calculations currently need to be combined \"manually\".\n    For details see: $HTMLPATH#check_pb_ivcap_one_five_manually\n");
      }
      if((exists $r_del->{"RADIOPT"} && defined $r_del->{"RADIOPT"} && $r_del->{"RADIOPT"} == 0) &&
	 (exists $r_del->{"INP"} && defined $r_del->{"INP"} && $r_del->{"INP"} == 2)){
	die("    If INP == 2 then RADIOPT must be set to 1.\n    For details see: $HTMLPATH#check_pb_inp_radiopt\n");
      }
      if($r_gen->{"DC"} && (exists $r_del->{"INP"} && defined $r_del->{"INP"} && $r_del->{"INP"} != 0)){
	die("    For decomposition with pbsa INP must be set to 0.\n    For details see: $HTMLPATH#check_pb_dec_inp\n");
      }
    }
    elsif($r_del->{"PROC"} == 3){
      # Test iapbs parameters
      &test_par(\@RL_PAR,[$r_del]);
      &test_fil(\@RL_FIL,[$r_del,$r_pro]);

      # No focussing calc for REFE > 0
      #if($r_del->{"REFE"} > 0){
	#die("    PB calculation with REFE > 0 doesn't work with PROC equals 2\n");
      #}

      # No energy decomp with iapbs
      if($r_gen->{"DC"} > 0){
        die("    PB energy decomp does not work with iapbs\n    For details see: $HTMLPATH#check_pb_no_pbdec_iapbs\n");
      }

    }
    else{
      die("    PROC parameter in PB section must be 1, 2 or 3\n    For details see: $HTMLPATH#check_pb_proc_one_two_three\n");
    }

    # Test if DIELC and INDI agree (if REFE > 0, no test is performed)
    if($r_del->{"REFE"} == 0 && $r_del->{"INDI"} != $r_san->{"DIELC"}){
      die("    INDI != DIELC\n    For details see: $HTMLPATH#check_pb_indi_unequ_dielc\n");
    }

    # Test if INDI > 1.0 if REFE > 0
    if($r_del->{"REFE"} > 0 && $r_del->{"INDI"} <= 1.0){
      die("    INDI <= 1.0 does not make sense for REFE > 0\n    For details see: $HTMLPATH#check_pb_indi_lt_one\n");
    }

    # No ISTRNG > 0.0 for REFE > 0
    if($r_del->{"REFE"} > 0 && exists $r_del->{"ISTRNG"} && $r_del->{"ISTRNG"} > 0.0){
      die("    ISTRNG > 0.0 does not (yet) work with REFE > 0\n    For details see: $HTMLPATH#check_pb_instrng_gt_zero\n");
    }
  }

  # Test GB parameters
  if($r_gen->{"GB"}){
    print "    Checking GB\n";
    &test_par(\@GB_PAR,[$r_gbo]);
    &test_fil(\@GB_FIL,[$r_gbo,$r_pro]);

    if($r_gbo->{"INTDIEL"} > 1.0 && $r_san->{"DIELC"} > 1.0){
      die("    DIELC must be 1.0 if INTDIEL > 1.0\n    For details see: $HTMLPATH#check_gb_dielc_one_intdiel_one\n");
    }
  }

  # Test MS parameters
  if($r_gen->{"MS"}){
    print "    Checking MS\n";
    &test_par(\@MS_PAR,[$r_mol]);
    &test_fil(\@MS_FIL,[$r_mol,$r_pro]);
  }
  return;
}

sub get_res_no(){
#################

  # Parameters: $str
  my $str = shift;

  my @tmp = split/ +/,$str;
  my $res_no = 0;
  my $tmp;
  foreach $tmp (@tmp){
    if(($tmp =~ /^\d+$/ && $tmp == 0) ||
       ($tmp =~ /0-0/)){
      last;
    }
    else{
      if($tmp =~ /(\d+)-(\d+)/){
        $res_no += $2 - $1 + 1;
      }
      else{
        die("    Wrong format for decomp: $str\n    For details see: $HTMLPATH#get_resno_wrong_format\n");
      }
    }
  }
  return $res_no;
}

sub test_par(){
###############

  # Parameters: \@,\@
  my $r_must = shift;
  my $r_is   = shift;

  my $must;
  foreach $must (@$r_must){
    my $is;

    my $notexists = 1;
    foreach $is (@$r_is){ 
      if(exists $is->{$must}){
        $notexists = 0;
        last;
      }
    }
    if($notexists){
      die("Param $must does not exist\nFor details see: $HTMLPATH#testpar_param_not_existent\n");
    }
    
    my $notdefined = 1;
    foreach $is (@$r_is){
      if(defined $is->{$must}){
        $notdefined = 0;
        last;
      }
    }
    if($notdefined){
      die("Param $must is not defined\nFor details see: $HTMLPATH#testpar_param_not_defined\n");
    }
  }
  return;
} 

sub test_fil(){
###############

  # Parameters: \@,\@
  my $r_must = shift;
  my $r_is   = shift;

  my $must;
  foreach $must (@$r_must){
    my $is;
  
    my $notfound = 1;
    foreach $is (@$r_is){
      if(exists $is->{$must} && -e $is->{$must}){
        $notfound = 0;
        last;
      }
    }
    if($notfound){
      die("File $must is not found\nFor details see: $HTMLPATH#test_file_not_found\n");
    } 
  }
  return;
}

1; # Necessary for package function
