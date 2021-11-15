use strict;

# Genration of copies of pdb-files used for system setup in topo-folder and
# deletion of crystal water from these files -> Required for later setup of
# topologies for calculations with mm_pbsa.pl
sub delete_cryst_wat{
	my $topo_dir = shift;
	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $pdb_file;
	my $wat_found = 0;

	for my $tag ("rec", "com", "lig"){
	
		$pdb_file = $c_in{'root_path'}."/leap/".$struct_b."/".$struct_b."_".$tag.".pdb";

		my $new_pdb_file = $topo_dir."/".$struct_b."_".$tag.".pdb";
		my $pdb_no_wat = $topo_dir."/".$struct_b."_".$tag."_noWAT.pdb";

		if($tag eq "lig"){
			copy($pdb_file, $pdb_no_wat) || die "Cannot copy $pdb_file to $pdb_no_wat.\n";
		}

		if($tag ne "lig"){
			copy($pdb_file, $new_pdb_file) || die "Cannot copy $pdb_file to $new_pdb_file.\n";

			open(PDB, $new_pdb_file) || die "Cannot open pdb-file $new_pdb_file.\n";
			open(PDB_NO_WAT, ">$pdb_no_wat") || die "Cannot open pdb-file $pdb_no_wat for writing.\n";

			$wat_found = 0;
			
			while(my $pdb_line = <PDB>){
				chomp($pdb_line);
				my @pdb_line = split(/\s+/, $pdb_line);
				
				if(($pdb_line[3] eq "WAT")||($pdb_line[3] eq "HOH")){
					$wat_found = 1;
					next;
				}
				elsif(($pdb_line[0] eq "TER")&&($wat_found == 1)){
					$wat_found = 0;
					next;
				}
				else{
					print PDB_NO_WAT $pdb_line."\n";
				}
			}
			close PDB_NO_WAT;
			close PDB;
			unlink($new_pdb_file);
		}
	}
}

# For MM-PBSA calculations with implicit membrane:
# Generation of copies of pdb-files used for system setup in topo-folder and
# deletion of all residues except the receptor residues (and the ligand in case
# of the complex) from these files -> Required for later setup of topologies
# for calculations with mm_pbsa.pl
sub delete_lipids_ions_wat{
	my $topo_dir = shift;
	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $pdb_file;

	for my $tag ("rec", "com", "lig"){
		$pdb_file = $c_in{'root_path'}."/leap/".$struct_b."/".$struct_b."_".$tag.".pdb";

		my $new_pdb_file = $topo_dir."/".$struct_b."_".$tag.".pdb";
		my $pdb_no_mem = $topo_dir."/".$struct_b."_".$tag."_noMem.pdb";

		if($tag eq "lig"){
			copy($pdb_file, $pdb_no_mem) || die "Cannot copy $pdb_file to $pdb_no_mem.\n";
		}

		if($tag ne "lig"){
			copy($pdb_file, $new_pdb_file) || die "Cannot copy $pdb_file to $new_pdb_file.\n";

			open(PDB, $new_pdb_file) || die "Cannot open pdb-file $new_pdb_file.\n";
			open(PDB_NO_MEM, ">$pdb_no_mem") || die "Cannot open pdb-file $pdb_no_mem for writing.\n";
			my $res_count = 0;
			my $prev_resno = "";
			my $prev_resname = "";

			while(my $pdb_line = <PDB>){
				chomp($pdb_line);
				if(($pdb_line =~ m/ATOM/)||($pdb_line =~ m/HETATM/)){
						my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z) =
					unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8", $pdb_line);
					$resno =~ s/^\s+//g;
					$resno =~ s/\s+$//g;
					$resname =~ s/^\s+//g;
					$resname =~ s/\s+$//g;
		
					if(($resname ne $prev_resname)||($resno ne $prev_resno)){
						$res_count++;
					}
		
					if(($resno <= $c_in{'rec_res'})&&($res_count <= $c_in{'rec_res'})){
						print PDB_NO_MEM $pdb_line."\n";
					}
					if((($resno == ($c_in{'rec_res'}+1))&&($tag eq "com"))||(($res_count == ($c_in{'rec_res'}+1))&&($tag eq "com"))){
						print PDB_NO_MEM $pdb_line."\n";
					}
					if((($resno > ($c_in{'rec_res'}+1))&&($tag eq "com"))||(($resno > $c_in{'rec_res'})&&($tag eq "rec"))
					||(($res_count > ($c_in{'rec_res'}+1))&&($tag eq "com"))||(($res_count > $c_in{'rec_res'})&&($tag eq "rec"))){
						last;
					}
					$prev_resname = $resname;
					$prev_resno = $resno;
				}
				if($pdb_line =~ m/^TER/){
					print PDB_NO_MEM $pdb_line."\n"; 
				}
			}
			close PDB;
			close PDB_NO_MEM;
			unlink($new_pdb_file);
		}
	}
}
#START LW_03_2018
# For protein-protein MM-PBSA calculations with implicit membrane:
# Generation of copies of pdb-files used for system setup in topo-folder and
# deletion of all except for receptor/ligand/complex from these files ->
# Required for later setup of topologies for calculations with mm_pbsa.pl
# Rename caps if necessary for correct leap input
sub prot_delete_lipids_ions_wat{
    my $struct_path = shift;
	my $topo_dir = shift;
	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $pdb_file;
	
	my ($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
	my ($lig_start, $lig_end) = split(/-/, $c_in{'lig_range'});

	for my $tag ("rec", "com", "lig"){
        $pdb_file = $struct_path."/".$struct_b.".pdb";

        my $new_pdb_file = $topo_dir."/".$struct_b."_".$tag.".pdb";
        my $pdb_no_mem = $topo_dir."/".$struct_b."_".$tag."_noMem.pdb";

        copy($pdb_file, $new_pdb_file) || die "Cannot copy $pdb_file to $new_pdb_file.\n";

        open(PDB, $new_pdb_file) || die "Cannot open pdb-file $new_pdb_file.\n";
        open(PDB_NO_MEM, ">$pdb_no_mem") || die "Cannot open pdb-file $pdb_no_mem for writing.\n";
        my $res_count = 0;
        my $prev_resno = "";
        my $prev_resname = "";

        while(my $pdb_line = <PDB>){
            chomp($pdb_line);
            if(($pdb_line =~ m/ATOM/)||($pdb_line =~ m/HETATM/)){
                    my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z) =
                unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8", $pdb_line);
                $resno =~ s/^\s+//g;
                $resno =~ s/\s+$//g;
                $resname =~ s/^\s+//g;
                $resname =~ s/\s+$//g;

                if(($resname ne $prev_resname)||($resno ne $prev_resno)){
                    $res_count++;
                }

                if(($resno <= $rec_end)&&($res_count <= $rec_end)&&($tag ne "lig")){
                    # Don't use ACE/NME atoms not recognized by leap
                    if((($resname eq "ACE")&&($atom =~ m/HH3/))||(($resname eq "NME")&&($atom =~ m/^H/))){
                    }
                    else{
                        print PDB_NO_MEM $pdb_line."\n";
                    }
                }
                if((($resno > ($lig_end))&&($tag eq "com"))||(($resno > $rec_end)&&($tag eq "rec")||(($resno > $lig_end)&&($tag eq "lig")))
                ||(($res_count > ($lig_end))&&($tag eq "com"))||(($res_count > $rec_end+1)&&($tag eq "rec"))||(($res_count > $lig_end+1)&&($tag eq "lig"))){
                    last;
                }
                if((($resno >= ($rec_end+1))&&($tag eq "com"))||(($res_count >= ($rec_end+1))&&($tag eq "com"))
                ||(($resno >= ($rec_end+1))&&($tag eq "lig"))||(($res_count >= ($rec_end+1))&&($tag eq "lig"))){
                    # Don't use ACE/NME atoms not recognized by leap
                    if((($resname eq "ACE")&&($atom =~ m/HH3/))||(($resname eq "NME")&&($atom =~ m/^H/))){
                    }
                    else{
                        print PDB_NO_MEM $pdb_line."\n";
                    }
                }
                $prev_resname = $resname;
                $prev_resno = $resno;
            }
            if($pdb_line =~ m/^TER/){
                print PDB_NO_MEM $pdb_line."\n"; 
            }
        }
        close PDB;
        close PDB_NO_MEM;
        unlink($new_pdb_file);
    }
}
#END LW_03_2018

# Setup leap script for topology and coordinate file generation for MM-PB(GB)SA calculations
sub gen_leap_in_top_setup{
	my $topo_dir = shift;
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	# Copy input script used for MD topology setup to current directory
	my $leap_vac_MD = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct."/cryst/leap_script_1.in";
	my $new_leap_in = $topo_dir."/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".in";
	my $mod_library = $c_in{'root_path'}."/leap/".$struct."/".$struct."_".$c_in{'chrg_meth'}."_mod.lib";
	my $frcmod_file = $c_in{'root_path'}."/leap/".$struct."/".$struct.".frcmod";
	my $first_lib = 0;
	my $first_frcmod = 0;

	open(LEAP_MD, $leap_vac_MD) || die "Cannot open leap-script $leap_vac_MD.\n";
	open(NEW_LEAP, ">$new_leap_in") || die "Cannot open leap-script $new_leap_in for writing.\n";

	while(my $l = <LEAP_MD>){
		chomp($l);
		# Ensure backward compartibility due to
		# force field support in AMBER14
		if(($l =~ m/source leaprc.ff99SB/)||($l =~ m/source leaprc.ff12SB/)){
			if($c_in{'backwards'} == 0){
				$l =~ s/leaprc.ff99SB/leaprc.protein.ff14SB/;
				$l =~ s/leaprc.ff12SB/leaprc.protein.ff14SB/;
				print NEW_LEAP $l."\n";
		                if($c_in{'water_model'} eq "TIP3P"){
					print NEW_LEAP "source leaprc.water.tip3p\n";
				}
				if($c_in{'water_model'} eq "OPC"){
					print NEW_LEAP "source leaprc.water.opc\n";
				}
			}
			else{
				print NEW_LEAP "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff99SB\n";
                                # need to source tip3p because we are removing tip3p parameters from the protein ff
                                # It already makes the assumption in this program that tip3p is default so keeping it same
	                        print NEW_LEAP "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";

			}
		}
		if(($c_in{'backwards'} == 0)&&($l =~ m/loadAmberParams frcmod.ionsjc_tip3p/)){
			next;
		}	
		# Set radii
		elsif($l =~ m/source leaprc.gaff/){
			print NEW_LEAP $l."\n";
			if($c_in{'GB'} == 1){
				print NEW_LEAP "set default PBradii mbondi\n";
			}
			if(($c_in{'GB'} == 2)||($c_in{'GB'} == 5)){
				print NEW_LEAP "set default PBradii mbondi2\n";
			}
			if($c_in{'GB'} == 6){
				print NEW_LEAP "set default PBradii bondi\n";
			}
			if(($c_in{'PB'} == 1)&&($c_in{'GB'} == 0)){
				print NEW_LEAP "set default PBradii mbondi\n";
			}
			if($c_in{'PB'} == 2){
				print NEW_LEAP "set default PBradii mbondi\n";
			}
			if($c_in{'PB'} == 3){
				print NEW_LEAP "set default PBradii parse\n";
			}
			if(($c_in{'PB'} == 4)&&($c_in{'GB'} == 0)){
				print NEW_LEAP "set default PBradii mbondi\n";
			}
		}
		# Reading of library files
		elsif(($first_lib == 0)&&($l =~ m/loadoff/)){
			if(-e $mod_library){
				print NEW_LEAP "loadoff $mod_library\n";
				$first_lib = 1;
			}
			else{
				print NEW_LEAP "loadoff ".$c_in{'root_path'}."/leap/".$struct."/".$struct."_".$c_in{'chrg_meth'}.".lib\n";
				$first_lib = 1;
			}
		}
		elsif(($first_lib == 1)&&($l =~ m/loadoff/)&&($c_in{'add_lib'})&&($c_in{'add_lib'} ne "")){
			print NEW_LEAP "loadoff ".$c_in{'add_lib'}."\n";
		}
		elsif(($first_lib == 1)&&($l =~ m/loadoff/)&&((!$c_in{'add_lib'})||($c_in{'add_lib'} eq ""))){
			print "\nWARNING: As you provided no information about the additional library file\n";
			print "in the command file, it will be assumed that this file can be found in the\n";
			print "same location as in the setup of the MD simulations.\n";
			print NEW_LEAP $l."\n";
		}
		# Reading of frcmod files
		elsif(($first_frcmod == 0)&&($l =~ m/loadAmberParams/)&&(!($l =~ m/frcmod.ionsjc_tip3p/))){
			if(-e $frcmod_file){
				print NEW_LEAP "loadAmberParams $frcmod_file\n";
				$first_frcmod = 1;
			}
		}
		elsif(($first_frcmod == 1)&&($l =~ m/loadAmberParams/)&&($c_in{'add_frcmod'})&&($c_in{'add_frcmod'} ne "")){
			print NEW_LEAP "loadAmberParams ".$c_in{'add_frcmod'}."\n";
		}
		elsif(($first_frcmod == 1)&&($l =~ m/loadAmberParams/)&&((!$c_in{'add_frcmod'})||($c_in{'add_frcmod'} eq ""))){
			print "\nWARNING: As you provided no information about the additional parameter (frcmod)\n";
			print "file in the command file, it will be assumed that this file can be found in the\n";
			print "same location as in the setup of the MD simulations.\n";
			print NEW_LEAP $l."\n";
		}
		
		# Case 1: No implicit membrane and no crystal water molecules
		elsif(($c_in{'ImplMem'} == 0)&&((! exists $c_in{'wat'})||($c_in{'wat'} == 0))&&($l =~ m/^COM/)){
			print NEW_LEAP "COM = loadpdb ".$c_in{'root_path'}."/leap/".$struct."/".$struct."_com.pdb\n";
		}
		elsif(($c_in{'ImplMem'} == 0)&&((! exists $c_in{'wat'})||($c_in{'wat'} == 0))&&($l =~ m/^REC/)){
			print NEW_LEAP "REC = loadpdb ".$c_in{'root_path'}."/leap/".$struct."/".$struct."_rec.pdb\n";
		}
		
		# Case 2: No implicit membrane and crystal water molecules present
		elsif(($c_in{'ImplMem'} == 0)&&((exists $c_in{'wat'})&&($c_in{'wat'}==1))&&($l =~ m/^COM/)){
			print NEW_LEAP "COM = loadpdb ".$topo_dir."/".$struct."_com_noWAT.pdb\n";
		}
		elsif(($c_in{'ImplMem'} == 0)&&((exists $c_in{'wat'})&&($c_in{'wat'}==1))&&($l =~ m/^REC/)){
			print NEW_LEAP "REC = loadpdb ".$topo_dir."/".$struct."_rec_noWAT.pdb\n";
		}
		
		# Case 3: Implicit membrane
		elsif(($c_in{'ImplMem'} == 1)&&($l =~ m/^COM/)){
			print NEW_LEAP "COM = loadpdb ".$topo_dir."/".$struct."_com_noMem.pdb\n";
		}
		elsif(($c_in{'ImplMem'} == 1)&&($l =~ m/^REC/)){
			print NEW_LEAP "REC = loadpdb ".$topo_dir."/".$struct."_rec_noMem.pdb\n";
		}
		
		# Case 4: No implicit membrane calculation requested, but membrane was present
		#         during the simulation and had to be removed for the MM-PB(GB)SA analysis
		elsif(($c_in{'ImplMem'} == 0)&&((exists $c_in{'memResNo'})&&($c_in{'memResNo'} != 0))&&($l =~ m/^COM/)){
			print NEW_LEAP "COM = loadpdb ".$topo_dir."/".$struct."_com_noMem.pdb\n";
		}
		elsif(($c_in{'ImplMem'} == 0)&&((exists $c_in{'memResNo'})&&($c_in{'memResNo'} != 0))&&($l =~ m/^REC/)){
			print NEW_LEAP "REC = loadpdb ".$topo_dir."/".$struct."_rec_noMem.pdb\n";
		}
	
		elsif($l =~ m/^LIG/){
			print NEW_LEAP "LIG = loadpdb ".$c_in{'root_path'}."/leap/".$struct."/".$struct."_lig.pdb\n";
		}
		elsif($l =~ m/^charge/){
			next;
		}
		elsif($l =~ m/^saveAmberParm LIG/){
			print NEW_LEAP "saveAmberParm LIG ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_lig.top ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_lig.crd\n";
		}
		elsif($l =~ m/^saveAmberParm REC/){
			print NEW_LEAP "saveAmberParm REC ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_rec.top ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_rec.crd\n";
		}
		elsif($l =~ m/^saveAmberParm COM/){
			print NEW_LEAP "saveAmberParm COM ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_com.top ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_com.crd\n";
		}
		elsif($l =~ m/^savepdb LIG/){
			print NEW_LEAP "savepdb LIG ".$topo_dir."/".$struct."_vac_lig.pdb\n";
		}
		elsif($l =~ m/^savepdb COM/){
			print NEW_LEAP "savepdb COM ".$topo_dir."/".$struct."_vac_com.pdb\n";
		}
		elsif($l =~ m/^savepdb REC/){
			print NEW_LEAP "savepdb REC ".$topo_dir."/".$struct."_vac_rec.pdb\n";
		}
		else{
			print NEW_LEAP $l."\n";
		}
	}
	close NEW_LEAP;
	close LEAP_MD;
}

#START LW_03_2018
# Create new leap input script for protein_ligand=1
# since it is post-processing
sub prot_gen_leap_in_top_setup{
	my $topo_dir = shift;
	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my %disulf;
	
	my $new_leap_in = $topo_dir."/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".in";
	my @unit_names = ("COM", "REC", "LIG");

	open(LEAP_IN, ">$new_leap_in") || die "Cannot open $new_leap_in";

	if($c_in{'backwards'} == 0){
		print LEAP_IN "source leaprc.protein.ff14SB\n";
		if($c_in{'water_model'} eq "TIP3P"){
			print LEAP_IN "source leaprc.water.tip3p\n";
		}
		if($c_in{'water_model'} eq "OPC"){
			print LEAP_IN "source leaprc.water.opc\n";
		}
	}
	else{
		print LEAP_IN "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff99SB\n";
                # need to source tip3p because we are removing tip3p parameters from the protein ff
	        print LEAP_IN "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";
	}
	print LEAP_IN "source leaprc.gaff\n";
			
	if(($c_in{'add_lib'})&&(-e $c_in{'add_lib'})){
		print LEAP_IN "loadoff $c_in{'add_lib'}\n";
	}
	if(($c_in{'add_frcmod'})&&(-e $c_in{'add_frcmod'})){
		print LEAP_IN "loadAmberParams $c_in{'add_frcmod'}\n";
	}
	if((exists $c_in{'cys_f'})&&(-e $c_in{'cys_f'})){
		my $disulf = read_disulf(\%c_in);
		%disulf = %{$disulf};
	}
	print LEAP_IN "LIG = loadpdb ".$topo_dir."/".$struct_b."_lig_noMem.pdb\n";
	print LEAP_IN "REC = loadpdb ".$topo_dir."/".$struct_b."_rec_noMem.pdb\n";
	print LEAP_IN "COM = loadpdb ".$topo_dir."/".$struct_b."_com_noMem.pdb\n";
	
	if(keys %disulf){
		foreach my $unit_name (@unit_names){
			foreach my $k (keys %disulf){
				print LEAP_IN "bond ".$unit_name.".".$k.".SG ".$unit_name.".".$disulf{$k}.".SG\n";
			}
		}
	}
	print LEAP_IN "saveAmberParm LIG ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_lig.top ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_lig.crd\n";
	print LEAP_IN "savepdb LIG ".$topo_dir."/".$struct_b."_vac_lig.pdb\n";
	print LEAP_IN "saveAmberParm REC ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_rec.top ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_rec.crd\n";
    print LEAP_IN "savepdb REC ".$topo_dir."/".$struct_b."_vac_rec.pdb\n";
	print LEAP_IN "saveAmberParm COM ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_com.top ".$topo_dir."/".$struct_b."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_com.crd\n";
    print LEAP_IN "savepdb COM ".$topo_dir."/".$struct_b."_vac_com.pdb\n";
	print LEAP_IN "charge LIG\n";
	print LEAP_IN "charge REC\n";
	print LEAP_IN "charge COM\n";
	print LEAP_IN "quit\n";
	close LEAP_IN;
}

#END LW_03_2018


# Copying topology files with water for PB=2
sub copy_topWAT{
	my $topo_dir = shift;
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	for my $tag ("com", "rec", "lig"){
		# Copy topology file used for MD to current directory
		my $q_solv_top = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct."/cryst/".$struct."_solv_".$tag.".top";
		my $t_solv_top = $topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_solv_".$tag.".top";

		copy($q_solv_top, $t_solv_top) || die "Cannot copy $q_solv_top to $t_solv_top.\n";
	}
}


# Identify number of atoms of receptor and ligand from pdb file of complex generated during setup of MD-input
# Required for extraction of coordinates from trajectories with mm_pbsa.pl
sub identify_atom_no{
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	my $lig_atom_start;
	my $lig_atom_end;
	my $rec_atom_start;
	my $rec_atom_end;
	my $total_atoms;

	my $atom_count = 0;
	my $enum = 0;
	my $rec_res_ref = $c_in{'rec_res'} + 1;

	my $pdb_compl = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct."/cryst/".$struct."_solv_com.pdb";

	open(PDB_COM, $pdb_compl) || die "Cannot open file $pdb_compl.\n";

	while(my $compl_line = <PDB_COM>){
		chomp($compl_line);
		my @compl_line = split(/\s+/, $compl_line);

		if($compl_line[0] eq "ATOM"){

			# Assign first atom of receptor
			if($enum == 0){
				$rec_atom_start = $compl_line[1];
				$enum = 1;
			}
			
			# Assign last atom of receptor
			if(($enum == 1)&&($compl_line[4] == $rec_res_ref)){
				$rec_atom_end = $atom_count;
			}

			# Assign first atom of ligand
			if(($enum == 1)&&(($compl_line[3] eq "<1>")||($compl_line[3] eq "MOL")||($compl_line[3] eq "UNN"))){
				$lig_atom_start = $compl_line[1];
				$enum = 2;
			}

			# Assign last atom of ligand
			if(($enum == 2)&&($compl_line[4] > $rec_res_ref)){
				$lig_atom_end = $atom_count;
				$enum = 3;
			}
			
			$atom_count++;
		}

		if($compl_line[0] eq "END"){
			$total_atoms = $atom_count;
		}
	}

	return ($lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, $total_atoms);
}


#START LW_03_2018
# Identify number of atoms of receptor and ligand from pdb file of complex given by the user
# Required for extraction of coordinates from trajectories with mm_pbsa.pl
# 
sub prot_identify_atom_no{
    my $struct_path = shift;
	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	my $lig_atom_start;
	my $lig_atom_end;
	my $rec_atom_start;
	my $rec_atom_end;
	my $total_atoms;

	my $atom_count = 0;
	my $enum = 0;
	
    my ($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
	my ($lig_start, $lig_end) = split(/-/, $c_in{'lig_range'});
	
	$rec_end +=1;

    my $pdb_compl = $struct_path."/".$struct_b.".pdb";
    
	open(PDB_COM, $pdb_compl) || die "Cannot open file $pdb_compl.\n";

	while(my $compl_line = <PDB_COM>){
		chomp($compl_line);
		my @compl_line = split(/\s+/, $compl_line);

		if($compl_line[0] eq "ATOM"){

			# Assign first atom of receptor
			if($enum == 0){
				$rec_atom_start = $compl_line[1];
				$enum = 1;
			}
			
			# Assign last atom of receptor
			if(($enum == 1)&&($compl_line[4] == $rec_end)){
				$rec_atom_end = $atom_count;
				$enum = 2;
			}

			# Assign first atom of ligand
			if($enum == 2){
				$lig_atom_start = $compl_line[1];
				$enum = 3;
			}
			
			# Assign last atom of ligand
			if(($enum == 3)&&($compl_line[4] == $lig_end)){
				$lig_atom_end =  $atom_count+1;
			}
			$atom_count++;
		}
		if($compl_line[0] eq "END"){
			$total_atoms = $atom_count;
        }
	}

	return ($lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, $total_atoms);
}
#END LW_03_2018

# Determine total atom number based on pdb-files of solvated systems
sub determine_total_atom_no{
	my $tag = shift;
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $atom_count = 0;
	my $pdb_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct."/cryst/".$struct."_solv_".$tag.".pdb";

	open(PDB, $pdb_file) || die "Cannot open $pdb_file for reading atom number.\n";

	while(my $pdb_l = <PDB>){
		chomp($pdb_l);

		if($pdb_l =~ m/^ATOM/){
			$atom_count++;
		}
	}
	close PDB;
	return $atom_count;
}


# Image receptor and ligand molecules to origin
sub image_traj{
	my $traj_file = shift;
	my $file_count = shift;
	my $res_no = shift;
	my $traj_path = shift;
	my $top = shift;
	my $backwards_compartibility = shift;

	# Setup cpptraj script for imaging
	my @traj_file_s = split(/\./, $traj_file);
	open(CPPTRAJ, ">".$traj_path."/cpptraj_img_".$file_count) || die "Cannot open file ".$traj_path."/cpptraj_img_".$file_count." for writing.\n";

	print CPPTRAJ "trajin ".$traj_path."/".$traj_file."\n";
	if($backwards_compartibility == 0){
		print CPPTRAJ "center :1-".$res_no." origin\n";
	}
	else{
		print CPPTRAJ "center :1-".$res_no." origin mass\n";
	}
	print CPPTRAJ "image origin center\n";
	print CPPTRAJ "trajout ".$traj_path."/".$traj_file_s[0]."_img.mdcrd\n";
	close CPPTRAJ;
	
	my $cpptraj_log = $traj_path."/cpptraj_img_".$file_count.".log";
	my $cpptraj_command = "cpptraj -p ".$top." -i ".$traj_path."/cpptraj_img_".$file_count;
	$cpptraj_command .= " > ".$cpptraj_log;
	
	my $cpptraj_file = "$traj_path/cpptraj_img_$file_count";
	my $cpptraj_out = $traj_path."/".$traj_file_s[0]."_img.mdcrd";
	return ($cpptraj_command, $cpptraj_log, $cpptraj_file, $cpptraj_out);
}


# Changing command-file for extraction of coordinates with mm_pbsa.pl
sub modify_extract_snaps{
	my $snaps_dir = shift;
	my $lig_atom_start = shift;
	my $lig_atom_end = shift;
	my $rec_atom_start = shift;
	my $rec_atom_end = shift;
	my $total_atoms_ref = shift;
	my %total_atoms = %{$total_atoms_ref};
	my $struct = shift;
	my $traj_files_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $traj = 0;
	my %traj_path;
	my @tags;
	my %top;

	if($c_in{'traj'} == 1){
		push(@tags, "com");
	}
	if($c_in{'traj'} == 3){
		@tags = ("com", "lig", "rec");
	}

	# Set path to trajectory files
	foreach my $tag (@tags){
		if($tag eq "rec"){
			$traj_path{$tag} = "../../../MD_".$c_in{'chrg_meth'}."/".$tag."/prod";
		}
		else{
			$traj_path{$tag} = "../../../MD_".$c_in{'chrg_meth'}."/".$struct."/".$tag."/prod";
		}
	}

	# Determine topology files - Same location for 1- and 3-trajectory approach
	for my $tag ("lig", "com", "rec"){
		$top{$tag} = "../".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_".$tag.".top";
	}
	

	# Change command file
	foreach my $tag (@tags){
		open(COORD_TEMPL, $c_in{'coord_templ'}) || die "Cannot open template file ".$c_in{'coord_templ'}.".\n";
		open(COORD_IN, ">$snaps_dir/extract_coordinates_$tag.in") || die "Cannot open file $snaps_dir/extract_coordinates_$tag.in!\n";

		while(my $coord = <COORD_TEMPL>){
			chomp($coord);

			if($coord =~ m/^PREFIX/){
				$coord =~ s/Structure/$struct/;
			}
			if($coord =~ m/^PATH/){
				$coord =~ s/\/Path\/to\/snapshot\/dir/\.\//;
			}
			# Set to 1 only in case of 1-trajectory approach, in 3-trajectory
			# apprach complex is treated as receptor
			if($coord =~ m/^COMPLEX/){
				if($c_in{'traj'} == 1){
					$coord =~ s/0/1/;
				}
			}
		 	# Set to 1 if 1-trajectory approach and for complex, receptor, and
			# ligand in case of 3-trajecotry approach
			if($coord =~ m/^RECEPTOR/){
				$coord =~ s/0/1/;
			}
			# Set to 1 in case of 1-trajectory approch. In 3-trajectory approach
			# ligand is treated as receptor
			if($coord =~ m/^LIGAND/){
				if($c_in{'traj'} == 1){
					$coord =~ s/0/1/;
				}
			}
			if($coord =~ m/^COMPT/){
				$coord =~ s/Complex\.top/$top{'com'}/;
			}
			if($coord =~ m/^RECPT/){
				$coord =~ s/Receptor\.top/$top{'rec'}/;
			}
			if($coord =~ m/^LIGPT/){
				$coord =~ s/Ligand\.top/$top{'lig'}/;
			}
			if($coord =~ m/^NTOTAL/){
				$coord =~ s/0/$total_atoms{$tag}/;
			}
			if($coord =~ m/^NSTART/){
				$coord =~ s/0/$c_in{'Nstart'}/;
			}
			if($coord =~ m/^NSTOP/){
				$coord =~ s/0/$c_in{'Nstop'}/;
			}
			if($coord =~ m/^NFREQ/){
				$coord =~ s/1/$c_in{'Nfreq'}/;
			}
			# Set to zero in all cases for traj=3 since
			# Ligand treated as receptor
			if($coord =~ m/^NUMBER_LIG_GROUPS/){
				if($c_in{'traj'} == 3){
					$coord =~ s/1/0/;
				}
			}
			if($coord =~ m/^LSTART/){
				# Not required for 3-trj, but incorporated for consistency
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					$coord =~ s/0/1/;
				}
				else{
					$coord =~ s/0/$lig_atom_start/;
				}
			}
			if($coord =~ m/^LSTOP/){
				# Not required for 3-trj, but incorporated for consistency
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					my $end = ($lig_atom_end - $lig_atom_start) + 1;
					$coord =~ s/0/$end/;
				}
				else{
					$coord =~ s/0/$lig_atom_end/;
				}
			}
			if($coord =~ m/^RSTART/){
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					$coord =~ s/0/1/;
				}
				else{
					$coord =~ s/0/$rec_atom_start/;
				}
			}
			if($coord =~ m/^RSTOP/){
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					my $end = ($lig_atom_end - $lig_atom_start) + 1;
					$coord =~ s/0/$end/;
				}
				elsif(($c_in{'traj'} == 3)&&($tag eq "com")){
					$coord =~ s/0/$lig_atom_end/;
				}
				else{
					$coord =~ s/0/$rec_atom_end/;
				}
			}
			if($traj == 1){
				my $check_comment = substr($coord, 0, 1);

				if($check_comment ne "#"){
					foreach my $f (@traj_files){
						my $traj_f;
						if($c_in{'image'} == 1){
							my @file_base = split(/\./, $f);
							$traj_f = $traj_path{$tag}."/".$file_base[0]."_img.mdcrd.gz";
						}
						else{
							$traj_f = $traj_path{$tag}."/".$f;
						}

						print COORD_IN "TRAJECTORY          $traj_f\n";
					}
					$traj = 0;
					next;
				}
			}
			if($coord eq "\@TRAJECTORY"){
				$traj = 1;
			}
			print COORD_IN $coord."\n";

		}	# For loop
	}	# End reading template command file
	close COORD_IN;
	close COORD_TEMPL;
}


# Changing command-file for extraction of snapshots with water for PB=2
sub modify_extract_snaps_solv{
	my $solv_snaps_dir = shift;
	my $total_atoms_ref = shift;
	my %total_atoms = %{$total_atoms_ref};
	my $struct = shift;
	my $tag_ref = shift;
	my @tags = @{$tag_ref};
	my $traj_path_ref = shift;
	my %traj_path = %{$traj_path_ref};
	my $traj_files_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $traj = 0;
	my %top;

	# Determine topology files
	foreach my $tag (@tags){
		$top{$tag} = "../../".$struct."_gb".$c_in{'GB'}."_solv_".$tag.".top";
	}

	# Change command file
	foreach my $tag (@tags){
		open(COORD_TEMPL, $c_in{'coord_templ'}) || die "Cannot open template file ".$c_in{'coord_templ'}.".\n";
		open(COORD_IN, ">$solv_snaps_dir/$tag/extract_coordinates_$tag.in") || die "Cannot open file $solv_snaps_dir/$tag/extract_coordinates_$tag.in!\n";

		while(my $coord = <COORD_TEMPL>){
			chomp($coord);

			if($coord =~ m/^PREFIX/){
				$coord =~ s/Structure/$struct/;
			}
			if($coord =~ m/^PATH/){
				$coord =~ s/\/Path\/to\/snapshot\/dir/\.\//;
			}
		 	# Set to 1 for complex, receptor, and ligand
			if($coord =~ m/^RECEPTOR/){
				$coord =~ s/0/1/;
			}
			if($coord =~ m/^COMPT/){
				$coord =~ s/Complex\.top/$top{'com'}/;
			}
			if($coord =~ m/^RECPT/){
				$coord =~ s/Receptor\.top/$top{'rec'}/;
			}
			if($coord =~ m/^LIGPT/){
				$coord =~ s/Ligand\.top/$top{'lig'}/;
			}
			if($coord =~ m/^NTOTAL/){
				$coord =~ s/0/$total_atoms{$tag}/;
			}
			if($coord =~ m/^NSTART/){
				$coord =~ s/0/$c_in{'Nstart'}/;
			}
			if($coord =~ m/^NSTOP/){
				$coord =~ s/0/$c_in{'Nstop'}/;
			}
			if($coord =~ m/^NFREQ/){
				$coord =~ s/1/$c_in{'Nfreq'}/;
			}
			# Set to zero in all cases since ligand is treated as receptor
			if($coord =~ m/^NUMBER_LIG_GROUPS/){
				if($c_in{'traj'} == 3){
					$coord =~ s/1/0/;
				}
			}
			# NUMBER_REC_GROUPS remains 1 in all cases
			if($coord =~ m/^RSTART/){
				$coord =~ s/0/1/;
			}
			if($coord =~ m/^RSTOP/){
				$coord =~ s/0/$total_atoms{$tag}/;
			}
			if($traj == 1){
				my $check_comment = substr($coord, 0, 1);

				if($check_comment ne "#"){
					foreach my $f (@traj_files){
						my $traj_f;
						if($c_in{'image'} == 1){
							my @file_base = split(/\./, $f);
							$traj_f = $traj_path{$tag}."/".$file_base[0]."_img.mdcrd.gz";
						}
						else{
							$traj_f = $traj_path{$tag}."/".$f;
						}

						print COORD_IN "TRAJECTORY          $traj_f\n";
					}
					$traj = 0;
					next;
				}
			}
			if($coord eq "\@TRAJECTORY"){
				$traj = 1;
			}
			print COORD_IN $coord."\n";

		}	# For loop
	}	# End reading template command file
	close COORD_IN;
	close COORD_TEMPL;
}

#START LW_03_2018
sub prot_modify_extract_snaps{
	my $snaps_dir = shift;
	my $lig_atom_start = shift;
	my $lig_atom_end = shift;
	my $rec_atom_start = shift;
	my $rec_atom_end = shift;
	my $total_atoms_ref = shift;
	my %total_atoms = %{$total_atoms_ref};
	my $struct = shift;
	my $traj_files_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $traj = 0;
	my %traj_path;
	my @tags;
	my %top;

	if($c_in{'traj'} == 1){
		push(@tags, "com");
	}

	# Set path to trajectory files
	foreach my $tag (@tags){
		$traj_path{$tag} = "";
	}

	# Determine topology files - Same location for 1- and 3-trajectory approach
	for my $tag ("lig", "com", "rec"){
		$top{$tag} = "../".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_".$tag.".top";
	}
	

	# Change command file
	foreach my $tag (@tags){
		open(COORD_TEMPL, $c_in{'coord_templ'}) || die "Cannot open template file ".$c_in{'coord_templ'}.".\n";
		open(COORD_IN, ">$snaps_dir/extract_coordinates_$tag.in") || die "Cannot open file $snaps_dir/extract_coordinates_$tag.in!\n";

		while(my $coord = <COORD_TEMPL>){
			chomp($coord);

			if($coord =~ m/^PREFIX/){
				$coord =~ s/Structure/$struct/;
			}
			if($coord =~ m/^PATH/){
				$coord =~ s/\/Path\/to\/snapshot\/dir/\.\//;
			}
			# Set to 1 only in case of 1-trajectory approach, in 3-trajectory
			# apprach complex is treated as receptor
			if($coord =~ m/^COMPLEX/){
				if($c_in{'traj'} == 1){
					$coord =~ s/0/1/;
				}
			}
		 	# Set to 1 if 1-trajectory approach and for complex, receptor, and
			# ligand in case of 3-trajecotry approach
			if($coord =~ m/^RECEPTOR/){
				$coord =~ s/0/1/;
			}
			# Set to 1 in case of 1-trajectory approch. In 3-trajectory approach
			# ligand is treated as receptor
			if($coord =~ m/^LIGAND/){
				if($c_in{'traj'} == 1){
					$coord =~ s/0/1/;
				}
			}	
            if($coord =~ m/^COMPT/){
				$coord =~ s/Complex\.top//;
			}
			if($coord =~ m/^RECPT/){
				$coord =~ s/Receptor\.top//;
			}
			if($coord =~ m/^LIGPT/){
				$coord =~ s/Ligand\.top//;
			}
			if($coord =~ m/^NTOTAL/){
				$coord =~ s/0/$total_atoms{$tag}/;
			}
			if($coord =~ m/^NSTART/){
				$coord =~ s/0/$c_in{'Nstart'}/;
			}
			if($coord =~ m/^NSTOP/){
				$coord =~ s/0/$c_in{'Nstop'}/;
			}
			if($coord =~ m/^NFREQ/){
				$coord =~ s/1/$c_in{'Nfreq'}/;
			}
			# Set to zero in all cases for traj=3 since
			# Ligand treated as receptor
			if($coord =~ m/^NUMBER_LIG_GROUPS/){
				if($c_in{'traj'} == 3){
					$coord =~ s/1/0/;
				}
			}
			if($coord =~ m/^LSTART/){
				# Not required for 3-trj, but incorporated for consistency
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					$coord =~ s/0/1/;
				}
				else{
					$coord =~ s/0/$lig_atom_start/;
				}
			}
			if($coord =~ m/^LSTOP/){
				# Not required for 3-trj, but incorporated for consistency
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					my $end = ($lig_atom_end - $lig_atom_start) + 1;
					$coord =~ s/0/$end/;
				}
				else{
					$coord =~ s/0/$lig_atom_end/;
				}
			}
			if($coord =~ m/^RSTART/){
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					$coord =~ s/0/1/;
				}
				else{
					$coord =~ s/0/$rec_atom_start/;
				}
			}
			if($coord =~ m/^RSTOP/){
				if(($c_in{'traj'} == 3)&&($tag eq "lig")){
					my $end = ($lig_atom_end - $lig_atom_start) + 1;
					$coord =~ s/0/$end/;
				}
				elsif(($c_in{'traj'} == 3)&&($tag eq "com")){
					$coord =~ s/0/$lig_atom_end/;
				}
				else{
					$coord =~ s/0/$rec_atom_end/;
				}
			}
			if($traj == 1){
				my $check_comment = substr($coord, 0, 1);

				if($check_comment ne "#"){
					foreach my $f (@traj_files){
						my $traj_f;
						
						#TODO:
						#change img-function and see how path has to be changed here
# 						if($c_in{'image'} == 1){
# 							my @file_base = split(/\./, $f);
# 							$traj_f = $traj_path{$tag}."/".$file_base[0]."_img.mdcrd.gz";
# 						}
# 						else{
							$traj_f = $f;
						

						print COORD_IN "TRAJECTORY          $traj_f\n";
					}
					$traj = 0;
					next;
				}
			}
			if($coord eq "\@TRAJECTORY"){
				$traj = 1;
			}
			print COORD_IN $coord."\n";

		}	# For loop
	}	# End reading template command file
	close COORD_IN;
	close COORD_TEMPL;
}
#END LW_03_2018

# Generate cpptraj-script for extraction of snapshots from trajectories
# imaged to receptor
sub generate_cpptraj_rec_img_pqr{
	my $pqr_dir = shift;
	my $traj_file = shift;
	my $traj_no = shift;
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $first_mem_residue = $c_in{'rec_res'} + 2;
	my $last_mem_residue = $c_in{'rec_res'} + 2 + $c_in{'memResNo'};
	
	my $pqr_cpptraj = $pqr_dir."/pqr_".$traj_no.".ptraj";
	open(PQR_CPPTRAJ, ">".$pqr_cpptraj) || die "Cannot open file $pqr_cpptraj for writing.\n";
	print PQR_CPPTRAJ "trajin ".$c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct."/com/prod/".$traj_file."\n";
	print PQR_CPPTRAJ "center :1-".$c_in{'rec_res'}." origin\n";
	print PQR_CPPTRAJ "image origin center\n\n";
	print PQR_CPPTRAJ "strip :Na+\n";
	print PQR_CPPTRAJ "strip :Cl-\n";
	print PQR_CPPTRAJ "strip :WAT\n";
	print PQR_CPPTRAJ "strip :".$first_mem_residue."-".$last_mem_residue."\n\n";
	print PQR_CPPTRAJ "trajout tmp/snap.pqr pdb multi parse\n";
	close PQR_CPPTRAJ;
	
	return $pqr_cpptraj;
}

#START LW_03_2018
sub prot_generate_cpptraj_rec_img_pqr{
    my $snaps_dir =shift;
    my $topo_dir =shift;
	my $pqr_dir = shift;
	my $struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	my ($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
	my ($lig_start, $lig_end) = split(/-/, $c_in{'lig_range'});
	my $snaps_count = 0;
	
	# set traj_no to 1
	my $pqr_cpptraj = $pqr_dir."/pqr_1.ptraj";
	open(PQR_CPPTRAJ, ">".$pqr_cpptraj) || die "Cannot open file $pqr_cpptraj for writing.\n";
	print PQR_CPPTRAJ "parm ".$topo_dir."/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_com.top\n";
	
	for(my $n=1; $n<=$c_in{'Nstop'}; $n++){
		if($n <= $c_in{'Nstop'}){
			if(($n >= $c_in{'Nstart'})&&(($n % $c_in{'Nfreq'}) == 0)){
				$snaps_count = $snaps_count + 1;
				print PQR_CPPTRAJ "trajin ".$snaps_dir."/".$struct."_com.crd.".$snaps_count."\n";
			}
		}
		else{
			last;
		}
	}
	print PQR_CPPTRAJ "center :".$rec_start."-".$lig_end." origin\n";
	print PQR_CPPTRAJ "image origin center\n\n";
	print PQR_CPPTRAJ "trajout tmp/snap.pqr pdb multi parse\n";
	close PQR_CPPTRAJ;
	
	return $pqr_cpptraj;
}
#END LW_03_2018

# Extraction of receptor and ligand from complex PQR
sub gen_rec_lig_pqr{
	my $pqr_dir = shift;
	my $struct = shift;
	my $com_pqr = shift;
	my $rec_res = shift;
	my $struct_count = shift;
	
	my $rec_str = "";
	my $lig_str = "";
	my $read_rec = 1;
	
	# Read complex structure information
	open(COM_PQR, $com_pqr) || die "Cannot open file $com_pqr for reading.\n";
	while(my $com_pqr_l = <COM_PQR>){
		chomp($com_pqr_l);
	
		if(($com_pqr_l =~ m/ATOM/)||($com_pqr_l =~ m/HETATM/)){
	    		my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z) =
			unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8", $com_pqr_l);
	
			$resno =~ s/^\s+//g;
			$resno =~ s/\s+$//g;
		
			if($resno <= $rec_res){
				$rec_str .= $com_pqr_l."\n";
			}
			else{
				$lig_str .= $com_pqr_l."\n";
				$read_rec = 0;
			}
		}
	
		if(($com_pqr_l =~ m/TER/)&&($read_rec == 1)){
			$rec_str .= "TER\n";
		}
		if(($com_pqr_l =~ m/TER/)&&($read_rec == 0)){
			$lig_str .= "TER\n";
		}
	}
	close COM_PQR;
	
	# Write receptor structure
	my $rec_pqr = $pqr_dir."/".$struct."_rec.pqr.".$struct_count;
	open(REC_PQR, ">$rec_pqr") || die "Cannot open file $rec_pqr for writing.\n";
	print REC_PQR $rec_str;
	close REC_PQR;	
	
	# Write ligand structure
	my $lig_pqr = $pqr_dir."/".$struct."_lig.pqr.".$struct_count;
	open(LIG_PQR, ">$lig_pqr") || die "Cannot open file $lig_pqr for writing.\n";
	print LIG_PQR $lig_str;
	close LIG_PQR;	
}

#START LW_03_2018
# For prot_prot = 1, rename lig pqr atom/res numbers
sub prot_gen_rec_lig_pqr{
	my $pqr_dir = shift;
	my $struct = shift;
	my $com_pqr = shift;
	my $rec_res = shift;
	my $struct_count = shift;
	my $rec_str = "";
	my $lig_str = "";
	my $read_rec = 1;
	
	my $lig_first_atom;
	my $lig_first_res;
	my $read_first = 1;
	
	# Read complex structure information
	open(COM_PQR, $com_pqr) || die "Cannot open file $com_pqr for reading.\n";
	while(my $com_pqr_l = <COM_PQR>){
		chomp($com_pqr_l);
	
		if(($com_pqr_l =~ m/ATOM/)||($com_pqr_l =~ m/HETATM/)){
	    		my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z,   $charge, $radius,$ele) =
			unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8  x  a8        a8     a8", $com_pqr_l);
			
			$resno =~ s/^\s+//g;
			$resno =~ s/\s+$//g;
 
			if($resno <= $rec_res){
				$rec_str .= $com_pqr_l."\n";
			}
			else{
				$read_rec = 0;

				# subtract first atom and res number for numbering to start at 1
				if ($read_first == 1){
					$lig_first_atom = $atnum-1;
					$lig_first_res = $resno-1;
					$read_first = 0;
				}
				$atnum = $atnum - $lig_first_atom;
				$resno = $resno - $lig_first_res;

				$lig_str .= sprintf("%6s%5s%5s%4s%6s%12s%8s%8s%9s%8s%8s\n",$card,$atnum,$atom,$resname,$resno,$x,$y,$z,$charge,$radius,$ele);
			}
		}
	
		if(($com_pqr_l =~ m/TER/)&&($read_rec == 1)){
			$rec_str .= "TER\n";
		}
		if(($com_pqr_l =~ m/TER/)&&($read_rec == 0)){
			$lig_str .= "TER\n";
		}
	}
	close COM_PQR;
	
	# Write receptor structure
	my $rec_pqr = $pqr_dir."/".$struct."_rec.pqr.".$struct_count;
	open(REC_PQR, ">$rec_pqr") || die "Cannot open file $rec_pqr for writing.\n";
	print REC_PQR $rec_str;
	close REC_PQR;	
	
	# Write ligand structure
	my $lig_pqr = $pqr_dir."/".$struct."_lig.pqr.".$struct_count;
	open(LIG_PQR, ">$lig_pqr") || die "Cannot open file $lig_pqr for writing.\n";
	print LIG_PQR $lig_str;
	close LIG_PQR;	
}
#END LW_03_2018


# Changing command-file for MM-PBSA/MM-GBSA calculation
sub modify_mmpbsa_in{

	my $calc_dir = shift;
	my $mmpbsa_in = shift;
	my $struct = shift;
	my $inter = shift;
	my $traj_inter_ref = shift;
	my %traj_inter = %{$traj_inter_ref};
	my $traj_files_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	# Set parameters for Decomposition
	my $rec_res_tot;
	my $lig_res_in_compl;

	#START LW_03_2018	
	# Set parameters for protein-ligand
	my $rec_start;
	my $rec_end;
	my $lig_start;
	my $lig_end;
	#END LW_03_2018

	if(($c_in{'Decomp'} > 0)&&(exists $c_in{'rec_res'})){
		$rec_res_tot = $c_in{'rec_res'};
		$lig_res_in_compl = $rec_res_tot+1;
	}
	#START LW_03_2018
	elsif(($c_in{'Decomp'} > 0)&&($c_in{'prot_prot'} == 1)){
		($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
		($lig_start, $lig_end) = split(/-/, $c_in{'lig_range'});
	}

	# Determine topology files
	my %top;

	for my $tag ("lig", "com", "rec"){
		$top{$tag} = "topo/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_".$tag.".top";
	}

	# Change command file
	my $pb_surften = 0;
	my $pb_surfoff = 0;

	open(MMPBSA_TEMPL, $c_in{'mmpbsa_templ'}) || die "Cannot open template file ".$c_in{'mmpbsa_templ'}.".\n";
	open(MMPBSA_IN, ">$calc_dir/mmpbsa.in") || die "Cannot open file $calc_dir/mmpbsa.in!\n";

	while(my $mmpbsa_c = <MMPBSA_TEMPL>){
		chomp($mmpbsa_c);

		if($mmpbsa_c =~ m/^PREFIX/){
			$mmpbsa_c =~ s/Structure/$struct/;
		}
		if($mmpbsa_c =~ m/^PARALLEL/){
			if($c_in{'parallel'} ne ""){
				$mmpbsa_c =~ s/0/$c_in{'parallel'}/;
			}
		}
		if($mmpbsa_c =~ m/^PATH/){
			$mmpbsa_c =~ s/\/path\/to\/snapshot\/dir/.\/snapshots/;
		}
		if($mmpbsa_c =~ m/^COMPLEX/){
			$mmpbsa_c =~ s/0/1/;
		}
		if($mmpbsa_c =~ m/^RECEPTOR/){
			$mmpbsa_c =~ s/0/1/;
		}
		if($mmpbsa_c =~ m/^LIGAND/){
			$mmpbsa_c =~ s/0/1/;
		}
		if($mmpbsa_c =~ m/^COMPT/){
			$mmpbsa_c =~ s/Complex\.top/$top{'com'}/;
		}
		if($mmpbsa_c =~ m/^RECPT/){
			$mmpbsa_c =~ s/Receptor\.top/$top{'rec'}/;
		}
		if($mmpbsa_c =~ m/^LIGPT/){
			$mmpbsa_c =~ s/Ligand\.top/$top{'lig'}/;
		}
		if($mmpbsa_c =~ m/^START/){
			my $start = $traj_inter{$inter}->{'Start'};
			$mmpbsa_c =~ s/1/$start/;
		}
		if($mmpbsa_c =~ m/^STOP/){
			my $stop = $traj_inter{$inter}->{'Stop'};
			$mmpbsa_c =~ s/1/$stop/;
		}
		if($mmpbsa_c =~ m/^OFFSET/){
			my $offset = $traj_inter{$inter}->{'Offset'};
			$mmpbsa_c =~ s/1/$offset/;
		}
		
		# Set PB to 0 if no Poisson-Boltzmann calculation is requested
		if($c_in{'PB'} == 0){
			if($mmpbsa_c =~ m/^PB /){
				$mmpbsa_c =~ s/1/$c_in{'PB'}/;
			}
		}
		
		# Set GB to 0, if no Generalized Born calculation is requested, or
		# to ensure method consistency
		if($c_in{'GB'} == 0){
			if($mmpbsa_c =~ m/^GB /){
				$mmpbsa_c =~ s/1/$c_in{'GB'}/;
			}
		}
		
		# Implicit solvent model : PB
		
		# PB 1 : Solvation free energy calculated according to the method 
		#        of Tan et al., J. Phys. Chem. B, 2007, 111, 12263-12274.
		#        Radii used for calculation of polar contribution: Tan&Luo radii for
		#        standard residues and mbondi for non-standard residues, i.e. ligands.
		if($c_in{'PB'} == 1){
			# PB remains 1
			# Set MS to zero
			if($mmpbsa_c =~ m/^MS/){
				$mmpbsa_c =~ s/1/0/;
			}
			# Set INP to 2
			if($mmpbsa_c =~ m/^INP/){
				$mmpbsa_c =~ s/0/2/;
			}
			# Set RADIOPT to 1
			if($mmpbsa_c =~ m/^RADIOPT/){
				$mmpbsa_c =~ s/0/1/;
			}
			# set SURFTEN to 0.03780
			if(($mmpbsa_c =~ m/^SURFTEN/)&&($pb_surften == 0)){
				$mmpbsa_c =~ s/1\.0/0\.03780/;
			}
			# set SURFOFF to -0.5692
			if(($mmpbsa_c =~ m/^SURFOFF/)&&($pb_surfoff == 0)){
				$mmpbsa_c =~ s/0\.0/-0\.5692/;
			}
		}
		# PB 3 : Non-polar part of solvation free energy calculated with molsurf
		#        Polar part of solvation free energy calculated with PBSA
		#        Radii used for calculation of polar contribution: Parse radii
		if($c_in{'PB'} == 3){
			# MS remains 1
			# RADIOPT remains 0
			# IVCAP remains 0
			# INP not regarded
			# set SURFTEN to 0.00542
			if(($mmpbsa_c =~ m/^SURFTEN/)&&($pb_surften == 0)){
				$mmpbsa_c =~ s/1\.0/0\.00542/;
			}
			# set SURFOFF to 0.92
			if(($mmpbsa_c =~ m/^SURFOFF/)&&($pb_surfoff == 0)){
				$mmpbsa_c =~ s/0\.0/0\.92/;
			}
		}	
		# PB 4 : Non-polar part of solvation free energy calculated with molsurf
		#        Polar part of solvation free energy calculated with PBSA
		#        Radii used for calculation of polar contribution: mbondi
		if($c_in{'PB'} == 4){
			# PB remains 1
			# MS remains 1
			# RADIOPT remains 0
			# IVCAP remains 0
			# INP not regarded
			# set SURFTEN to 0.0072
			if(($mmpbsa_c =~ m/^SURFTEN/)&&($pb_surften == 0)){
				$mmpbsa_c =~ s/1\.0/0\.0072/;
			}
			# SURFOFF remains 0.0
		}

		# Implicit solvent model : GB
		if(($c_in{'GB'} == 1)||($c_in{'GB'} == 2)||($c_in{'GB'} == 5)){
			# GB remains 1
			# Set IGB to corresponding method
			if($mmpbsa_c =~ m/^IGB/){
				$mmpbsa_c =~ s/1/$c_in{'GB'}/;
			}
			# Set MS to 0 to use LCPO method
			if($mmpbsa_c =~ /^MS/){
				$mmpbsa_c =~ s/1/0/;
			}
		}
		if($c_in{'GB'} == 6){
			if($mmpbsa_c =~ /^GBSA/){
				$mmpbsa_c =~ s/1/6/;
			}
			if($mmpbsa_c =~ /^IGB/){
				$mmpbsa_c =~ s/1/66/;
			}
			if($mmpbsa_c =~ /^EXTDIEL/){
				$mmpbsa_c =~ s/80.0/78.5/;
			}
		}
		if(($c_in{'GB'} == 2)||($c_in{'GB'} == 5)||($c_in{'GB'} == 6)){
			if(($mmpbsa_c =~ m/^SURFTEN/)&&($pb_surften == 1)){
				$mmpbsa_c =~ s/0\.0072/0\.005/;
			}
			# SURFOFF remains 0.0
		}
			
		if(($c_in{'Decomp'} > 0)&&(($c_in{'GB'} != 0)||($c_in{'ImplMem'} == 1))&&($c_in{'prot_prot'} != 1)){
			if($mmpbsa_c =~ m/^DC /){
				$mmpbsa_c =~ s/0/1/;
			}
			if($mmpbsa_c =~ m/^MS/){
				$mmpbsa_c =~ s/1/0/;
			}
			if($mmpbsa_c =~ m/^DCTYPE/){
				$mmpbsa_c =~ s/0/$c_in{'Decomp'}/;
			}
			if($mmpbsa_c =~ m/^COMREC/){
				$mmpbsa_c =~ s/0/1\-$rec_res_tot/;
			}
			if($mmpbsa_c =~ m/^COMLIG/){
				$mmpbsa_c =~ s/0/$lig_res_in_compl\-$lig_res_in_compl/;
			}
			if($mmpbsa_c =~ m/^COMPRI/){
				my $tmp_str = "1-".$rec_res_tot." ".$lig_res_in_compl."-".$lig_res_in_compl."\n";
				$mmpbsa_c =~ s/0/$tmp_str/;
			}
			if($mmpbsa_c =~ m/^RECRES/){
				$mmpbsa_c =~ s/0/1\-$rec_res_tot/;
			}
			if($mmpbsa_c =~ m/^RECPRI/){
				$mmpbsa_c =~ s/0/1\-$rec_res_tot/;
			}
			if($mmpbsa_c =~ m/^RECMAP/){
				$mmpbsa_c =~ s/0/1\-$rec_res_tot/;
			}
			if($mmpbsa_c =~ m/^LIGMAP/){
				$mmpbsa_c =~ s/0/$lig_res_in_compl\-$lig_res_in_compl/;
			}
			if($mmpbsa_c =~ m/^GBSA/){
				$mmpbsa_c =~ s/1/2/;
			}
		}
	
		#START LW_03_2018	
		# for prot_prot = 1
		if(($c_in{'Decomp'} > 0)&&(($c_in{'GB'} != 0)||($c_in{'ImplMem'} == 1))&&($c_in{'prot_prot'} == 1)){
			if($mmpbsa_c =~ m/^DC /){
				$mmpbsa_c =~ s/0/1/;
			}
			if($mmpbsa_c =~ m/^MS/){
				$mmpbsa_c =~ s/1/0/;
			}
			if($mmpbsa_c =~ m/^DCTYPE/){
				$mmpbsa_c =~ s/0/$c_in{'Decomp'}/;
			}
			if($mmpbsa_c =~ m/^COMREC/){
				$mmpbsa_c =~ s/0/1\-$rec_end/;
			}
			if($mmpbsa_c =~ m/^COMLIG/){
				$mmpbsa_c =~ s/0/$lig_start\-$lig_end/;
			}
			if($mmpbsa_c =~ m/^COMPRI/){
				my $tmp_str = "1-".$rec_end." ".$lig_start."-".$lig_end."\n";
				$mmpbsa_c =~ s/0/$tmp_str/;
			}
			if($mmpbsa_c =~ m/^RECRES/){
				$mmpbsa_c =~ s/0/1\-$rec_end/;
			}
			if($mmpbsa_c =~ m/^RECPRI/){
				$mmpbsa_c =~ s/0/1\-$rec_end/;
			}
			if($mmpbsa_c =~ m/^LIGRES/){
				my $lig_sol_start = $lig_start - $rec_end;
				my $lig_sol_end = $lig_end - $rec_end;
				$mmpbsa_c =~ s/1-1/$lig_sol_start-$lig_sol_end/;
			}
			if($mmpbsa_c =~ m/^LIGPRI/){
				my $lig_sol_start = $lig_start - $rec_end;
				my $lig_sol_end = $lig_end - $rec_end;
				$mmpbsa_c =~ s/1-1/$lig_sol_start-$lig_sol_end/;
			}
			if($mmpbsa_c =~ m/^RECMAP/){
				$mmpbsa_c =~ s/0/1\-$rec_end/;
			}
			if($mmpbsa_c =~ m/^LIGMAP/){
				$mmpbsa_c =~ s/0/$lig_start\-$lig_end/;
			}
			if($mmpbsa_c =~ m/^GBSA/){
				$mmpbsa_c =~ s/1/2/;
			}
		}
		#END LW_03_2018
		
		# Set program executable
		if(($c_in{'sander_exe'})&&(-e $c_in{'sander_exe'})){
			if($mmpbsa_c =~ m/Additional program executable/){
				$mmpbsa_c = $mmpbsa_c."\n".$c_in{'sander_exe'};
			}
		}
		
		# Check pass of PB surften and surfoff
		if($mmpbsa_c =~ m/^SURFTEN/){
			$pb_surften = 1;
		}			
		if($mmpbsa_c =~ m/^SURFOFF/){
			$pb_surfoff = 1;
		}	
	
		if($c_in{'ImplMem'} == 1){
			if($mmpbsa_c =~ m/^PB/){
				$mmpbsa_c =~ s/1/0/;
			}
			if($mmpbsa_c =~ m/^GB/){
				if($c_in{'Decomp'} == 0){
					# In case of normal analysis set GB to 0.
				}
				if($c_in{'Decomp'} > 0){
					# In case of implicit membrane decomposition analysis use GB>0,
					# otherwise electrostatic interactions, van der Waals energies, and SASA
					# are not calculated on a per-residue basis.
					$mmpbsa_c =~ s/0/1/;
				}
			}
		}		

		print MMPBSA_IN $mmpbsa_c . "\n";

	} # End read mmpbsa-template file
	close MMPBSA_IN;
	close MMPBSA_TEMPL;
}


sub modify_mmpbsa_in_pb2{

	my $calc_dir = shift;
	my $pb2_tag = shift;
	my $struct = shift;
	my $inter = shift;
	my $traj_inter_ref = shift;
	my %traj_inter = %{$traj_inter_ref};
	my $traj_files_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	# Determine topology files
	my %top;

	for my $tag ("lig", "com", "rec"){
		if($pb2_tag eq "gas"){
			$top{$tag} = "topo/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_".$tag.".top";
		}
		else{
			$top{$tag} = "topo/".$struct."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_solv_".$tag.".top";
		}
	}

	for my $tag ("com", "rec", "lig"){

		# Change command file
		my $surften_last = 0;
		my $surfoff_last = 0;

		open(MMPBSA_TEMPL, $c_in{'mmpbsa_templ'}) || die "Cannot open template file ".$c_in{'mmpbsa_templ'}.".\n";
		open(MMPBSA_IN, ">$calc_dir/$tag/mmpbsa.in") || die "Cannot open file $calc_dir/$tag/mmpbsa_$tag.in!\n";

		while(my $mmpbsa_c = <MMPBSA_TEMPL>){
			chomp($mmpbsa_c);

			if($mmpbsa_c =~ m/^PREFIX/){
				$mmpbsa_c =~ s/Structure/$struct/;
			}
			if(($mmpbsa_c =~ m/^PATH/)&&($pb2_tag eq "hyd")){
				$mmpbsa_c =~ s/\/path\/to\/snapshot\/dir/.\/snapshots\/$tag/;
			}
			if(($mmpbsa_c =~ m/^PATH/)&&($pb2_tag eq "gas")){
				$mmpbsa_c =~ s/\/path\/to\/snapshot\/dir/.\/snapshots/;
			}
			if($pb2_tag eq "gas"){
				if(($mmpbsa_c =~ m/^COMPLEX/)&&($tag eq "com")){
					$mmpbsa_c =~ s/0/1/;
				}
				if(($mmpbsa_c =~ m/^RECEPTOR/)&&($tag eq "rec")){
					$mmpbsa_c =~ s/0/1/;
				}
				if(($mmpbsa_c =~ m/^LIGAND/)&&($tag eq "lig")){
					$mmpbsa_c =~ s/0/1/;
				}
				if($mmpbsa_c =~ m/^COMPT/){
					$mmpbsa_c =~ s/Complex\.top/$top{'com'}/;
				}
				if($mmpbsa_c =~ m/^RECPT/){
					$mmpbsa_c =~ s/Receptor\.top/$top{'rec'}/;
				}
				if($mmpbsa_c =~ m/^LIGPT/){
					$mmpbsa_c =~ s/Ligand\.top/$top{'lig'}/;
				}
			}
			if($pb2_tag eq "hyd"){
				if($mmpbsa_c =~ m/^RECEPTOR/){
					$mmpbsa_c =~ s/0/1/;
				}
				if($mmpbsa_c =~ m/^RECPT/){
					$mmpbsa_c =~ s/Receptor\.top/$top{$tag}/;
				}
			}
			if($mmpbsa_c =~ m/^START/){
				my $start = $traj_inter{$inter}->{'Start'};
				$mmpbsa_c =~ s/1/$start/;
			}
			if($mmpbsa_c =~ m/^STOP/){
				my $stop = $traj_inter{$inter}->{'Stop'};
				$mmpbsa_c =~ s/1/$stop/;
			}
			if($mmpbsa_c =~ m/^OFFSET/){
				my $offset = $traj_inter{$inter}->{'Offset'};
				$mmpbsa_c =~ s/1/$offset/;
			}
			if($pb2_tag eq "gas"){
				# set PB to 0
				if($mmpbsa_c =~ m/^PB/){
					$mmpbsa_c =~ s/1/0/;
				}
				# set MS to 0
				if($mmpbsa_c =~ m/^MS/){
					$mmpbsa_c =~ s/1/0/;
				}
				# MM remains 1
			}
			if($pb2_tag eq "hyd"){
				# set MM to 0
				if($mmpbsa_c =~ m/^MM/){
					$mmpbsa_c =~ s/1/0/;
				}
				# Set IVCAP to 5
				if($mmpbsa_c =~ m/^IVCAP/){
					$mmpbsa_c =~ s/0/5/;
				}
				# Set CUTCAP to 50
				if($mmpbsa_c =~ m/^CUTCAP/){
					$mmpbsa_c =~ s/\-1\.0/50\.0/;
				}
				# set PROBE to 1.4
				if($mmpbsa_c =~ m/^PROBE/){
					$mmpbsa_c =~ s/0\.0/1\.4/;
				}
			}

			# set GB to 0
			if($mmpbsa_c =~ m/^GB/){
				$mmpbsa_c =~ s/1/0/;
			}
			# set SURFTEN to 0.069
			if(($mmpbsa_c =~ m/^SURFTEN/)&&($surften_last == 0)){
				$mmpbsa_c =~ s/1\.0/0\.069/;
				$surften_last = 1;
			}
			# SURFOFF remains 0.0
			if(($mmpbsa_c =~ m/^SURFOFF/)&&($surfoff_last == 0)){
				$surfoff_last = 1;
			}
			# Set INP to 1
			if($mmpbsa_c =~ m/^INP/){
				$mmpbsa_c =~ s/0/1/;
			}
			# set RADIOPT to 1
			if($mmpbsa_c =~ m/^RADIOPT/){
				$mmpbsa_c =~ s/0/1/;
			}
			
			# Set program executable
			if(($c_in{'sander_exe'})&&(-e $c_in{'sander_exe'})){
				if($mmpbsa_c =~ m/Additional program executable/){
					$mmpbsa_c = $mmpbsa_c."\n".$c_in{'sander_exe'};
				}
			}
			
			print MMPBSA_IN $mmpbsa_c . "\n";

		} # End read mmpbsa-template file
	} # End tag
	close MMPBSA_IN;
	close MMPBSA_TEMPL;
}

# Generate input file for mmpbsa_FEWmem
sub gen_mmpbsa_few_mem_input{
	my $calc_dir = shift;
	my $currentFEW_dir = shift;
	my $struct = shift;
	my $start = shift;
	my $stop = shift;
	my $offset = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	open(APBS_FEW_IN, ">".$calc_dir."/mmpbsa_FEWmem.in") || die "Cannot open file $calc_dir/mmpbsa_FEWmem.in for writing.\n";
	print APBS_FEW_IN "MMPBSAexe\t".$c_in{'mmpbsa_pl'}."\n";	# Location of MM-PBSA executable
	print APBS_FEW_IN "MMPBSAstat\t"."\$AMBERHOME/bin/mm_pbsa_statistics.pl\n"; # Executable for statics analysis
	print APBS_FEW_IN "APBSexe\t".$c_in{'apbs_exe'}."\n";		# Location of APBS executable
	print APBS_FEW_IN "DRAWmem\t".$currentFEW_dir."/miscellaneous/draw_membrane2\n";	# Location of program for
																					# membrane incorporation
	print APBS_FEW_IN "CALCdir\t".$calc_dir."\n";		# Current calculation directory
	print APBS_FEW_IN "STRUCT\t".$struct."\n";			# Name of current ligand structure
	print APBS_FEW_IN "START\t".$start."\n";			# First snapshot to regard
	print APBS_FEW_IN "STOP\t".$stop."\n";				# Last snapshot to regard
	print APBS_FEW_IN "OFFSET\t".$offset."\n";			# Offset between snapshots
	print APBS_FEW_IN "PDIE\t$c_in{'indi'}\n";			# Protein dielectric
	print APBS_FEW_IN "SDIE\t80\n"; 					# Solvent dielectric
	print APBS_FEW_IN "ZMEM\t$c_in{'zmem'}\n";			# Lower boundary of membrane slab
	print APBS_FEW_IN "TMEM\t$c_in{'tmem'}\n";			# Membrane thickness
	print APBS_FEW_IN "DIELC_MEM\t$c_in{'dielc_mem'}\n";		# Dielectric constant of membrane
	print APBS_FEW_IN "T_SEC_SLAB\t$c_in{'t_sec_slab'}\n";		# Thickness of second slab
	print APBS_FEW_IN "DIELC_SEC_SLAB\t$c_in{'dielc_sec_slab'}\n";	# Dielectric constant of the second slab
	print APBS_FEW_IN "T_THIRD_SLAB\t$c_in{'t_third_slab'}\n";		# Thickness of third slab
	print APBS_FEW_IN "DIELC_THIRD_SLAB\t$c_in{'dielc_third_slab'}\n";	# Dielectric constant of third slab
	print APBS_FEW_IN "IONCONC\t$c_in{'ionconc'}\n";	# Salt concentration [M]
	print APBS_FEW_IN "R_TOP\t$c_in{'r_top'} \n";		# Upper membrane exclusion radius
	print APBS_FEW_IN "R_BOTTOM\t$c_in{'r_bottom'}\n";	# Lower membrane exclusion radius
	
	if($c_in{'focus'} == 0){
		print APBS_FEW_IN "FOCUS\t0\n";
		if(( exists $c_in{'glen_s'})&&($c_in{'glen_s'} != 0)){
			print APBS_FEW_IN "GLEN_S\t".$c_in{'glen_s'}."\n";
		}
		else{
			print APBS_FEW_IN "GLEN_S\t100\n";			# Default grid length
		}
		if(($c_in{'dime'})&&($c_in{'dime'} != 0)){
			print APBS_FEW_IN "DIME\t".$c_in{'dime'}."\n";
		}
		else{
			print APBS_FEW_IN "DIME\t161\n";			# Dimensions / No. of grid points
		}
		#START LW_03_2018
		if($c_in{'prot_prot'} == 1){
			print APBS_FEW_IN "SOLV\t".$currentFEW_dir."/miscellaneous/apbs_mem_solv_prot.in\n";
			print APBS_FEW_IN "DUMMY\t".$currentFEW_dir."/miscellaneous/apbs_mem_dummy_prot.in\n";	# TEMPLATE for generation of maps
		}
		else{
			print APBS_FEW_IN "SOLV\t".$currentFEW_dir."/miscellaneous/apbs_mem_solv.in\n";	# TEMPLATE for solvation free energy calculations
			print APBS_FEW_IN "DUMMY\t".$currentFEW_dir."/miscellaneous/apbs_mem_dummy.in\n";	# TEMPLATE for generation of maps
		}
		#END LW_03_2018
	}
	else{
		print APBS_FEW_IN "FOCUS\t1\n";
		if((( exists $c_in{'glen_l'})&&($c_in{'glen_l'} != 0))
		 &&(( exists $c_in{'glen_m'})&&($c_in{'glen_m'} != 0))
		 &&(( exists $c_in{'glen_s'})&&($c_in{'glen_s'} != 0))){
			print APBS_FEW_IN "GLEN_L\t$c_in{'glen_l'}\n";		# Large grid length
			print APBS_FEW_IN "GLEN_M\t$c_in{'glen_m'}\n";		# Medium grid length
			print APBS_FEW_IN "GLEN_S\t$c_in{'glen_s'}\n";		# Small grid length
		}
		else{
			print APBS_FEW_IN "GLEN_L\t300\n";
			print APBS_FEW_IN "GLEN_M\t200\n";
			print APBS_FEW_IN "GLEN_S\t100\n";
		}
		#START LW_03_2018
		if($c_in{'prot_prot'} == 1){
			print APBS_FEW_IN "SOLV\t".$currentFEW_dir."/miscellaneous/apbs_mem_solv_focus_prot.in\n";
			print APBS_FEW_IN "DUMMY\t".$currentFEW_dir."/miscellaneous/apbs_mem_dummy_focus_prot.in\n";	# TEMPLATE for generation of maps
		}
		else{
			print APBS_FEW_IN "SOLV\t".$currentFEW_dir."/miscellaneous/apbs_mem_solv_focus.in\n";	# TEMPLATE for solvation free energy calculations
			print APBS_FEW_IN "DUMMY\t".$currentFEW_dir."/miscellaneous/apbs_mem_dummy_focus.in\n";	# TEMPLATE for generation of maps
		}
		#END LW_03_2018

		if((exists $c_in{'dime'})&&($c_in{'dime'} != 0)){
			print APBS_FEW_IN "DIME\t".$c_in{'dime'}."\n";
		}
		else{
			print APBS_FEW_IN "DIME\t97\n";			# Default dimensions / No. of grid points
		}
	}
	
	if($c_in{'Decomp'} == 1){
		print APBS_FEW_IN "DEC\t1\n";
	}
	#START LW_03_2018
	if($c_in{'prot_prot'} == 1){
		print APBS_FEW_IN "PROT_PROT\t1\n";
	}
	if($c_in{'npsolv'} == 1){
		print APBS_FEW_IN "NPSOLV\t1\n";
	}
	#END LW_03_2018
	print APBS_FEW_IN "VERBOSE\t0\n";
	close APBS_FEW_IN;
}

# Generate MMPBSA-input script
sub gen_pbs_script{

	my $calc_dir = shift;
	my $mmpbsa_in = shift;
	my $rel_path = shift;
	my $struct = shift;
	my $job_name = shift;
	my $rel_snaps = shift;
	my $rel_pqr_dir = shift;
	my $rel_topo = shift;
	my $currentFEWdir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	open(PBS_TEMPL, $c_in{'pbs_mmpbsa'}) || die "Cannot open ".$c_in{'pbs_mmpbsa'}.".\n";
	open(PBS_FILE, ">$calc_dir/run_mmpbsa.pbs") || die "Cannot open $calc_dir/run_mmpbsa.pbs for writing.\n";

	while(my $pbs_l = <PBS_TEMPL>){
		chomp($pbs_l);

		if($pbs_l =~ m/ \-N/){
			$pbs_l = $pbs_l." ".$job_name;
		}
		if($pbs_l eq "#PBS -o"){
			$pbs_l = $pbs_l." ".$job_name.".OU";
		}
		if($pbs_l eq "#PBS -e"){
			$pbs_l = $pbs_l." ".$job_name.".ER";
		}
		#START - For SLUM script support
		if($pbs_l =~ m/ --job-name=/){
			$pbs_l = $pbs_l.$job_name;
		}
		if($pbs_l =~ m/ --error=/){
			$pbs_l = $pbs_l.$job_name.".ER";
		}
		if($pbs_l =~ m/ --output=/){
			$pbs_l = $pbs_l.$job_name.".OU";
		}
		#END - For SLUM script support
		if($pbs_l eq "set SY="){
			$pbs_l = $pbs_l.$struct;
		}
		if($pbs_l eq "set TOPO="){
			$pbs_l = $pbs_l.$c_in{'pbs_m_path'}."/".$rel_topo;
		}
		if($pbs_l eq "set SNAPS="){
			$pbs_l = $pbs_l.$c_in{'pbs_m_path'}."/".$rel_snaps;
			if($c_in{'ImplMem'} == 1){
				$pbs_l .= "\nset PQR_DIR=".$c_in{'pbs_m_path'}."/".$rel_pqr_dir;
			}
		}
		if($pbs_l =~ m/set CALC/){
			$pbs_l = $pbs_l.$job_name."_pb".$c_in{'PB'}."_gb".$c_in{'GB'};
		}
		if($pbs_l eq "set SCRIPT="){
			$pbs_l = $pbs_l.$c_in{'pbs_m_path'}."/".$rel_path."/run_mmpbsa.pbs";
		}
		if(($pbs_l eq "set PATH=")||($pbs_l eq "set WORKDIR=")){
			$pbs_l = $pbs_l.$c_in{'pbs_m_path'}."/".$rel_path;
		}
		if(($pbs_l eq "set PARAMS=")&&($c_in{'ImplMem'} == 0)){
			$pbs_l = $pbs_l."mmpbsa.in";
		}
		if(($pbs_l eq "set PARAMS=")&&($c_in{'ImplMem'} == 1)){
			$pbs_l = $pbs_l."mmpbsa_FEWmem.in";
		}
		if(($pbs_l eq "set MMPBSA=")&&($c_in{'ImplMem'} == 0)){
			$pbs_l = $pbs_l.$c_in{'mmpbsa_pl'};
		}
		if(($pbs_l eq "set MMPBSA=")&&($c_in{'ImplMem'} == 1)){
			if($c_in{'focus'} == 0){
				$pbs_l = $pbs_l.$currentFEWdir."/miscellaneous/mmpbsa_FEWmem.pl";
			}
			else{
				$pbs_l = $pbs_l.$currentFEWdir."/miscellaneous/mmpbsa_FEWmem_focus.pl";
			}
		}
		if(($pbs_l eq "ln -s \$SNAPS snapshots")&&($c_in{'ImplMem'} == 1)){
			$pbs_l = $pbs_l."\nln -s \$PQR_DIR pqr_snaps";
		}

		print PBS_FILE $pbs_l."\n";
	}
	close PBS_TEMPL;
	close PBS_FILE;
}

1;
