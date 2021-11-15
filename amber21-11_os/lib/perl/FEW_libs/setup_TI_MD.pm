# Subroutines needed for setup of the TI simulations
use prepare_top_crd_MD;
use Normality 'shapiro_wilk_test';
use strict;

# Modifying the blocking groups in order to enable backwards compartibility
# to AMBER versions < 18.
sub modify_blocking_groups{
	my $org_file = shift;
	my $out_str = "";

	open(PDB, $org_file) || die "Cannot open pdb-file $org_file.\n";
	while(my $pdb_line = <PDB>){
		chomp($pdb_line);
		if(($pdb_line =~ m/ATOM/)||($pdb_line =~ m/HETATM/)){
			my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z) =
			unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8", $pdb_line);
			$resname =~ s/^\s+//g;
			$resname =~ s/\s+$//g;

			# Don't use ACE/NME atoms not recognized by leap
			if((($resname eq "ACE")&&($atom =~ m/HH3/))||(($resname eq "NME")&&($atom =~ m/^H/))){
			}
			else{
				$out_str .= $pdb_line."\n";
			}
		}
		else{
			$out_str .= $pdb_line."\n";
		}
	}
	close PDB;
	unlink($org_file);

        open(PDB_MOD, ">$org_file") || die "Cannot open pdb-file $org_file for writing.\n";
	print PDB_MOD $out_str;
        close PDB_MOD;
}

# Inserting TER between chains of PDB file according to user definition
sub insert_ter{
	my $pdb_com_tmp = shift;
	my $pdb_com_start = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $ter_site_found = 0;
	my $ter_site = 0;
	my $out_str = "";
	my @TER_sites = split(/\,/, $c_in{'chain_termini'});
	
	open(PDB_IN, $pdb_com_tmp) || die "Cannot open PDB file $pdb_com_tmp for reading\n";
	while(my $pdb_in = <PDB_IN>){
		chomp($pdb_in);
		my @pdb_in = split(/\s+/, $pdb_in);
		my $next_res = $ter_site + 1;
		if(($ter_site_found == 1)&&($pdb_in[4] == $next_res)){
			$pdb_in = "TER\n$pdb_in";
			$ter_site_found = 0;
		}
		if(grep {$pdb_in[4] == $_} @TER_sites){
			$ter_site_found = 1;
			$ter_site = $pdb_in[4];
		}
		$out_str .= $pdb_in."\n";
	}
	close PDB_IN;
	
	open(PDB_OUT, ">$pdb_com_start") || die "Cannot open PDB file $pdb_com_start for writing\n";
	print PDB_OUT $out_str;
	close PDB_OUT;
}


# Change ligand residue name in PDB-files
sub change_res_name{
	my $file = shift;
	my $o_dir = shift;
	my $alias = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	open(ORG_FILE, $file) || die "Cannot open file $file.\n";
	open(TMP_FILE, ">$o_dir/tmp.pdb") || "Cannot open temporary file $o_dir/tmp.pdb for writing.\n";
	
	while(my $org_line = <ORG_FILE>){
		chomp($org_line);
		
		$org_line =~ s/<1>/$alias/;
		$org_line =~ s/MOL/$alias/;
		$org_line =~ s/UNN/$alias/;
		
		if($alias eq $c_in{'Lalias_end'}){
			$org_line =~ s/$c_in{'Lalias_start'}/$c_in{'Lalias_end'}/;
		}
		
		print TMP_FILE $org_line."\n";
	}
	close ORG_FILE;
	close TMP_FILE;
	
	unlink($file);
	copy("$o_dir/tmp.pdb", $file);
	unlink("$o_dir/tmp.pdb");
}


# Remove superfuous atoms from PDB-files of state V1
sub del_atoms{
	my $file = shift;
	my $o_dir = shift;
	my $atoms_to_del = shift;
	my $c_in = shift;
	my %c_in = %{$c_in};
	
	my @atoms = split(/\,/, $atoms_to_del);
	
	open(PDB, $file) || die "Cannot open pdb-file $file.\n";
	open(TMP, ">$o_dir/tmp.pdb") || die "Cannot open temporary file $o_dir/tmp.pdb for writing.\n";
	
	while(my $pdb_l = <PDB>){
		chomp($pdb_l);
		
		if($pdb_l =~ m/$c_in{'Lalias_end'}/){
			my @pdb_l = split(/\s+/, $pdb_l);
			if(grep{$_ eq $pdb_l[2]}@atoms){
				next;
			} 
			else{
				print TMP "$pdb_l\n";
			}
		}
		else{
			print TMP "$pdb_l\n";
		}
	}
	close PDB;
	close TMP;
	
	unlink($file);
	copy("$o_dir/tmp.pdb", $file);
	unlink("$o_dir/tmp.pdb");	
}


# Shift position of atoms in mask to end of mol2-file to avoid mismatches of atoms
# between state V0 and state V1
sub shift_position{
	my $org_mol2_file = shift;
	my $mask = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	
	my @mask = split(/\@/, $c_in{$mask});
	my @atoms_in_mask = split(/\,/, @mask[1]);
	my $state = $mask[0];
	
	# Case: State contains nothing since soft core in one of the two states is equal to ''
	my $file_name = "";
	if($state eq ""){
		my @path_components = split(/\//, $org_mol2_file);
		my @name = split(/\./, $path_components[$#path_components]);
		$file_name = $name[0];
		
		if($file_name eq $c_in{'Lname_start'}){
			$state = $c_in{'Lalias_start'};
		}
		if($file_name eq $c_in{'Lname_end'}){
			$state = $c_in{'Lalias_end'};
		} 
	}
	
	# Generate PDB-file from mol2-file with antechamber
	my $q_pdb_file = substr($org_mol2_file, 0, -4);
	$q_pdb_file .= "pdb";
	my $create_pdb_log = $c_in{'root_path'}."/antechamber_create_pdb.log";
	my $c_create_pdb = "antechamber -i $org_mol2_file -fi mol2 -o $q_pdb_file -fo pdb";
	$c_create_pdb .= " > ".$create_pdb_log." 2>&1";
	call_prog($c_create_pdb, $create_pdb_log);
	unlink($create_pdb_log);
	
	# Read Position and re-name information
	my %position;
	my %rename;
	my @atom_names = ();
	my @position;
	
	if(-e $c_in{'match_list'}){
		open(MATCH, $c_in{'match_list'}) || die "Cannot open file ".$c_in{'match_list'}." for reading.\n";
		while(my $match = <MATCH>){
			my $check_comment = substr($match, 0, 1);
			if($check_comment ne "#"){
				chomp($match);
				my @match = split("\t", $match);
				$position{$match[0]} = $match[2];
				$rename{$match[2]} = $match[1];
				push(@atom_names, $match[2]);
				push(@position, $match[0]);
			}
		}
	}
	else{
		print "The match-file you specified does not exist. Please make sure that you specified\n";
		print "the correct path and name of the match-file.\n";
		exit;
	}
	
	my @sorted_position = sort {$a <=> $b} @position;
	
	# Read PDB-file and store information about atoms separately for later rearragement
	my %str_out;
	my $str_out = "";
	my $mask_atom_str = "";
	
	open(Q_PDB, $q_pdb_file) || die "Cannot open file $q_pdb_file for reading.\n";
	while(my $q_pdb = <Q_PDB>){
		chomp($q_pdb);
		if($q_pdb =~ m/^ATOM/){
			my @q_pdb = split(/\s+/, $q_pdb);
			if(grep {$_ eq $q_pdb[2]} @atoms_in_mask){
				$mask_atom_str .= $q_pdb."\n";
			}
			else{
				if($state eq $c_in{'Lalias_end'}){
					$str_out{$q_pdb[2]} = $q_pdb."\n";
					if(grep {$_ eq $q_pdb[2]} @atom_names)
					{}
					else{
						print "At least one of the atom names you specified for state V1 in your matching-list\n";
						print "was not found in the mol2-file you provided for state V1. Please ensure that the\n";
						print "list contains all none soft core atoms of V1 and correct the naming.\n";
						print "Look for atom $q_pdb[2] in your matching list.\n";
						exit; 
					}
				}

				if($state eq $c_in{'Lalias_start'}){
					$str_out .= $q_pdb."\n";
				}
			}
		}
	}
	close Q_PDB;
	
	if($state eq $c_in{'Lalias_end'}){
		foreach my $p (@sorted_position){
			$str_out .= $str_out{$position{$p}};
		}
	}
	
	# Create new PDB-file
	my $t_pdb_file = substr($org_mol2_file, 0, -5);
	$t_pdb_file .= "_shifted.pdb";
	open(T_PDB, ">$t_pdb_file") || die "Cannot open file $t_pdb_file for writing.\n";
	print T_PDB $str_out;
	print T_PDB $mask_atom_str;
	close T_PDB;
	
	# Convert modified PDB-file into mol2-file
	my $t_mol2_file = substr($t_pdb_file, 0, -3);
	$t_mol2_file .= "mol2";
	my $convert_mod_log = $c_in{'root_path'}."/antechamber_convert_mol2.log";
	my $c_convert_mod = "antechamber -i $t_pdb_file -fi pdb -o $t_mol2_file -fo mol2";
	$c_convert_mod .= " > ".$convert_mod_log." 2>&1";
	call_prog($c_convert_mod, $convert_mod_log);
	unlink($convert_mod_log);
	
	# Read charges, atomtypes from original mol2-file
	my %charges;
	my %atom_types;
	my $start_reading = 0;
	my $stop_reading = 0;
	
	open(Q_MOL, $org_mol2_file) || die "Cannot open file $org_mol2_file for reading.\n";
	while(my $q_mol = <Q_MOL>){
		chomp($q_mol);
		if($q_mol eq "@<TRIPOS>BOND"){
			$stop_reading = 1;
		}
		if(($start_reading == 1)&&($stop_reading == 0)){
			my @q_mol = split(/\s+/, $q_mol);
			$atom_types{$q_mol[2]} = $q_mol[6];
			$charges{$q_mol[2]} = $q_mol[9];
		}
		if($q_mol eq "@<TRIPOS>ATOM"){
			$start_reading = 1;
		}
	}
	close Q_MOL;
	
	# Exchange atomtypes, charges, and atom names in new mol2-file
	$start_reading = 0;
	$stop_reading = 0;
	my $new_file_str = "";
	
	open(T_MOL, $t_mol2_file) || die "Cannot open file $t_mol2_file for reading.\n";
	while(my $t_mol = <T_MOL>){
		chomp($t_mol);
		if($t_mol =~ m/No Charge or Current Charge/){
			$t_mol =~ s/No Charge or Current Charge/resp/;
		}
		if($t_mol eq "@<TRIPOS>BOND"){
			$stop_reading = 1;
		}
		if(($start_reading == 1)&&($stop_reading == 0)){
			my @t_mol = split(/\s+/, $t_mol);
			
			# Replace atom types and atom names
			my @modes = ();
			if($state eq $c_in{'Lalias_end'}){
				@modes = ("atomtypes", "atomnames");
			}
			if($state eq $c_in{'Lalias_start'}){
				@modes = ("atomtypes");
			}
			
			for my $mode (@modes){
				my $length_t = 0;
				my $length_q = 0;
				my $str_t = "";
				my $str_q = "";
				
				if($mode eq "atomtypes"){
					$str_t = $t_mol[6];
					$str_q = $atom_types{$t_mol[2]};
				}
				if($mode eq "atomnames"){
					$str_t = $t_mol[2];
					$str_q = $rename{$t_mol[2]};
				}	
					
				$length_t = length($str_t);
				$length_q = length($str_q);	
				
				if($length_t == $length_q){
					$t_mol =~ s/$str_t/$str_q/;
				}
				if((($length_t == 2)&&($length_q == 1))||(($length_t == 3)&&($length_q == 2))){
					$t_mol =~ s/$str_t/$str_q /;
				}
				if((($length_t == 1)&&($length_q == 2))||(($length_t == 2)&&($length_q == 3))){
					$t_mol =~ s/$str_t /$str_q/;
				}
				if(($length_t == 1)&&($length_q == 3)){
					$t_mol =~ s/$str_t  /$str_q/;
				}
				if(($length_t == 3)&&($length_q == 1)){
					$t_mol =~ s/$str_t/$str_q  /;
				}
			}
			
			# Replace charges
			if($charges{$t_mol[2]} < 0){
				$t_mol =~ s/ $t_mol[9]/$charges{$t_mol[2]}/;
			}
			if($charges{$t_mol[2]} > 0){
				$t_mol =~ s/$t_mol[9]/$charges{$t_mol[2]}/;
			}
		}
		if($t_mol eq "@<TRIPOS>ATOM"){
			$start_reading = 1;
		}
		$new_file_str .= $t_mol."\n";
	} 
	close T_MOL;
	
	open(MOL_FINAL, ">$org_mol2_file") || die "Cannot open file $org_mol2_file for writing\t";
	print MOL_FINAL $new_file_str;
	close MOL_FINAL;
	
	unlink($t_mol2_file);
	unlink($t_pdb_file);
}


# Generate leap-input file for Library creation
sub setup_leap_for_lib_TI{
	my $lib_leap_in = shift;
	my $lib_file = shift;
	my $frcmod_file = shift;
	my $alias = shift;
	my $setup_dir = shift;
	my $end_mol2 = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};

	open(LEAP_LIB, ">$lib_leap_in") || die "Cannot open leap input file $lib_leap_in.\n"; 

	if($c_in{'backwards'} == 0){
		print LEAP_LIB "source leaprc.protein.ff14SB\n";
		if($c_in{'water_model'} eq "TIP3P"){
			print LEAP_LIB "source leaprc.water.tip3p\n";
		}
		if($c_in{'water_model'} eq "OPC"){
			print LEAP_LIB "source leaprc.water.opc\n";
		}
	}
	else{
		print LEAP_LIB "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff99SB\n";
	        print LEAP_LIB "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";
	}
	print LEAP_LIB "source leaprc.gaff\n";
	print LEAP_LIB "loadAmberParams ".$frcmod_file."\n";
	if($alias eq $c_in{'Lalias_end'}){
		print LEAP_LIB $alias." = loadmol2 ".$end_mol2."_amb\n";
	}
	else{
		print LEAP_LIB $alias." = loadmol2 ".$setup_dir."/".$c_in{'Lname_start'}.".mol2\n";
	}
	print LEAP_LIB "check ".$alias."\n";
	print LEAP_LIB "charge ".$alias."\n";
	print LEAP_LIB "saveoff ".$alias." ".$lib_file."\n";
	print LEAP_LIB "quit\n";
	close LEAP_LIB;
}


# Check charge and warnings
sub check_warn_chrg_leap{
	my $f_log_file = shift;
	my $library_file = shift;
	my $warn = 0;
	my $charge;
				
	open(LOG, $f_log_file) || die "Cannot open log-file $f_log_file\n";
				
	while(my $log_line = <LOG>){
		chomp($log_line);
					
		if(($log_line =~ m/Warning/)||($log_line =~ m/WARNING/)){
			if($log_line =~ m/perturbed charge/){
				next;
			}
                        elsif($log_line =~ m/Exiting LEaP: Errors/){
                                next;
                        }
			else{
				$warn = 1;
			}
		}
		
		if($log_line =~ m/Total unperturbed charge/){
			my @c = split(/\:\s+/, $log_line);
			$charge = $c[1];
		}
	}
			
	if($warn != 1){
		print "Terminated LEaP without warnings.\n";
	}
	if($warn == 1){
		print "LEaP run terminated with warnings. Please inspect the LEaP log file at ".$f_log_file."\n";
	}


	# Modify library of ligand if charge is marginally smaller than next integer value
	my $int = int($charge);
	my $res;
	
	if($charge > 0){
		$res = $charge - $int;
	}
	if($charge < 0){
		$res = $charge + $int;
	} 
	my $diff = 0;
			
	# Case charge slightly smaller than next integer
	if(abs($res>0.9)){
		if($charge>0){
			$diff = (($int+1)-$charge);
		}
		if($charge<0){
			$diff = (($int-1)+$charge);
		}
	}

	# Case charge slightly larger than next integer
	if(abs($res)<0.005){
		if($charge<0){
			$diff = abs($charge);
		}
		if($charge>0){
			$diff = (-1)*$diff;
		}
	}
			
	if((abs($diff)<0.005)&&($diff != 0)){
		$library_file = &change_charge_lib($library_file, $diff);
	}
}


# Generate leap input file for setup of coordinate 
# and topology files for start- and end-structures
sub gen_leap_crd_top_in{
	my $o_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $cys_bridges_ref = shift;
	my %cys_bridges;
	my $prog = shift;
	if($c_in{'cys_f'} ne ""){
		%cys_bridges = %{$cys_bridges_ref};
	}
	my $frcmod_start_file = $c_in{'root_path'}."/leap/".$c_in{'Lname_start'}."/".$c_in{'Lname_start'}.".frcmod";

	my $leap_file_name = "";
	if($prog eq "sander"){
		$leap_file_name = "$o_dir/leap_top_crd.in";
	}
	if($prog eq "pmemd"){
		$leap_file_name = "$o_dir/leap_top_crd_pmemd.in";
	}

	open(LEAP_CT, ">".$leap_file_name) || die "Cannot open leap_top_crd.in for writing.\n";
	if($c_in{'backwards'} == 0){
		print LEAP_CT "source leaprc.protein.ff14SB\n";
		if($c_in{'water_model'} eq "TIP3P"){
			print LEAP_CT "source leaprc.water.tip3p\n";
		}
		if($c_in{'water_model'} eq "OPC"){
			print LEAP_CT "source leaprc.water.opc\n";
		}
	}
	else{
		print LEAP_CT "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff99SB\n";
	        print LEAP_CT "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";

	}
	print LEAP_CT "source leaprc.gaff\n";

	# Load libraries for start and end ligands
	if(-e $o_dir."/".$c_in{'Lname_start'}."mod.lib"){
		print LEAP_CT "loadoff ".$o_dir."/".$c_in{'Lname_start'}."mod.lib\n";
	}
	else{
		print LEAP_CT "loadoff ".$o_dir."/".$c_in{'Lname_start'}.".lib\n";
	}
	
	if(-e $o_dir."/".$c_in{'Lname_end'}."_mod.lib"){
		print LEAP_CT "loadoff ".$o_dir."/".$c_in{'Lname_end'}."_mod.lib\n";
	}
	else{
		print LEAP_CT "loadoff ".$o_dir."/".$c_in{'Lname_end'}.".lib\n";
	}
	
	if(($c_in{'add_lib'})&&($c_in{'add_lib'} ne "")&&(-e $c_in{'add_lib'})){
		print LEAP_CT "loadoff ".$c_in{'add_lib'}."\n";
	}
	if(($c_in{'add_frcmod'})&&($c_in{'add_frcmod'} ne "")&&(-e $c_in{'add_frcmod'})){
		print LEAP_CT "loadAmberParams ".$c_in{'add_frcmod'}."\n";
	}

	# Load parameters for start and end ligands	
	for my $name ($c_in{'Lname_start'}, $c_in{'Lname_end'}){
		for my $tag ("_com", "_lig") {
			if($name.$tag eq $c_in{'Lname_start'}."_com"){
				print LEAP_CT "\n";
				print LEAP_CT "loadAmberParams ".$frcmod_start_file."\n";
			}
			if($name.$tag eq $c_in{'Lname_end'}."_com"){
				print LEAP_CT "\n";
				print LEAP_CT "loadAmberParams ".$o_dir."/".$name.".frcmod\n";
			}
		}
	}

	if($prog eq "pmemd"){
		my $combined = $c_in{'Lname_start'}."_".$c_in{'Lname_end'};
		for my $tag ("_com", "_lig") {
			print LEAP_CT $combined.$tag." = loadpdb ".$o_dir."/".$combined.$tag."_TIin_leap.pdb\n";
		
			# Include disulfide bridges
			if($tag eq "_com"){
				if($c_in{'cys_f'} ne ""){
					foreach my $k (keys %cys_bridges){
						print LEAP_CT "bond ".$combined.$tag.".".($k+2).".SG ".$combined.$tag.".".($cys_bridges{$k}+2).".SG\n";
					}
				}
			}

			print LEAP_CT "setBox ".$combined.$tag." vdw\n";
			print LEAP_CT "saveamberparm ".$combined.$tag." ".$o_dir."/".$combined.$tag."_PARMEDin.top ".$o_dir."/".$combined.$tag."_PARMEDin.crd\n";
			print LEAP_CT "savepdb ".$combined.$tag." ".$o_dir."/".$combined.$tag."_PARMEDin.pdb\n";
			
			if($combined.$tag eq $combined."_com"){
				print LEAP_CT "\n";
			}
		}
	}
	if($prog eq "sander"){
		for my $name ($c_in{'Lname_start'}, $c_in{'Lname_end'}){
			for my $tag ("_com", "_lig") {
				print LEAP_CT $name.$tag." = loadpdb ".$o_dir."/".$name.$tag.".pdb\n";
			
				# Modify to read automatic from input template later
				if($tag eq "_com"){
					if($c_in{'cys_f'} ne ""){
						foreach my $k (keys %cys_bridges){
							print LEAP_CT "bond ".$name.$tag.".".$k.".SG ".$name.$tag.".".$cys_bridges{$k}.".SG\n";
						}
					}
				}	
			
				print LEAP_CT "setBox ".$name.$tag." vdw\n";
				print LEAP_CT "saveamberparm ".$name.$tag." ".$o_dir."/".$name.$tag."_TIin.top ".$o_dir."/".$name.$tag."_TIin.crd\n";
				print LEAP_CT "savepdb ".$name.$tag." ".$o_dir."/".$name.$tag."_TIin.pdb\n";
			
				if(($name.$tag eq $c_in{'Lname_end'}."_com")||($name.$tag eq $c_in{'Lname_start'}."_com")){
					print LEAP_CT "\n";
				}
			}
		}
	}
	print LEAP_CT "quit\n";
	close LEAP_CT;
}


# Generation of input file for ParmED for TI setup with pmemd
sub generate_parmed_input{
	my $o_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $parmed_input_file_name = $o_dir."/parmed.in";

	open(PARMED, ">".$parmed_input_file_name) || die "Cannot open ".$o_dir."/parmed.in for writing.\n";
	for my $tag ("com", "lig"){
		print PARMED "parm ".$o_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_PARMEDin.top\n";
		print PARMED "loadCoordinates ".$o_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_PARMEDin.crd\n";
		print PARMED "timerge :1 :2 :1 :2 tol 0.5\n";
		print PARMED "outparm ".$o_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.top ".$o_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.crd\n";
		print PARMED "outPDB ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.pdb\n";
	}
}


# Generation of PDB files for setup of TI simulations with PMEMD
sub gen_PDBs_pmemd{
	my $o_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};


	for my $tag ("_com", "_lig"){

		# Read start ligand information and common information
		my $start_lig_str = "";
		my $common_info_str = "";
		my $start_res = 0;
		
		open(PDB_START, $o_dir."/".$c_in{'Lname_start'}.$tag."_TIin.pdb");

		while(my $pdb_start_line = <PDB_START>){
			chomp($pdb_start_line);

			if($pdb_start_line =~ m/CRYST/){
				next;
			}

			# Handle TER after start residue to avoid double TER entry
			if($pdb_start_line =~ m/TER/){
				if($start_res == 1){
					$start_lig_str .= "TER\n";
					$start_res = 0;
					next;
				}
				else{
					$common_info_str .= "TER\n";
					next;
				}
			}
			else{
				my $resname = substr($pdb_start_line, 17, 3);
				# Save string of start ligand			
				if($resname eq $c_in{'Lalias_start'}){
					$start_lig_str .= $pdb_start_line."\n";
					$start_res = 1;
					next;
				}
				# Save common string
				else{
					$common_info_str .= $pdb_start_line."\n";
					next;
				}
			}
		}
		close(PDB_START);

		# Read end ligand information
		my $end_lig_str = "";

		open(PDB_END, $c_in{'Lname_end'}.$tag."_TIin.pdb");

		while(my $pdb_end_line = <PDB_END>){
			chomp($pdb_end_line);
		
			if($pdb_end_line =~ m/CRYST/){
				next;
			}
	
			my $resname = substr($pdb_end_line, 17, 3);

			# Save string of start ligand			
			if($resname eq $c_in{'Lalias_end'}){
				$end_lig_str .= $pdb_end_line."\n";
			}
		}
		close(PDB_END);

		# Add final TER card
		$end_lig_str .= "TER\n";

		# Write new combined PDB file
		open(COMBINED_PDB, ">".$o_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}.$tag."_TIin.pdb");
		print COMBINED_PDB $start_lig_str;
		print COMBINED_PDB $end_lig_str;
		print COMBINED_PDB $common_info_str;
		close(COMBINED_PDB);
	}
}


# Replace BoxInfo in topology files according to Box-Info in initial
# RESTRT-files -> Box size after equilibration
# This is required since the box info included by running xleap with 
# the "setBox" command is not equivalent to the box size of the system 
# after equilibration, since the system is not back-imaged into the
# initial box.
sub setBox_info{
	my $o_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $prog = shift;

	my @prefix;
	if($prog eq "sander"){
		@prefix = ($c_in{'Lname_start'}, $c_in{'Lname_end'});
	}
	if($prog eq "pmemd"){
		@prefix = ($c_in{'Lname_start'}."_".$c_in{'Lname_end'});
	}
	for my $name (@prefix){
		for my $tag ("lig", "com"){;
			my $top_file = $o_dir."/".$name."_".$tag."_TIin.top";
			my $crd_file = $o_dir."/".$name."_".$tag."_TIin.crd";
		
			my $rst_file;
			if($tag eq "lig"){
				$rst_file = $c_in{'rst_lig'};
			}
			else{
				$rst_file = $c_in{'rst_com'};
			}
		
			# Replace box info
			my $rst = File::ReadBackwards->new($rst_file) || "Cannot open file $rst_file for reading.\n";
			my $rst_line = $rst->readline;
			my @rst_line = split(/\s+/, $rst_line);
			my $x = sprintf("%.8e", $rst_line[1]);
			my $y = sprintf("%.8e", $rst_line[2]);
			my $z = sprintf("%.8e", $rst_line[3]);
			$x =~ s/e/E/g;
			$y =~ s/e/E/g;
			$z =~ s/e/E/g;
		
			# In coordinate file
			my $last_line;
			my $crd_out_str = "";
			open(CRD, $crd_file) || die "Cannot open file $crd_file for reading.\n";
			while(my $crd = <CRD>){
				chomp($crd);
				if(eof){
					$crd_out_str .= $rst_line;
				}
				else{
					$crd_out_str .= $crd."\n";
				}
			}
			close CRD;
			
			open(CRD_NEW, ">$crd_file") || die "Cannot open file $crd_file for writing.\n";
			print CRD_NEW $crd_out_str;
			close CRD_NEW;
			
			# In Topology file
			my $box_info = 0;
			my $top_out_str = "";
			open(TOP, $top_file) || die "Cannot open file $top_file for reading.\n";
			while(my $top = <TOP>){
				chomp($top);
				if($box_info == 2){
					my @top = split(/\s+/, $top);
					$top = sprintf("%16s%16s%16s%16s", $top[1], $x, $y, $z);
					$box_info = 3;
				}
				if($box_info == 1){
					$box_info = 2;
				}
				if($top =~ m/BOX_DIMENSIONS/){
					$box_info = 1;
				}
				$top_out_str .= $top."\n";
			}
			close TOP;
			
			open(TOP_NEW, ">$top_file") || die "Cannot open file $top_file for writing.\n";
			print TOP_NEW $top_out_str;
			close TOP_NEW;
		}
	}
}


# Determine total number of residues in ligand and complex files
sub determine_no_of_res {

	my $name = shift;
	my $alias = shift;
	my $o_dir = shift;
	my $temp_res_no;
	my $lig_res_no;
	my $com_res_no;
	my $com_res_lig;
	
	for my $tag ("com", "lig"){
		my $file = $o_dir."/".$name."_".$tag."_TIin.pdb";
		
		open(PDB, $file) || die "Cannot open file $file\n";
		
		while(my $pdb_l = <PDB>){
			chomp($pdb_l);
			if($pdb_l =~ m/^ATOM/){
				my @pdb_l = split(/\s+/, $pdb_l);
				
				# Determine residue number of ligand residue in complex
				if(($pdb_l[3] eq $alias)&&($tag eq "com")){
					$com_res_lig = $pdb_l[4];
				}
				if($pdb_l[3] eq $alias){
					$temp_res_no = $pdb_l[4];
				}
			}
		}
		close PDB;
		
		if($tag eq "com"){
			$com_res_no = $temp_res_no;
		}
		else{
			$lig_res_no = $temp_res_no;
		}
	}
	return ($com_res_no, $lig_res_no, $com_res_lig);
}

# Determines new softcore mask for V1 after automatic renaming by leap
# and returns the new mask
sub determine_new_sc_str{
	my $setup_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	# Determine line in which soft core starts in ligand state V0
	my @org_mask0 = split(/\@/, $c_in{'mask0'});
	my @diff_atoms0 = split(/\,/, @org_mask0[1]);
	my $pdb_v0 = $setup_dir."/".$c_in{'Lname_start'}."_lig_TIin.pdb";
	my $sc_start_line_no;
	my $line_count = 0;
	
	open(PDB_VZERO, $pdb_v0) || die "Cannot open file $pdb_v0 for reading.\n";
	while(my $pdb_v0_line = <PDB_VZERO>){
		chomp($pdb_v0_line);
		
		if($pdb_v0_line =~ m/TER/){
			print "Error: The soft core you specified for state V0 was not found in PDB file of state V0.\n";
			exit;
		}
		
		if($pdb_v0_line =~ m/ATOM/){
			my @pdb = split(/\s+/, $pdb_v0_line);
			$line_count++;
			
			if(grep {$pdb[2] eq $_} @diff_atoms0){
				$sc_start_line_no = $line_count;
				last;
			}
		}
	} 

	# Read soft core for state V1 from ligand file V1
	my @org_mask1 = split(/\@/, $c_in{'mask1'});
	my @diff_atoms1 = split(/\,/, @org_mask1[1]);
	my %diff_atoms1;
	map { $diff_atoms1{$_} = 1 } @diff_atoms1;
	my $pdb_v1 = $setup_dir."/".$c_in{'Lname_end'}."_lig_TIin.pdb";
	my @new_sc1 = ();
	my $start_reading == 0;
	$line_count = 0;
	
	open(PDB_VONE, $pdb_v1) || die "Cannot open file $pdb_v1 for reading.\n";
	while(my $pdb_v1_line = <PDB_VONE>){
		chomp($pdb_v1_line);
		
		if($pdb_v1_line =~ m/TER/){
			last;
		}
		
		if($pdb_v1_line =~ m/ATOM/){
			my @pdb = split(/\s+/, $pdb_v1_line);
			$line_count++;
			
			if(($line_count == $sc_start_line_no)&&($diff_atoms1{$pdb[2]})&&(@diff_atoms1 >= 1)){
				$start_reading = 1;
			}
			if(($line_count == $sc_start_line_no)&&(! $diff_atoms1{$pdb[2]})&&(@diff_atoms1 <= 1)){
				$start_reading = 1;
			}
			if(($line_count == $sc_start_line_no)&&(! $diff_atoms1{$pdb[2]})&&(@diff_atoms1 > 1)){
				print "Error: Problem with soft core of state V1. The soft core of V1 was not found\n";
				print "at the expected position. Please make sure that the soft core specifications\n";
				print "for states V0 and V1 are correct.\n";
				exit;
			}
			
			if($start_reading == 1){
				push(@new_sc1, $pdb[2]);
			}
		}
	}
	
	# Generate appropriate format
	my $sc1 = join ",", map {$_} @new_sc1;
	my $new_sc1_str = $org_mask1[0]."@".$sc1;
	
	return $new_sc1_str;
}


# Setup equilibration files according to templates
sub equi_setup_TI{
	my $tag_dir = shift;
	my $tag = shift;
	my $lig_res_no = shift;
	my $com_res_no = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	# Parameters for equilibration: Currently a fixed equilibration procedure is 
	# implemented consisting of a single equilibration step according to provided
	# equilibration template followed by a 1 ns MD simulation that is considered
	# as part oft the equilibration phase 
	my $equi_step_no = 1; # Number of equilibration steps; currently one per default
	my $prod_time = 1; # Production time in [ns]
	my $requested_equi_time_templ = determine_requested_time("equilibration", $c_in{'equi_templ'}, \%c_in);
	my $requested_prod_time_templ = determine_requested_time("production", $c_in{'prod_templ'}, \%c_in);
	my $prod_no = $prod_time/$requested_prod_time_templ;
	my $first_lambda = 1;
	
	my @lambda = split(/\,/, $c_in{'equi_lambda'});
	

	if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		my $inpcrd_v0 = $c_in{'Lalias_start'}."_".$tag."_TIin.crd";
		my $inpcrd_v1 = $c_in{'Lalias_end'}."_".$tag."_TIin.crd";

		foreach my $lambda (@lambda){
	
			# Generate group-file for equilibration part
			for my $E (1 ... $equi_step_no){

				open(GROUP_EQUI, ">$tag_dir/group_equi0$E"."_l$lambda") || die "Cannot open file $tag_dir/group_equi0$E"."_l$lambda for writing\n";
				print GROUP_EQUI "-O -i $c_in{'Lalias_start'}"."_equi0$E"."_v0_l$lambda.in ". 
							"-o $c_in{'Lalias_start'}"."_equi0$E"."_v0_l$lambda.out ". 
							"-p $c_in{'Lalias_start'}"."_$tag"."_TIin.top ".
							"-c $inpcrd_v0 ".
							"-r $c_in{'Lalias_start'}"."_equi0$E"."_v0_l$lambda.rst ".
							"-x $c_in{'Lalias_start'}"."_equi0$E"."_v0_l$lambda.mdcrd\n";
								
				print GROUP_EQUI "-O -i $c_in{'Lalias_end'}"."_equi0$E"."_v1_l$lambda.in ". 
							"-o $c_in{'Lalias_end'}"."_equi0$E"."_v1_l$lambda.out ". 
							"-p $c_in{'Lalias_end'}"."_$tag"."_TIin.top ".
							"-c $inpcrd_v1 ".
							"-r $c_in{'Lalias_end'}"."_equi0$E"."_v1_l$lambda.rst ".
							"-x $c_in{'Lalias_end'}"."_equi0$E"."_v1_l$lambda.mdcrd\n";
				close GROUP_EQUI;
			
				$inpcrd_v0 = $c_in{'Lalias_start'}."_prod0$prod_no"."_v0_l$lambda".".rst";
				$inpcrd_v1 = $c_in{'Lalias_end'}."_prod0$prod_no"."_v1_l$lambda".".rst";
			}

			# Generate group-file for production part
			my $restrt = "equi01";
			for(my $p=1; $p<=$prod_no; $p++){
				group_file_prod_setup($tag_dir, $restrt, $tag, $p, $lambda, \%c_in);
				$restrt = "";		
			}
		}
	}

	
	# Generate equilibration file
	foreach my $lambda (@lambda){

		my $E = 1;
		my @equi_templ_file = split(/\,/, $c_in{'equi_templ'});
		my ($ntx, $irest);

		# Determine whether one or two equilibration input files are generated
		my @prefix;
		if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
			@prefix = ($c_in{'Lalias_start'}."_".$c_in{'Lalias_end'});
		}
		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			@prefix = ($c_in{'Lalias_start'}, $c_in{'Lalias_end'});
		}
		
		for my $equi_templ_file (@equi_templ_file){
		
			my $masks_found = 0;

			for my $pre (@prefix){
				my $str_out = "";
						
				open(EQUI_TEMPL, $equi_templ_file) || die "Cannot open file ".$equi_templ_file." for reading.\n";
				while(my $templ = <EQUI_TEMPL>){
					chomp($templ);
					
					# Specify coordinate reading mode
					if($first_lambda == 1){
						# Do not read the initial coordinates, velocities and box size from the inpcrd file
						# Check if ntx = 1 in template file otherwise set ntx = 1 for first equilibration step 
						if($templ =~ m/ntx /){
							my @templ_components = split(/\,/, $templ);
							foreach my $templ_component (@templ_components){
								if(($templ_component ne "")&&($templ_component =~ m/ntx /)){
									$templ_component =~ m/(\d+)/;
									if($1 != 1){
										print "Warning: ntx flag not set to 1 in equilibration template file.\n";
										print "The ntx flag will be set to 1 for the first equilibration step.\n";
										$templ =~ s/$1/1/g;
									}
								}
							}
						}
						# Check if irest = 0 in template file otherwise set irest = 0 for first equilibration step
						if($templ =~ m/irest/){
							my @templ_components = split(/\,/, $templ);
							foreach my $templ_component (@templ_components){
								if(($templ_component ne "")&&($templ_component =~ m/irest/)){
									$templ_component =~ m/(\d+)/;
									if($1 != 0){
										print "Warning: irest flag not set to 0 in equilibration template file.\n";
										print "The irest flag will be set to 0 for the first equilibration step.\n";
										$templ =~ s/$1/0/g;
									}
								}
							}
						}
					}
					else{
						# Read the initial coordinates, velocities and box size from the inpcrd file
						if($templ =~ m/ntx /){
							my @templ_components = split(/\,/, $templ);
							foreach my $templ_component (@templ_components){
								if(($templ_component ne "")&&($templ_component =~ m/ntx /)){
									$templ_component =~ m/(\d+)/;
									$templ =~ s/$1/5/g;
								}
							}
						}
						if($templ =~ m/irest/){
							my @templ_components = split(/\,/, $templ);
							foreach my $templ_component (@templ_components){
								if(($templ_component ne "")&&($templ_component =~ m/irest/)){
									$templ_component =~ m/(\d+)/;
									$templ =~ s/$1/1/g;
								}
							}
						}
					}
					
					# Specify lambda
					if($templ =~ m/clambda/){
						$templ .= "0.$lambda,";
					}
					
					# Specify mask
					if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
						if($templ =~ m/scmask/){
							$masks_found++;
							if($pre eq $c_in{'Lalias_start'}){
								$templ .= "':".$c_in{'mask0'}."',";
							}
							if($pre eq $c_in{'Lalias_end'}){
								$templ .= "':".$c_in{'mask1'}."',";
							}
						}
					}
					if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
						print "Hier2\n";
						if($templ =~ m/scmask1/){
							$masks_found++;
							#$templ .= "':".$c_in{'mask0'}."',";
							$templ .= "':".$c_in{'mask0'}."',"; 
						}
						if($templ =~ m/timask1/){
							#$templ .= "':".$c_in{'mask0'}."',";
							$templ .= "':1',";
						}
						if($templ =~ m/scmask2/){
							$masks_found++;
							#$templ .= "':".$c_in{'mask1'}."',";
							$templ .= "':".$c_in{'mask1'}."',"; 
						}
						if($templ =~ m/timask2/){
							#$templ .= "':".$c_in{'mask1'}."',";
							$templ .= "':2'";
						}
					}
					
					$str_out .= $templ."\n"; 
				}
				close EQUI_TEMPL;

				if($pre eq $c_in{'Lalias_start'}){		
					open(EQUI_IN, ">".$tag_dir."/".$pre."_equi0$E"."_v0_l$lambda.in") || die "Cannot open file $pre"."_equi0$E"."_v0_l$lambda.in\n";
				}
				elsif($pre eq $c_in{'Lalias_end'}){
					open(EQUI_IN, ">".$tag_dir."/".$pre."_equi0$E"."_v1_l$lambda.in") || die "Cannot open file $pre"."_equi$E"."_v1_l$lambda.in\n";
				}
				else{
					open(EQUI_IN, ">".$tag_dir."/md_equi0$E"."_l$lambda.in") || die "Cannot open file $pre"."_equi$E"."_l$lambda.in\n";
				}
				print EQUI_IN $str_out;
				close EQUI_IN;
			}
			# Check if two masks were found
			if($masks_found != 2){
				die("Required 'scmask' flag not found in TI equilibration input file template.\n");
			}

			
			$E++;
		}
		
		# Generation of production files for equilibration phase
		my $time = $requested_equi_time_templ*1000;
		
		for(my $p=1; $p<=$prod_no; $p++){		
			prod_setup_TI($tag_dir, $tag, $lig_res_no, $com_res_no, $time, $p, $lambda, \%c_in);
			$time = $time + ($requested_prod_time_templ*1000);
		}				
		$first_lambda = 0;
	}
	return $prod_no;
}


# Determines requested production time and returns requested time in [ns]
# Additionaly the correct setting of the shake options is checked
sub determine_requested_time{
	my $phase = shift;
	my $file_to_check = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $shake = 1;
	my $time_steps = 0;
	my $dt = 0;
	my $total_time = 0; # Total requested time in [ns]
	my $ntc = 0;
	
	if(-e $file_to_check){
		open(PROD_TEMPL, $file_to_check) || die "Cannot open file ".$file_to_check." for reading.\n";
		
		while(my $templ_line = <PROD_TEMPL>){
			chomp($templ_line);
			my @line_components = split(/\,/, $templ_line);
			
			foreach my $line_component (@line_components){
			
				if(($line_component =~ m/nstlim\s+\=/)||($line_component =~ m/nstlim\=/)){
					my @nstlim = split(/\=/, $line_component);
					$time_steps = $nstlim[$#nstlim]; 
				}
				if(($line_component =~ m/dt\s+\=/)||($line_component =~ m/dt\=/)){
					my @dt = split(/\=/, $line_component);
					$dt = $dt[$#dt];
				}
				if(($line_component =~ m/ntc\s+\=/)||($line_component =~ m/ntc\=/)){
					my @ntc = split(/\=/, $line_component);
					$ntc = $ntc[$#ntc];
				}
			}
		}
		$total_time = $time_steps * $dt / 1000;
		
		if(($c_in{'noshake'} == 1)&&($ntc != 1)){
			print "\'noshake\' option was selected, but ntc flag in $phase template is set to $ntc.\n";
			print "Please ensure that ntc and ntf are set to 1 if \'noshake\' is requested.\n";
			exit;
		}
		if(($c_in{'noshake'} == 0)&&($ntc != 2)){
			print "Shake on hydrogens bonds was requested (\'noshake\' was set to 0), but ntc flag in\n";
			print "$phase template is set to $ntc. Please ensure that ntc and ntf are set to 1 if\n";
			print "\'noshake\' is set to 0.\n";
			exit;
		}
		if(($c_in{'noshake'} == 1)&&($dt > 0.001)){
			print "The time step (dt) that is set in $phase template file is to large for a simulation\n";
			print "without shake. Please ensure consistency between the \'noshake\' setting and the parameters\n";
			print "in the $phase template file.\n";
			exit;
		}  
	}
	return($total_time);
}


# Prepare pbs-script for equilibration
sub generate_pbs_equi_TI{
	my $equi_tag_dir = shift;
	my $pbs_work_dir = shift;
	my $tag = shift;
	my $prod_no = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $out_str = "";
	
	my @lambda = split(/\,/, $c_in{'equi_lambda'});
	my $lambda_diff = $lambda[1] - $lambda[0];
	
	open(PBS_EQUI_T, $c_in{'pbs_equi_t'}) || die "Cannot open file ".$c_in{'pbs_equi_t'}." for reading.\n";
	while(my $pbs_l = <PBS_EQUI_T>){
		chomp $pbs_l;
		
		if($pbs_l =~ m/ -N/){
			$pbs_l .= " Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag;
		}
		if($pbs_l =~ m/#PBS -o/){
			$pbs_l .= " Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".OU";
		}
		if($pbs_l =~ m/#PBS -e/){
			$pbs_l .= " Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".ER";
		}
		#START - For SLUM script support
		if($pbs_l =~ m/ --job-name=/){
			$pbs_l .= "Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag;
		}
		if($pbs_l =~ m/ --error=/){
			$pbs_l .= "Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".OU";
		}
		if($pbs_l =~ m/ --output=/){
			$pbs_l .= "Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".ER";
		}
		#END - For SLUM script support
		if($pbs_l =~ m/^set SCRIPT/){
			$pbs_l .= $pbs_work_dir."/run.pbs";
		}
		if(($pbs_l =~ m/^set WORKDIR/)||($pbs_l =~ m/^set PATH/)){
			$pbs_l .= $pbs_work_dir;
		}
		if($pbs_l =~ m/^set DIFF/){
			$pbs_l .= $lambda_diff;
		}
		if($pbs_l eq "@ START="){
			$pbs_l .= $lambda[0];
		}
		if($pbs_l eq "@ END="){
			my $total_number = (@lambda * $lambda_diff)+1;
			$pbs_l .= $total_number;
		}
		if($pbs_l =~ m/END=\$START \+/){
			$pbs_l .= " ".$lambda_diff;
		}
		if($pbs_l =~ m/set INPCRD=/){
			$pbs_l .= $c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.crd";
		}
		if($pbs_l =~ m/^EQUILIBRATION/){
			if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
				$pbs_l .= "# Equilibration\n";
				$pbs_l .= "echo Equilibration_md_equi01_l".$lambda[0].".\n";
				$pbs_l = "\$DO_PARALLEL \$EXE -O -i md_equi01_l".$lambda[0].".in -p ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.top -c \${INPCRD} -o md_equi01_l".$lambda[0].".out -r md_equi01_l".$lambda[0].".rst -x md_equi01_l".$lambda[0].".mdcrd\n";
			}
		}
	
		if($pbs_l =~ m/^PRODUCTION/){
			if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
				$pbs_l = "";		
				for(my $p=1; $p<=$prod_no; $p++){
					$pbs_l .= "# Production ".$p."\n";
					$pbs_l .= "set PROD_".$p."=group_prod0".$p."_l".$lambda[0]."\n";
					$pbs_l .= "echo PROD_".$p." \$PROD"."_".$p."\n";
					$pbs_l .= "\$DO_PARALLEL \$EXE -O -ng 2 -groupfile \$PROD"."_".$p."\n";
				}
			}
			if((exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 1)){
				$pbs_l = "";
				for(my $p=1; $p<=$prod_no; $p++){
					$pbs_l .= "# Production ".$p."\n";
					$pbs_l .= "echo Prod_".$p." \n";
					my $restrt = "";
					if($p == 1){
						$restrt = "md_equi01_l".$lambda[0].".rst";
					}
					if($p == 2){
						$restrt = "md_prod01_l".$lambda[0].".rst";
					}
					if($p == 3){
						$restrt = "md_prod02_l".$lambda[0].".rst";
					}
					if($p == 4){
						$restrt = "md_prod03_l".$lambda[0].".rst";
					}
					$pbs_l .= "\$DO_PARALLEL \$EXE -O -i md_prod0".$p."_l".$lambda[0].".in -p ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_$tag"."_TIin.top -c ".$restrt." -o md_prod0".$p."_l".$lambda[0].".out -r md_prod0".$p."_l".$lambda[0].".rst -x md_prod0".$p."_l".$lambda[0].".mdcrd\n";
				}
			}
		}


		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			if($pbs_l =~ m/_lX/){
				$pbs_l =~ s/\_lX/\_l$lambda[0]/g;
			}
		
			if(($pbs_l =~ m/^set C2=/)&&($pbs_l =~ m/group_equi01_l/)){
				$pbs_l =~ s/COUNTEXP/COUNT/;
			}
		}

		if($pbs_l =~ m/^set C3=PRODUCTION/){
			$pbs_l = "";
			for(my $p=1; $p<=$prod_no; $p++){
				my $c = $p+2;
				my $lambda_length = length($lambda[0]);
				if($lambda_length == 1){
					if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
						$pbs_l .= "set C".$c."=\"s/^set PROD_".$p."=.*/set PROD_".$p."=group_prod0".$p."_l\"\${COUNT}\"/\"";
					}
					if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
						$pbs_l .= "set C".$c."=\"s/md_prod0".$p."_l\"\${LAMBDA_OLD}\"/md_prod0".$p."_l\"\${COUNT}\"/g\"";
					}
				}
				else{
					if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
						$pbs_l .= "set C".$c."=\"s/^set PROD_".$p."=.*/set PROD_".$p."=group_prod0".$p."_l\"\${COUNTEXP}\"/\"";
					}
					if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
						$pbs_l .= "set C".$c."=\"s/md_prod".$p."_l\"\${LAMBDA_OLD}\"/md_prod".$p."_l\"\${COUNT}\"/g\"";
					}
				}
				if($p != $prod_no){
					$pbs_l .= "\n";
				}
			}
			$pbs_l .= "\nset C7=\"s/^set INPCRD=.*/set INPCRD=md_prod04_l\"\${LAMBDA_OLD}\".rst/\"\n";
		}
		
		if($pbs_l =~ m/^cat \$SCRIPT/){
			for(my $p=1; $p<=($prod_no+1); $p++){
				my $c = $p+2;
				$pbs_l .= "| sed -e \"\$C".$c."\"";
			}
			$pbs_l .= " > tmp";
		}
		$out_str .= $pbs_l."\n";
	}
	close PBS_EQUI_T;

	# Write equilibration batch file	
	open(PBS_EQUI_N, ">".$equi_tag_dir."/run.pbs") || die "Cannot open file $equi_tag_dir/run.pbs for writing.\n";
	print PBS_EQUI_N $out_str;
	close PBS_EQUI_N;
	
	chmod 0755, "$equi_tag_dir/run.pbs";		
}	


# Check equilibration for completion of equilibrating phase
sub check_TIequi{
	my $transform_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $total_passed = 0;
	my %passed;
	my $check_TIequi = 0;
	my $equi_region = 0;
	my @sorted_lambdas;
	
	my @tags = ("com", "lig");

	# Set V and end alicas (state V1) for sander, not necessary 
	# for pmemd, since there only one state exits and pmemd_flag=1
        my $v = "";
	my $alias = "";
	my $pmemd_flag = 0;
	if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		$v = "\_v1";
		$alias = $c_in{'Lalias_end'};
		$pmemd_flag = 0;
	}
	if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
		$v = "";
		$alias = "md_";
		$pmemd_flag = 1;
	}

	print "Running equilibration quality check.\n";
	
	foreach my $tag (@tags){
		if($tag eq "com"){
			print "COMPLEX\n";
		}
		else{
			print "LIGAND\n";
		}
	
		# Identify MD output files of equilibration
		my @TIequi_out_files = <$transform_dir/equi/$tag/*.out>;
		
		my %TIequi_out_files;
		my @lambdas = ();
		my $last_equi_file_no = 0;
		my $last_prod_file_no = 0;
		
		foreach my $file (@TIequi_out_files){
			my @basename = split(/\//, $file);
		
			if($basename[$#basename] =~ m/$alias/){
	
				$basename[$#basename] =~ m/$v\_l(\d+)/;
				my $lambda = $1;
				my %lambdas;
			
				# Check if lambda was already stored
				map{$lambdas{$_} = 1} @lambdas;
				if(! exists $lambdas{$lambda}){
					push(@lambdas, $lambda);
				}
			
				# If part of production part of equilibration
				if($basename[$#basename] =~ m/prod/){
					$basename[$#basename] =~ m/prod(\d+)/;
					my $prod_no = $1;
					my $prod_str = "prod".$prod_no;
					$TIequi_out_files{$lambda}->{$prod_str} = $basename[$#basename];
					if($prod_no > $last_prod_file_no){
						$last_prod_file_no = $prod_no;
					}
				}
			
				# If part of equilibration part of equilibration
				if($basename[$#basename] =~ m/equi/){
					$basename[$#basename] =~ m/equi(\d+)/;
					my $equi_no = $1;
					my $equi_str = "equi".$equi_no;
					$TIequi_out_files{$lambda}->{$equi_str} = $basename[$#basename];
					if($equi_no > $last_equi_file_no){
						$last_equi_file_no = $equi_no;
					}
				}
			}
		}
		
		# Sort lambdas
		@sorted_lambdas = ();
		@sorted_lambdas = sort {$a <=> $b} @lambdas;
		
		# Read dVdL values from equilibration part of equilibration
		print "Reading dV/dL data ...\n";
		my $path_to_file = $transform_dir."/equi/".$tag."/";
		my ($dVdL_ref_equi, $time_inter_ref_equi, $dVdL_value_no_ref_equi) = read_dVdL_data($last_equi_file_no, "equi", $path_to_file, \@sorted_lambdas, \%TIequi_out_files, $pmemd_flag);
		my %dVdL_equi = %{$dVdL_ref_equi};
		my %time_inter_equi = %{$time_inter_ref_equi};
		my %dVdL_value_no_equi = %{$dVdL_value_no_ref_equi};
		my ($dVdL_ref_prod, $time_inter_ref_prod, $dVdL_value_no_ref_prod) = read_dVdL_data($last_prod_file_no, "prod", $path_to_file, \@sorted_lambdas, \%TIequi_out_files, $pmemd_flag);
		my %dVdL_prod = %{$dVdL_ref_prod};
		my %time_inter_prod = %{$time_inter_ref_prod};
		my %dVdL_value_no_prod = %{$dVdL_value_no_ref_prod};
		
		# Determine autocorrelation time based on last 1/5 of equilibration phase,
		# but at least 100 values.
		foreach my $l (@sorted_lambdas){
			print "Processing lambda: 0.".$l."\n";
			my $one_fifth = $dVdL_value_no_equi{$l} / 5;
			my @dVdL_data = ();
			$one_fifth = sprintf("%.u", $one_fifth);

			if($one_fifth >= 100){
				for(my $n=($dVdL_value_no_equi{$l}-$one_fifth+1); $n<=$dVdL_value_no_equi{$l}; $n++){
					my $time = $n * $time_inter_equi{$l};
					$time = sprintf("%.3f", $time);
					push(@dVdL_data, $dVdL_equi{$l}->{$time});
				}
			}
			else{
				
				# Considering the production part of the equilibration phase in the determination of
				# the autocorrelation time is avoided, since dt and ntpr can be different in the production
				# and usually dV/dL values are recorded in larger intervals in the production than in 
				# the equilibration (intervals of 0.5 - 1ps are common), which leads to a decrease in the 
				# accuracy of the determined autocorrelation time.
				
				if(($dVdL_value_no_equi{$l} < 100)&&($time_inter_equi{$l} >= 0.002)){
					print "WARNING: Your equilibration run should contain at least 100 dVdL values.\n";
					print "         If you do not wish to prolong the equilibration, write dVdL values\n";
					print "         more often, by decreasing 'ntpr'.\n";
					print "         Due to the small number of dVdL values the quality check of the\n";
					print "         equilibration could not be performed.\n";
					return 0;
				}
				
				else{
					my $percentage = 100 / ($dVdL_value_no_equi{$l} / 100);
					print "WARNING: The last 100 values of the equilibration used for the determination of the\n";
					print "         autocorrelation time, represent only ".$percentage." % of your total\n";
					print "         equilibration at lambda ".$l."\n";
				
					for(my $n=$dVdL_value_no_equi{$l}; $n>($dVdL_value_no_equi{$l}-100); $n--){
						my $time = $n * $time_inter_equi{$l};
						push(@dVdL_data, $dVdL_equi{$l}->{$time});
					}
				}
			}
			
			my @res;
			my $autocorr = 0;
			my @data1 = @dVdL_data;
			my @data2 = @dVdL_data;
			$autocorr = determine_autocorrelation(\@data1, \@data2, \@res);

			my $autocorr_time = $autocorr * $time_inter_equi{$l};
			if($autocorr_time <= 1){
				$autocorr = 1/$time_inter_equi{$l};
				$autocorr = sprintf("%.0f", $autocorr);
			}
			
			# Since at least 3 windows a 20 x autocorr shall be checked,
			# it is required that the total number of dVdL values is equal to or
			# larger than 3 x (20 x autocorr). -> Commonly >= 60 ps
			if($dVdL_value_no_equi{$l} < (3 * (20*$autocorr))){
				print "WARNING: The total number of dVdL values recorded in the equilibration of lambda ".$l."\n";
				print "         is smaller than 3 times (20 * autocorrelation time of ".$autocorr_time." ps)\n";
				print "         Therefore the quality check of the equilibration could not be performed.\n";
				return 0;
			}
			
			my ($passed, $equi_end) = check_TIequi_termination($l, $autocorr, $dVdL_value_no_equi{$l}, $tag, $time_inter_equi{$l}, 0, \%dVdL_equi);
			
			# Determine maximum time of equilibration phase
			# found for ligand and complex.
			my $end_equi_time = $equi_end * $time_inter_equi{$l};
			if($end_equi_time > $equi_region){
				$equi_region = $end_equi_time;
			}
			
			# If equilibration check was not passed consider also production phase of equilibration in check.
			# Problematic: Time intervals in which dV/dL data are recorded in production phase can differ 
			# from equilibration phase.
			# Therefore, the tau interval is set to a value larger the autocorrelation time which is a multiple
			# of the time interval in which the dV/dL data are recorded in the production.			
			if(($passed == 0)&&(exists $TIequi_out_files{$l}->{'prod01'})){
				my %dVdL_sum;
								
				# Determine new interval time (tau)
				my $time_inter_prod = $time_inter_prod{$l};

				while($time_inter_prod < ($autocorr*$time_inter_equi{$l})){
					$time_inter_prod = $time_inter_prod + $time_inter_prod{$l};
				}
				$time_inter_prod = sprintf("%.3f", $time_inter_prod);
				my $time_inter_prod_fs = $time_inter_prod * 1000; # Convert to femtoseconds
				$time_inter_prod_fs = sprintf("%.0f", $time_inter_prod_fs);
				
				# Generate new hash with combined values from both equilibration
				# and production phase of equilibration.
				# Do not consider dVdL values recorded during equilibration phase
				# that are saved more often than the new interval time.
				my $largest_t_equi = 0;
				my $count_data_equi = 0;
				foreach my $t_equi (keys %{$dVdL_equi{$l}}){
					my $fs_time = $t_equi * 1000;
					$fs_time = sprintf("%.0f", $fs_time);
					if($fs_time % $time_inter_prod_fs){
						# Do not store those that cannot be devided
						# by time interval of production that is larger
						# than autocorrelation time
					}
					else{
						$dVdL_sum{$l}->{$t_equi} = $dVdL_equi{$l}->{$t_equi};
						if($t_equi > $largest_t_equi){
							$largest_t_equi = $t_equi;
						}
						$count_data_equi++;
					}
				}

				foreach my $t_prod (keys %{$dVdL_prod{$l}}){
					my $fs_time = $t_prod * 1000;
					$fs_time = sprintf("%.0f", $fs_time);
					if($fs_time % $time_inter_prod_fs){
						# Do not store those that cannot be devided
						# by time interval of production that is larger
						# than autocorrelation time
					}
					else{
						my $t_prod_new = $t_prod + $largest_t_equi;
						$t_prod_new = sprintf("%.3f", $t_prod_new);
						$dVdL_sum{$l}->{$t_prod_new} = $dVdL_prod{$l}->{$t_prod};
					}
				}

				my @dVdL_sum = keys %{$dVdL_sum{$l}};
				my $total_values_sum = @dVdL_sum;
				
				my $start_increase = $count_data_equi;
					
				($passed, $equi_end) = check_TIequi_termination($l, 1, $total_values_sum, $tag, $time_inter_prod, $start_increase, \%dVdL_sum);					

				# Determine maximum time of equilibration phase
				# found for ligand and complex.
				my $end_equi_time = $equi_end * $time_inter_prod;
				if($end_equi_time > $equi_region){
					$equi_region = $end_equi_time;
				}
			}

			if($passed == 1){
				print "Passed!\n";
				$passed{$tag}->{$l} = $passed;
				$total_passed++;
			}
		} # End: Do for each lambda
	} # End: Do for complex and receptor
	
	# Equilibration check failed, if the confidence criterion is not fulfilled for less than 50%
	# of all simulations. If not all simulations fulfill the confidence criterion, provide a note to the user.
	if(($total_passed >= (@sorted_lambdas * 2)/2)&&($total_passed < (@sorted_lambdas * 2))){
		print "Equilibration check was successful, but please have a closer look at equilibrations:\n";
		$check_TIequi = 1;
		foreach my $tag (@tags){
			foreach my $l (@sorted_lambdas){
				if($passed{$tag}->{$l} == 0){
					if($tag eq "lig"){
						print "  - Ligand ";
					}
					else{
						print "  - Complex ";
					}
					print "at lambda=0.$l\n";
				}
			}
		}
		print "\n";
		print "Estimate of equilibrating region: First $equi_region ps of equilibration simulation.\n\n";
	}
	if($total_passed == (@sorted_lambdas * 2)){
		$check_TIequi = 1;
		print "Estimate of equilibrating region: First $equi_region ps of equilibration simulation.\n\n";
	}
	return $check_TIequi;
}


# Reading of dV/dL data from file and return
# 1.) Hash with dVdL_data
# 2.) Hash of time intervals per lambda
# 3.) Hash of total dVdL values per lambda
sub read_dVdL_data{
	my $last_file_no = shift;
	my $file_id = shift;
	my $path_to_file = shift;
	my $sorted_lambdas_ref = shift;
	my @sorted_lambdas = @{$sorted_lambdas_ref};
	my $TI_out_files_ref = shift;
	my %TI_out_files = %{$TI_out_files_ref};
	my $pmemd_flag = shift;

	my %dVdL;
	my %time_inter;
	my %dVdL_value_no;

	foreach my $l (@sorted_lambdas){
		my $save_time = 0;
		my $time_offset = 0;
		my $dVdL_value_no = 0;
		my $time_interval_detected = 0;
		my $time_inter = 0;
		my $time = 0;
		my $delta_time = 0;
		
		for(my $i=1; $i<=$last_file_no; $i++){
			my $id = $file_id;
			if($i < 10){
				$id .= "0".$i;
			}
			else{
				$id .= $i;
			}
			my $section_two = 0;
			my $begin_read = 0;
			my $begin_read_two = 0;
			open(OUT_FILE, $path_to_file.$TI_out_files{$l}->{$id}) || die "Cannot open file ".$TI_out_files{$l}->{$id}."\n";
			while(my $f_line = <OUT_FILE>){
				chomp($f_line);
				
				if($f_line =~ m/2\.  CONTROL  DATA  FOR  THE  RUN/){
					$section_two = 1;	
				}

				if(($section_two == 1)&&($f_line =~ m/dt\s+=\s+(\d+)\.(\d+)/)){
					$delta_time = $1.".".$2;
					$section_two = 0;
				}

				if($f_line =~ m/4\.  RESULTS/){
					$begin_read = 1;
				}

				if($f_line =~ m/TI region  1/){
					$begin_read_two = 1;
				}

				if($f_line =~ m/A V E R A G E S/){
					$begin_read = 0;
					close OUT_FILE;
					last;
				}

				if($f_line =~ m/TI region  2/){
					$begin_read_two = 0;
				}
		
				if((($pmemd_flag == 0)&&($begin_read == 1)&&($f_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))
				||(($pmemd_flag == 1)&&($begin_read_two == 1)&&($f_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))){
					$time = $2.".".$3;
		
					if($save_time == 0){
						$f_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/;
						my $nstep = $1;
						$time_inter = $delta_time * $nstep;
						$time_offset = $time - $time_inter;
					}
					$time = $time - $time_offset;
					$time = sprintf("%.3f", $time);
						
					$save_time = $time;	
				}

				if((($pmemd_flag == 0)&&($begin_read == 1)&&($f_line =~ m/^ DV\/DL/))
				||(($pmemd_flag == 1)&&($begin_read_two == 1)&&($f_line =~ m/^ DV\/DL/))){
					my @value = split(/\s+/, $f_line);
					if($time != 0){
						$dVdL{$l}->{$time} = $value[3];
						$dVdL_value_no++;
					}
				}
			} # End: Reding of output file
			close OUT_FILE;
		} # End: Processing output files of one lambda
		$dVdL_value_no{$l} = $dVdL_value_no;
		$time_inter{$l} = $time_inter;
	} # End: Processing all lambdas
	return(\%dVdL, \%time_inter, \%dVdL_value_no);
}


# Check termination of equilibration
# This is fulfilled, if the equilibrating period is complete
# Check is based on the procedure described in:
# W. Yang, R. Bitetti-Putzer, M. Karplus, Free energy simulations: Use of
# reverse cummulative averaging to determine the equilibrated region and 
# the time required for convergence. J. Chem. Phys. 120(6), 2618-2628.
# The confidence is determined for increasing sections of the equilibration
# in reverse order, starting from the last recoreded dVdL value.
# Sections of dVdL values regarded are increased in intervals of 
# 5 x (20 x autocorrelation time).
sub check_TIequi_termination{
	my $lambda = shift;
	my $initial_autocorr = shift;
	my $dVdL_value_no = shift;
	my $tag = shift;
	my $time_inter = shift;
	my $start_increase = shift;
	my $dVdL_ref = shift;
	my %dVdL = %{$dVdL_ref};
	my $initial_tau = ($initial_autocorr * $time_inter) * 1000;
	my $count_x = 0;
	my $count_y = 0;
	my $y = 0;
	my %y;
	my $x_start = $dVdL_value_no % $initial_autocorr;
	my @ystep;
	
	# Prepare coarse grained dataset Y
	# Y data are read from the start of the simulation,
	# but the first data points that cannot be devided by
	# the autocorrelation time are skipped.
	for(my $i=($x_start+1); $i<=$dVdL_value_no; $i++){
		my $time = $i * $time_inter;
		$time = sprintf("%.3f", $time);
		my $fs_time = $time * 1000;
		$fs_time = sprintf("%.0f", $fs_time);
		if($fs_time % $initial_tau){
			 $y = $y + $dVdL{$lambda}->{$time};
			 $count_x++;
		}
		else{
			$y = $y + $dVdL{$lambda}->{$time};
			$count_x++;
			$count_y++;
			$y = $y / $count_x;			
			push(@ystep, $count_y);
			$y{$count_y} = $y;
			$y = 0;
			$count_x = 0;
		}
	}
	
	# Analysis
	my $total_y_steps = @ystep;
	my $inter_steps = 20;
	my $inter_count = 0;
	my $stop_iterative_increase = 0;
	my $curr_nsteps_no = 0;
	my $check_passed = 0;
	my $equi_end;
	my $passed_without_equi = 0; # To ensure that search is carried on as long as no
								 # equilibration interval is detected, information that
								 # check was passed once is stored.

	while($stop_iterative_increase == 0){

		$inter_count++;

		# Prepare output
		$curr_nsteps_no = $start_increase + ($inter_count * $inter_steps);
		print "Processing interval ".$inter_count." with steps: ".$curr_nsteps_no."\n";

		# Perform Shapiro-Wilk-Test for k-subsets		
		($check_passed, $equi_end) = run_shapiro_wilk_test_for_k_subsets($curr_nsteps_no, $inter_steps, $initial_autocorr, \@ystep, \%y);

		if(($check_passed == 1)&&($equi_end == 0)){
			$passed_without_equi = 1;
		}

		if(($check_passed == 1)&&($equi_end > 0)){
			last;
		}
		
		# Terminate iteration, if interval is larger than remaining y-steps
		if(($total_y_steps - $curr_nsteps_no) < $inter_steps){
			$stop_iterative_increase = 1;
		}
	}
	
	# Ensure that check is passed, if it was passed once, even if equilibration
	# period was not detected.
	if($passed_without_equi){
		$check_passed = 1;
	}
	
	# If dVdL data were skipped, because the total data could not be devided
	# evenly into subsets of X times tau, add the time frame of the skipped
	# data to the time of the equilibrating region
	if($x_start != 0){
		$equi_end = $equi_end + ($x_start * $time_inter); 
	}
	
	return ($check_passed, $equi_end);
} 


# Determination of autocorrelation time
# Original time correlation function written by H. Gohlke
sub determine_autocorrelation{
	my $autocorrelation_time;
	my $r_data1 = shift;
	my $r_data2 = shift;
	my $r_res   = shift;
	@$r_res = ();

	my $len1 = scalar(@$r_data1);
	my $len2 = scalar(@$r_data2);
	my $len;
	if ($len1 > $len2){
		$len = $len1;
	}
	else{
		$len = $len2;
	}

	# data1 lags data2 (i.e. is shifted to the right of it)
	# -> peak in positive lags (i.e. 0 .. $len-1 in @res)
	# data2 lags data1 (i.e. is shifted to the right of it)
	# -> peak in negative lags (i.e. $len .. 2*$len-1 in @res)
	#
	# * data1/2 are centered
	# * correlation values are finally normalized
	my $i;
	for($i = 0;$i < 2*$len; $i++){
		push @$r_res,0.0;
	}
	my $mean1 = 0.0;
	for($i=0; $i < $len1; $i++){
		$mean1 += $r_data1->[$i];
	}
	$mean1 /= $len1;

	my $mean2 = 0.0;
	for($i=0; $i < $len2; $i++){
		$mean2 += $r_data2->[$i];
	}
	$mean2 /= $len2;

	my $lag;
	my $norm1 = 0;
	my $norm2 = 0;
	my $norm;
	for($lag = 0; $lag < $len; $lag++){
		for($i = 0; $i < $len; $i++){
			$r_res->[$lag] +=	($i+$lag >= $len1 ? 0.0 : $r_data1->[$i+$lag]-$mean1) *
								($i      >= $len2 ? 0.0 : $r_data2->[$i]-$mean2);
			$r_res->[$len+$lag] +=	($i      >= $len1 ? 0.0 : $r_data1->[$i]-$mean1) *
									($i+$lag >= $len2 ? 0.0 : $r_data2->[$i+$lag]-$mean2);
			if($lag == 0){
				$norm1 += ($r_data1->[$i]-$mean1) * ($r_data1->[$i]-$mean1);
				$norm2 += ($r_data2->[$i]-$mean2) * ($r_data2->[$i]-$mean2);
			}
		}
		if($lag == 0){
			$norm = sqrt($norm1 * $norm2);
		}
		$r_res->[$lag] /= $norm;
		$r_res->[$len+$lag] /= $norm;
	}
	
	my $e_factor = 1/exp(1);
	for($lag = 0; $lag < $len; $lag++){
		if($r_res->[$lag] < $e_factor){
			$autocorrelation_time = $lag;
			last;
		}
	}
	
	return $autocorrelation_time;
}


# Calculate W and confidence for increasing k
# starting from k=6, since this is the lower limit for data input that can
# be handled by the perl module Statistics::Normality::shapiro_wilk_test
sub run_shapiro_wilk_test_for_k_subsets{
	
	my $curr_nsteps_no = shift;
	my $inter_steps = shift;
	my $initial_autocorr = shift;
	my $ystep_ref = shift;
	my $y_ref = shift;
	my @ystep = @{$ystep_ref};
	my %y = %{$y_ref};
	my $check_passed = 0;
	my $count_conf = 0;	# Variable for counting the number of successive k
						# with confidence above the threshold of 80%
	my $equi_end = 0; # Variable for storing the first drop of confidence
					  # after equilibration termination criterion is passed.
	my $tmp_lower = 0;
	my $last_conv = 0;
		
	for(my $k=6; $k<=$curr_nsteps_no; $k++){
	
		# Generate dataset for the respective k
		my %data_by_k;
		my $avg_y = 0;
		for(my $i=$curr_nsteps_no; $i>($curr_nsteps_no-$k); $i--){
			$data_by_k{$i} = $y{$ystep[$i-1]};
			$avg_y = $avg_y + $y{$ystep[$i-1]};
		}
		$avg_y = $avg_y / $k;
	
		# Sort data
		my @sorted = sort { $data_by_k{$a} <=> $data_by_k{$b} } keys %data_by_k;
		my @sorted_data_by_k = ();
		foreach my $s (@sorted){
			push(@sorted_data_by_k, $data_by_k{$s});
		}
	
		# Perform Shapiro-Wilk test for normality and determine confidence
		my ($conf, $W) = Statistics::Normality::shapiro_wilk_test(\@sorted_data_by_k);

		# Check if confidence threshold is passed
		if($conf >= 0.8){
			$count_conf++;
		}
		if($conf < 0.8){
			if($check_passed == 1){
				$equi_end = (($curr_nsteps_no-$k-1)*$initial_autocorr);
				if($equi_end < 0){
					$equi_end = 0;
				}
				return ($check_passed, $equi_end);

			}
			if(($count_conf > 0)&&($tmp_lower >= 2)){
				$count_conf = 0;
				$tmp_lower = 0;
			}
			if(($count_conf > 0)&&($tmp_lower < 2)){
				$tmp_lower++;
			}
		}
			
		# Check if equilibration criterion is fulfilled
		# Fulfilled if a) Number of subsets with confidence >=80 % at leat 20
		#          and b) to ensure that high confidence interval is not
		#                 accidentally found directly at the beginning
		#                 the number of k-subsets left should be larger than 10.
		if(($count_conf == $inter_steps)&&($curr_nsteps_no-$k > 10)){
			$check_passed = 1;
		}
	
		# Prepare for next round
		foreach my $i (keys %data_by_k){
			delete $data_by_k{$i};
		}
	}
	return ($check_passed, $equi_end);
}


# Identify mask for cpptraj run - containing all atoms that are
# part of both V0 and V1
sub identify_mask{
	my $setup_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $mask = "";
	
	my @org_mask = split(/\@/, $c_in{'mask1'});
	my @diff_atoms = split(/\,/, @org_mask[1]);
	
	open(V_ONE_PDB, $setup_dir."/".$c_in{'Lname_end'}."_lig.pdb") || die "Cannot open file ".$setup_dir."/".$c_in{'Lname_end'}."_lig.pdb for reading 9\n";
	while(my $pdb_l = <V_ONE_PDB>){
		my @pdb_l = split(/\s+/, $pdb_l);
		
		if($pdb_l[3] eq $c_in{'Lalias_end'}){		
			if(grep {$pdb_l[2] eq $_} @diff_atoms)
			{}
			else{
				if($mask ne ""){
					$mask .= ",$pdb_l[2]";
				}
				if($mask eq ""){
					$mask .= $pdb_l[2];
				}
			}
		}
	}
	close V_ONE_PDB;
	return $mask;
} 


# Generate cpptraj template file for imaging
sub generate_cpptraj_img{
	my $MD_dir = shift;
	my $tag = shift;
	my $com_res_lig = shift;
	my $mask = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	open(CPPTRAJ_T, ">".$MD_dir."/cpptraj_img_template.ptrj") || die "Cannot open file ".$MD_dir."/cpptraj_img_template.ptrj for writing.\n";
	if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		print CPPTRAJ_T "trajin ./PREFIX_prodN_v_lX.rst\n";
	}
	if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
		print CPPTRAJ_T "trajin ./PREFIX_prodN_lX.rst\n";
	}
	if($tag eq "lig"){
		print CPPTRAJ_T "center :1@".$mask." origin\n";
	}
	if($tag eq "com"){
		my $rec_res_no = $com_res_lig - 1;
		print CPPTRAJ_T "center :1-".$rec_res_no.",:".$com_res_lig."@".$mask." origin\n";
	}
	print CPPTRAJ_T "image origin center familiar\n"; # Currently no separate imaging for different chains implemented
	if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		print CPPTRAJ_T "trajout PREFIX_prodN_v_lX_img.rst restrt\n";
	}
	if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
		print CPPTRAJ_T "trajout PREFIX_prodN_lX_img.rst restrt\n";
	}
}


# Determine basename of last coordinate file from TI equilibration
# that shall be used as input for TI production.
sub detect_last_equi_basename{
	my $equi_path_and_name = shift;
	my $prod_lambda_ref = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my @prod_lambda = @{$prod_lambda_ref};
	my @equi_rst_files = <$equi_path_and_name*.rst>;

	# Set V
	my $v = "";
	my $vr = "";
	if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		$v = "\_v1";
		$vr = "_v1";
        }
        if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
                $v = "";
		$vr = "";
        }

	# Read file information
	my @rst_no;
	my %rst_no;
	my %lambdas;
	foreach my $rst_file (@equi_rst_files){
		my @filename_components = split(/\//, $rst_file);

		if($filename_components[$#filename_components] =~ m/prod(\d+)$v\_l(\d+)\.rst/){
			my $curr_no = $1;
			my $lambda = $2;
			if(! $rst_no{$curr_no}){
				push(@rst_no, $curr_no);
				$rst_no{$curr_no} = 1;
			}
			if(! $lambdas{$lambda}){
				$lambdas{$lambda} = 1;
			}
		}
	}
	
	if(@rst_no == 0){
		print "\nERROR: No coordinate output files from TI equilibration phase were found.\n";
		print "Before setup of the TI production simulations the TI equilibration simulations\n";
		print "have to be completed. Please ensure that the output files from the TI\n";
		print "equilibration are present under the default location and can be accessed.\n";
	}
	
	my @sorted_no = sort {$a <=> $b} @rst_no;
	my $last_rst_file_base = "prod".$sorted_no[$#sorted_no];
	
	# Check if last file exists for all lambda values requested
	# to be considered in production setup.
	foreach my $l (keys %lambdas){
		if(! -e $equi_path_and_name."_".$last_rst_file_base.$vr."_l".$l.".rst"){
			print "\nERROR: Equilibration run probably not completely finished.\n";
			print "The final coordinate file of the TI equilibration that is required\n";
			print "for setup of TI production was not found for lamda 0.".$l.".\n";
			print "Please check the existence and accessibility of the restart-file\n";
			print $equi_path_and_name."_".$last_rst_file_base.$vr."_l".$l.".rst\n";
			exit;
		}
	}
	
	return $last_rst_file_base;
}

# Setup of group file for production
sub group_file_prod_setup{
	my $MD_dir = shift;
	my $restrt = shift;
	my $tag = shift;
	my $p = shift;
	my $lambda = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	
	# Ensure correct naming for large number of production runs 
	my $p_alias = "";
	if($p < 10){
		$p_alias = "0$p";
	}
	if($p >= 10){
		$p_alias = $p;
	}
	
	# Previous production for restrt specification
	my $prev_prod = $p - 1;
	if($prev_prod < 10){
		$prev_prod = "0$prev_prod";
	}
	if($prev_prod >= 10){
		$prev_prod = "$prev_prod";
	}
	
	# Define restrt-files
	my $restrt_v0;
	my $restrt_v1;
	
	if($restrt ne ""){
		$restrt_v0 = $c_in{'Lalias_start'}."_".$restrt."_v0_l".$lambda.".rst";
		$restrt_v1 = $c_in{'Lalias_end'}."_".$restrt."_v1_l".$lambda.".rst";	
	}
	else{
		if($p == 1){
			$restrt_v0 = "equi_v0_l$lambda.rst";
			$restrt_v1 = "equi_v1_l$lambda.rst";
		}
		if($p > 1){
			$restrt_v0 = $c_in{'Lalias_start'}."_prod$prev_prod"."_v0_l".$lambda.".rst";
			$restrt_v1 = $c_in{'Lalias_end'}."_prod$prev_prod"."_v1_l".$lambda.".rst";
		}
	}
			
	open(GROUP_PROD, ">$MD_dir/group_prod$p_alias"."_l$lambda") || die "Cannot open file $MD_dir/group_prod$p_alias"."_l$lambda for writing\n";
	print GROUP_PROD "-O -i $c_in{'Lalias_start'}"."_prod$p_alias"."_v0_l$lambda".".in ".
	 			        "-o $c_in{'Lalias_start'}"."_prod$p_alias"."_v0_l$lambda".".out ".
					    "-p $c_in{'Lalias_start'}"."_$tag"."_TIin.top ".
						"-c $restrt_v0 ".
						"-r $c_in{'Lalias_start'}"."_prod$p_alias"."_v0_l$lambda".".rst ".
						"-x $c_in{'Lalias_start'}"."_prod$p_alias"."_v0_l$lambda".".mdcrd\n";
								
	print GROUP_PROD "-O -i $c_in{'Lalias_end'}"."_prod$p_alias"."_v1_l$lambda".".in ".
	 			        "-o $c_in{'Lalias_end'}"."_prod$p_alias"."_v1_l$lambda".".out ".
					    "-p $c_in{'Lalias_end'}"."_$tag"."_TIin.top ".
						"-c $restrt_v1 ".
						"-r $c_in{'Lalias_end'}"."_prod$p_alias"."_v1_l$lambda".".rst ".
						"-x $c_in{'Lalias_end'}"."_prod$p_alias"."_v1_l$lambda".".mdcrd\n";
	close GROUP_PROD;
}


# Setup MD-production files according to templates
sub prod_setup_TI{
	my $MD_dir = shift;
	my $tag = shift;
	my $lig_res_no = shift;
	my $com_res_no = shift;
	my $time = shift;
	my $p = shift;
	my $lambda = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	
	# Ensure correct naming for large number of production runs 
	my $p_alias = "";
	if($p < 10){
		$p_alias = "0$p";
	}
	if($p >= 10){
		$p_alias = "$p";
	}

	my @prefix;
	if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
		@prefix = ($c_in{'Lalias_start'}, $c_in{'Lalias_end'});
	}
	if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
		@prefix = ($c_in{'Lalias_start'}."_".$c_in{'Lalias_end'});
	}
			
	my $masks_found = 0;
	for my $pre (@prefix){
		my $str_out = "";
						
		open(PROD_TEMPL, $c_in{'prod_templ'}) || die "Cannot open file ".$c_in{'prod_templ'}." for reading.\n";

		while(my $templ = <PROD_TEMPL>){
			chomp($templ);
					
			# Specify lambda
			if($templ =~ m/clambda/){
				$templ .= "0.$lambda,";
			}
					
			# Specify mask
			if((!exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
				if($templ =~ m/scmask/){
					$masks_found++;
					if($pre eq $c_in{'Lalias_start'}){
						$templ .= "':".$c_in{'mask0'}."',";
					}
					if($pre eq $c_in{'Lalias_end'}){
						$templ .= "':".$c_in{'mask1'}."',";
					}
				}
					
				if(($templ =~ m/^RES/)||($templ =~ m/^LRES/)){
					if($tag eq "lig"){
						$templ .= " $lig_res_no";
					}
					if($tag eq "com"){
						$templ .= " $com_res_no";
					}
				}
			}
			if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
				if($templ =~ m/scmask1/){
					$masks_found++;
					$templ .= "':".$c_in{'mask0'}."',";
				}
				if($templ =~ m/timask1/){
					$templ .=  "':1',";
				}
				if($templ =~ m/scmask2/){
					$masks_found++;
					$templ .= "':".$c_in{'mask1'}."',";
				}
				if($templ =~ m/timask2/){
					$templ .=  "':2',";
				}
				if(($templ =~ m/^RES/)||($templ =~ m/idecomp/)){
					die("Per residue decomposition cannot be performed with pmemd.\nPlease use sander for the calculation.\n");
				}
			}
			
			if(($time != 0)&&($templ =~ m/ t /)){
				$templ =~ s/0\.0/$time/g;
			}
					
			$str_out .= $templ."\n"; 
		}
		close PROD_TEMPL;
				
		if($pre eq $c_in{'Lalias_start'}){	
			open(PROD_IN, ">".$MD_dir."/".$pre."_prod$p_alias"."_v0_l$lambda.in") || die "Cannot open file $pre"."_prod$p_alias"."_v0_l$lambda.in\n";
		}
		elsif($pre eq $c_in{'Lalias_end'}){
			open(PROD_IN, ">".$MD_dir."/".$pre."_prod$p_alias"."_v1_l$lambda.in") || die "Cannot open file $pre"."_prod$p_alias"."_v1_l$lambda.in\n";
		}
		else{
			open(PROD_IN, ">".$MD_dir."/md_prod$p_alias"."_l$lambda.in") || die "Cannot open file $pre"."_prod$p_alias"."_l$lambda.in\n";
		}
		print PROD_IN $str_out;
		close PROD_IN;
	}
	if($masks_found != 2){
		die("Required 'scmask' flag not found in TI equilibration/production input file template.\n");
	}
}


# Prepare pbs-script for running MD-production
sub generate_pbs_prod_TI{
	my $MD_dir = shift;
	my $pbs_work_dir = shift;
	my $current_dir = shift;
	my $tag = shift;
	my $lambda = shift;
	my $prod_no = shift;
	my $requested_prod_time_templ = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $out_str = "";
	my $end = $prod_no;

	my $delta_time = $requested_prod_time_templ * 1000;

	open(PBS_TEMPL, $c_in{'pbs_prod_t'}) || die "Cannot open file $c_in{'pbs_prod_t'}\n";
	while(my $pbs_l = <PBS_TEMPL>){
		chomp($pbs_l);
		if($pbs_l =~ m/ -N/){
			$pbs_l .= " ".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag."_l".$lambda;
		}
		if($pbs_l =~ m/#PBS -o/){
			$pbs_l .= " ".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag."_l".$lambda.".OU";
		}
		if($pbs_l =~ m/#PBS -e/){
			$pbs_l .= " ".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag."_l".$lambda.".ER";
		}
                if($pbs_l =~ m/#PBS -o/){
                        $pbs_l .= " Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".OU";
                }
                if($pbs_l =~ m/#PBS -e/){
                        $pbs_l .= " Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".ER";
                }
		#START - For SLUM script support
		if($pbs_l =~ m/ --job-name=/){
			$pbs_l .= $c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag."_l".$lambda;
		}
		if($pbs_l =~ m/ --error=/){
                        $pbs_l .= "Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".ER";
		}
		if($pbs_l =~ m/ --output=/){
                        $pbs_l .= "Equi_".$c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}."_".$tag.".OU";
		}
		#END - For SLUM script support
		if($pbs_l =~ m/@ END/){
			$pbs_l .= $end+1;
		}
		if(($pbs_l =~ m/^set WORKDIR/)||($pbs_l =~ m/^set PATH/)){
			$pbs_l .= $pbs_work_dir;
		}
		if(($pbs_l =~ m/^set SCRIPT/)||($pbs_l =~ m/^set GROUP_FILE/)||($pbs_l =~ m/^gzip/)){
			$pbs_l =~ s/X/$lambda/g;
		}
		if($pbs_l =~ m/tmp_rst_v/){
			$pbs_l =~ s/lX/l$lambda/;
		}
		if($pbs_l =~ m/^set CPPTRAJ_TEMPL/){
			$pbs_l .= "./cpptraj_img_template.ptrj";
		}
		if($pbs_l =~ m/^set LAMBDA/){
			$pbs_l .= $lambda;
		}
		if($pbs_l =~ m/^set TAG/){
			$pbs_l .= $tag;
		}
		if($pbs_l =~ m/^set V0/){
			$pbs_l .= $c_in{'Lalias_start'};
		}
		if($pbs_l =~ m/^set V1/){
			$pbs_l .= $c_in{'Lalias_end'};
		}
		if($pbs_l =~ m/^set VV/){
			$pbs_l .= $c_in{'Lname_start'}."_".$c_in{'Lname_end'};
		}
		if($pbs_l =~ m/^set RESTART/){
			if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
				$pbs_l .= $c_in{'Lalias_start'}."_".$c_in{'Lalias_end'}.".rst";
			}
			if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
				$pbs_l .= "equi_l".$lambda.".rst";
			}
		}
		if($pbs_l =~ m/^PRODUCTION/){
			$pbs_l = "\$DO_PARALLEL \$EXE -O -i md_prod01_l".$lambda.".in -p ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.top -c \$RESTART -o md_prod01_l".$lambda.".out -r md_prod01_l".$lambda.".rst -x md_prod01_l".$lambda.".mdcrd\n";
		}
		if($pbs_l =~ m/^set CALC_STDERR/){
			if((! exists $c_in{'conv_exe'})||($c_in{'conv_exe'} eq "")){
				$pbs_l .= "$current_dir/miscellaneous/convergenceCheck.pl";
			}
			else{
				$pbs_l .= $c_in{'conv_exe'};
			}
		}
		if($pbs_l =~ m/^set ERROR_LIMIT/){
			if(($c_in{'conv_meth'} == 1)&&(! $c_in{'error_limit'})||($c_in{'error_limit'} eq "")){
				$pbs_l .= "0.01";
			}
			elsif(($c_in{'conv_meth'} == 2)&&(! $c_in{'error_limit'})||($c_in{'error_limit'} eq "")){
				$pbs_l .= "0.2";
			}
			else{
				$pbs_l .= $c_in{'error_limit'};
			}
		}
		if($pbs_l =~ m/^set CONV_METH/){
			$pbs_l .= $c_in{'conv_meth'};
		}
		if($pbs_l =~ m/@ DELTATIME/){
			$pbs_l .= $delta_time;
		}
		
		$out_str .= $pbs_l."\n";
	}
	close PBS_TEMPL;
		
	open(PBS, ">".$MD_dir."/run_prod.pbs.l".$lambda) || die "Cannot open file run_prod.pbs.l$lambda for writing.\n";
	print PBS $out_str."\n";
	close PBS;
	
	chmod 0755, "$MD_dir./run_prod.pbs.l.$lambda";
}

1;
