use strict;

# Write first leap-input script for setup of topology and coordinate files
# of system in vacuum
sub create_leap_in_top_crd{
	my $cryst_dir = shift;
	my $curr_leap_dir = shift;
	my $struct_b = shift;
	my $chrg_meth = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my %disulf;
	
	my @unit_names = ("COM", "REC");

	open(LEAP_IN, ">$cryst_dir/leap_script_1.in") || die "Cannot open $cryst_dir/leap_script_1.in";

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
	        #print LEAP_IN "source leaprc.water.tip3p\n"; 
	        print LEAP_IN "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";
	}
	print LEAP_IN "source leaprc.gaff\n";
	print LEAP_IN "loadoff ".$curr_leap_dir."/".$struct_b."_".$chrg_meth.".lib\n";
			
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
	print LEAP_IN "LIG = loadpdb ".$curr_leap_dir."/".$struct_b."_lig.pdb\n";
	print LEAP_IN "REC = loadpdb ".$curr_leap_dir."/".$struct_b."_rec.pdb\n";
	print LEAP_IN "COM = loadpdb ".$curr_leap_dir."/".$struct_b."_com.pdb\n";
	
	if(keys %disulf){
		foreach my $unit_name (@unit_names){
			foreach my $k (keys %disulf){
				print LEAP_IN "bond ".$unit_name.".".$k.".SG ".$unit_name.".".$disulf{$k}.".SG\n";
			}
		}
	}
	print LEAP_IN "loadAmberParams ".$curr_leap_dir."/".$struct_b.".frcmod\n";
	print LEAP_IN "saveAmberParm LIG $cryst_dir/".$struct_b."_vac_lig.top $cryst_dir/".$struct_b."_vac_lig.crd\n";
	print LEAP_IN "savepdb LIG $cryst_dir/".$struct_b."_vac_lig.pdb\n";
	print LEAP_IN "saveAmberParm REC $cryst_dir/".$struct_b."_vac_rec.top $cryst_dir/".$struct_b."_vac_rec.crd\n";
	print LEAP_IN "savepdb REC $cryst_dir/".$struct_b."_vac_rec.pdb\n";
	print LEAP_IN "saveAmberParm COM $cryst_dir/".$struct_b."_vac_com.top $cryst_dir/".$struct_b."_vac_com.crd\n";
	print LEAP_IN "savepdb COM $cryst_dir/".$struct_b."_vac_com.pdb\n";
	print LEAP_IN "charge LIG\n";
	print LEAP_IN "charge REC\n";
	print LEAP_IN "charge COM\n";
	print LEAP_IN "quit\n";
	close LEAP_IN;
}


# Write first leap-input script for setup of topology and coordinate files
# of system to determine total charge of membrane system
sub create_mem_leap_in_top_crd{
	my $cryst_dir = shift;
	my $curr_leap_dir = shift;
	my $struct_b = shift;
	my $chrg_meth = shift;
	my $curr_FEW_dir = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my %disulf;
	
	my @unit_names = ("COM", "REC");

	open(LEAP_IN, ">$cryst_dir/leap_script_1.in") || die "Cannot open $cryst_dir/leap_script_1.in";

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
		print LEAP_IN "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff10\n";
		print LEAP_IN "loadAmberParams frcmod.ionsjc_tip3p\n";
	}
	print LEAP_IN "source leaprc.gaff\n";
	if(($c_in{'Lipid_ff'})&&($c_in{'Lipid_ff'} == 1)){
		print LEAP_IN "source leaprc.lipid14\n";
	}
	if(($c_in{'Lipid_gaff'})&&($c_in{'Lipid_gaff'} == 1)){
		print LEAP_IN "source ".$curr_FEW_dir."/miscellaneous/GAFFlipid.dat\n";
	}
	print LEAP_IN "loadoff ".$curr_leap_dir."/".$struct_b."_".$chrg_meth.".lib\n";
			
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
	print LEAP_IN "LIG = loadpdb ".$curr_leap_dir."/".$struct_b."_lig.pdb\n";
	print LEAP_IN "REC = loadpdb ".$curr_leap_dir."/".$struct_b."_rec.pdb\n";
	print LEAP_IN "COM = loadpdb ".$curr_leap_dir."/".$struct_b."_com.pdb\n";
	
	if(keys %disulf){
		foreach my $unit_name (@unit_names){
			foreach my $k (keys %disulf){
				print LEAP_IN "bond ".$unit_name.".".$k.".SG ".$unit_name.".".$disulf{$k}.".SG\n";
			}
		}
	}
	print LEAP_IN "loadAmberParams ".$curr_leap_dir."/".$struct_b.".frcmod\n";
	print LEAP_IN "saveAmberParm LIG $cryst_dir/".$struct_b."_vac_lig.top $cryst_dir/".$struct_b."_vac_lig.crd\n";
	print LEAP_IN "savepdb LIG $cryst_dir/".$struct_b."_vac_lig.pdb\n";
	print LEAP_IN "saveAmberParm REC $cryst_dir/".$struct_b."_vac_rec.top $cryst_dir/".$struct_b."_vac_rec.crd\n";
	print LEAP_IN "savepdb REC $cryst_dir/".$struct_b."_vac_rec.pdb\n";
	print LEAP_IN "saveAmberParm COM $cryst_dir/".$struct_b."_vac_com.top $cryst_dir/".$struct_b."_vac_com.crd\n";
	print LEAP_IN "savepdb COM $cryst_dir/".$struct_b."_vac_com.pdb\n";
	print LEAP_IN "charge LIG\n";
	print LEAP_IN "charge REC\n";
	print LEAP_IN "charge COM\n";
	print LEAP_IN "quit\n";
	close LEAP_IN;
}


# Check charge and warnings in leap-log file and
# return reference to hash containing the total charge of 
# complex, receptor,and ligand
sub check_charge_and_warnings{
	my $cryst_dir = shift;
	my $warn = 0;
	my %charge;
	my $no_charge_calls = 0;
	my $f_log_file = "$cryst_dir/leap_1.log";

	open(LOG, $f_log_file) || die "Cannot open log-file $f_log_file\n";

	my $line_count_after_warning = 0;
	my $detected_warning = 0;
	while(my $log_line = <LOG>){
		chomp($log_line);
		if($detected_warning == 1){
			$line_count_after_warning++;
		}
		if($line_count_after_warning == 2){
			$detected_warning = 0;
			$line_count_after_warning = 0;
		}
		if(($detected_warning == 1)&&($line_count_after_warning == 1)){
			if(($log_line =~ m/charge/)||($log_line =~ m/check/)){
				next;
			}
			elsif($log_line =~ m/name change in pdb file residue/){
				next;
			}
			elsif($log_line =~ m/residue name to PDB format/){
				next;
			}
			elsif($log_line =~ m/Exiting LEaP: Errors/){
				next;
			}
			else{
				$warn = 1;
			}
			# Reset for next occurence
			$detected_warning = 0;
			$line_count_after_warning = 0;
		}
		if(($log_line =~ m/Warning/)||($log_line =~ m/WARNING/)){
			$detected_warning = 1;
			if($log_line =~ m/charge/){
				next;
			}
			elsif($log_line =~ m/name change in pdb file residue/){
				next;
			}
			elsif($log_line =~ m/residue name to PDB format/){
				next;
			}
			elsif($log_line =~ m/teLeap: Warning!/){
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
			$no_charge_calls++;
			my @charge = split(/\:\s+/, $log_line);

			if($no_charge_calls == 1){
				$charge{'lig'} = $charge[1];
			}
			if($no_charge_calls == 2){
				$charge{'rec'} = $charge[1];
			}
			if($no_charge_calls == 3){
				$charge{'com'} = $charge[1];
			}
		}
	}
	if($warn == 1){
		print "LEaP run terminated with warnings. Please inspect the LEaP log file at ".$f_log_file."\n";
	}
	if($warn != 1){
		print "Terminated LEaP without warnings.\n";
	}
	return \%charge;
}


# Modify libray file to ensure that the correct number of ions is added
# and small rounding errors do not lead to an addion of a smaller number
# of ions than required.
sub change_charge_lib{
	my $library = shift;
	my $diff = shift;
	my $entry = 0;
	my $new_lib = substr($library, 0, -4);
	$new_lib = $new_lib . "_mod.lib";

	open(LIB, $library) || die "Cannot open Library file $library.\n";
	open(NEW_LIB, ">$new_lib") || die "Cannot open new Library file $new_lib.\n";

	while(my $lib_line = <LIB>){
		chomp($lib_line);

		if($entry == 1){
			my @lib_line = split(/\s+/, $lib_line);
			$lib_line[$#lib_line] = $lib_line[$#lib_line]+$diff;

			foreach (0 ..($#lib_line-1)){
				print NEW_LIB $lib_line[$_] . " ";
			}
			print NEW_LIB $lib_line[$#lib_line]."\n";

			$entry = 0;
			next;
		}


		if($lib_line =~ m/unit\.atoms table/){
			$entry = 1;
		}

		print NEW_LIB $lib_line."\n";
	}
	close LIB;
	close NEW_LIB;

	return $new_lib;
}


# Write input script for second leap run for setup of solvated system
sub create_leap_in_setup_solv{
	my $cryst_dir = shift;
	my $library = shift;
	my $curr_leap_dir = shift;
	my $struct_b = shift;
	my $ref_charge = shift;
	my %charge = %{$ref_charge};
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my %disulf;
	my @tags;
	
	push(@tags, "COM");

	open(LEAP_IN, ">$cryst_dir/leap_script_2.in") || die "Cannot open $cryst_dir/leap_script_2.in";
	
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
	        #print LEAP_IN "source leaprc.water.tip3p\n"; 

	}
	print LEAP_IN "source leaprc.gaff\n";
	print LEAP_IN "loadoff ".$library."\n";
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
	print LEAP_IN "COM = loadpdb ".$curr_leap_dir."/".$struct_b."_com.pdb\n";

	if($c_in{'traj'} == 3){
		print LEAP_IN "REC = loadpdb ".$curr_leap_dir."/".$struct_b."_rec.pdb\n";
		print LEAP_IN "LIG = loadpdb ".$curr_leap_dir."/".$struct_b."_lig.pdb\n";
		push(@tags, "REC");
	}
	print LEAP_IN "loadAmberParams ". $curr_leap_dir."/".$struct_b.".frcmod\n";

	if(keys %disulf){
		foreach my $tag (@tags){
			foreach my $k (keys %disulf){
				print LEAP_IN "bond ".$tag.".".$k.".SG ".$tag.".".$disulf{$k}.".SG\n";
			}
		}
	}

	my $ions = define_counter_ions(\%charge, "com");
	if($ions ne ""){
		print LEAP_IN "addIons COM $ions\n";
	}
	print LEAP_IN "saveAmberParm COM ".$cryst_dir."/".$struct_b."_cio_com.top ".$cryst_dir."/".$struct_b."_cio_com.crd\n";
        if($c_in{'water_model'} eq "TIP3P"){
		print LEAP_IN "solvatebox COM TIP3PBOX 11.0\n";
	}
	if($c_in{'water_model'} eq "OPC"){
		print LEAP_IN "solvatebox COM OPCBOX 11.0\n";
	}
	print LEAP_IN "savepdb COM ".$cryst_dir."/".$struct_b."_solv_com.pdb\n";
	print LEAP_IN "saveAmberParm COM ".$cryst_dir."/".$struct_b."_solv_com.top ".$cryst_dir."/".$struct_b."_solv_com.crd\n";
	if($c_in{'traj'} == 3){
		my $ions = define_counter_ions(\%charge, "rec");
		if($ions ne ""){
			print LEAP_IN "addIons REC $ions\n";
		}
		my $ions = define_counter_ions(\%charge, "lig");
		if($ions ne ""){
			print LEAP_IN "addIons LIG $ions\n";
		}
		print LEAP_IN "saveAmberParm REC ".$cryst_dir."/".$struct_b."_cio_rec.top ".$cryst_dir."/".$struct_b."_cio_rec.crd\n";
		print LEAP_IN "saveAmberParm LIG ".$cryst_dir."/".$struct_b."_cio_lig.top ".$cryst_dir."/".$struct_b."_cio_lig.crd\n";
		if($c_in{'water_model'} eq "TIP3P"){
			print LEAP_IN "solvatebox REC TIP3PBOX 11.0\n";
			print LEAP_IN "solvatebox LIG TIP3PBOX 11.0\n";
		}
		if($c_in{'water_model'} eq "OPC"){
			print LEAP_IN "solvatebox REC OPCBOX 11.0\n";
			print LEAP_IN "solvatebox LIG OPCBOX 11.0\n";
		}
		print LEAP_IN "saveAmberParm REC ".$cryst_dir."/".$struct_b."_solv_rec.top ".$cryst_dir."/".$struct_b."_solv_rec.crd\n";
		print LEAP_IN "savepdb REC ".$cryst_dir."/".$struct_b."_solv_rec.pdb\n";
		print LEAP_IN "saveAmberParm LIG ".$cryst_dir."/".$struct_b."_solv_lig.top ".$cryst_dir."/".$struct_b."_solv_lig.crd\n";
		print LEAP_IN "savepdb LIG ".$cryst_dir."/".$struct_b."_solv_lig.pdb\n";
	}
	print LEAP_IN "quit";
}


# Write input script for second leap run for setup of solvated explicit membrane system
sub create_mem_leap_in_setup_solv{
	my $cryst_dir = shift;
	my $library = shift;
	my $curr_leap_dir = shift;
	my $struct_b = shift;
	my $curr_FEW_dir = shift;
	my $ref_charge = shift;
	my %charge = %{$ref_charge};
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my %disulf;
	my @tags;
	
	push(@tags, "COM");

	open(LEAP_IN, ">$cryst_dir/leap_script_2.in") || die "Cannot open $cryst_dir/leap_script_2.in";

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
		print LEAP_IN "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff10\n";
		print LEAP_IN "loadAmberParams frcmod.ionsjc_tip3p\n";
	}
	print LEAP_IN "source leaprc.gaff\n";
	if(($c_in{'Lipid_ff'})&&($c_in{'Lipid_ff'} == 1)){
		print LEAP_IN "source leaprc.lipid14\n";
	}
	if(($c_in{'Lipid_gaff'})&&($c_in{'Lipid_gaff'} == 1)){
		print LEAP_IN "source ".$curr_FEW_dir."/miscellaneous/GAFFlipid.dat\n";
	}
	print LEAP_IN "loadoff ".$library."\n";
	if(($c_in{'add_lib'})&&(-e $c_in{'add_lib'})){
		print LEAP_IN "loadoff $c_in{'add_lib'}\n";
	}
	if(($c_in{'add_frcmod'})&&(-e $c_in{'add_frcmod'})){
		print LEAP_IN "loadAmberParams $c_in{'add_frcmod'}\n";
	}
	if(($c_in{'cys_f'} ne "")&&(-e $c_in{'cys_f'})){
		my $disulf = read_disulf(\%c_in);
		%disulf = %{$disulf};
	}
	print LEAP_IN "COM = loadpdb ".$curr_leap_dir."/".$struct_b."_com.pdb\n";
	print LEAP_IN "loadAmberParams ". $curr_leap_dir."/".$struct_b.".frcmod\n";

	if(keys %disulf){
		foreach my $tag (@tags){
			foreach my $k (keys %disulf){
				print LEAP_IN "bond ".$tag.".".$k.".SG ".$tag.".".$disulf{$k}.".SG\n";
			}
		}
	}

	print LEAP_IN "\n";
	print LEAP_IN "setBox COM vdw\n";
	print LEAP_IN "savepdb COM ".$cryst_dir."/".$struct_b."_solv_com.pdb\n";
	print LEAP_IN "saveAmberParm COM ".$cryst_dir."/".$struct_b."_solv_com.top ".$cryst_dir."/".$struct_b."_solv_com.crd\n";
	print LEAP_IN "quit";
}


# Read disulfide bridge definition
sub read_disulf{
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my %disulf;

	open(CYS, $c_in{'cys_f'}) || die "Cannot open disulfide bridge definition file".$c_in{'cys_f'}.".\n";

	while(my $cys_line = <CYS>){
		chomp($cys_line);
		if($cys_line =~ m/(\d+)/){
			 my @cys_line = split(/\t/, $cys_line);
			 $disulf{$cys_line[0]} = $cys_line[1];
		}
	}
	close CYS;
	return \%disulf;
}


# Define counter Ions
sub define_counter_ions{
	my $ref_charge = shift;
	my %charge = %{$ref_charge};
	my $tag = shift;
	my $ions;

	if($charge{$tag} > 0.05){
		$ions = "Cl- 0";
	}
	elsif($charge{$tag} < -0.05){
		$ions = "Na+ 0";
	}
	elsif(($charge{$tag} >= -0.05)&&($charge{$tag} <= 0.05)){
		$ions = "";
	}
	return $ions;
}


1;
