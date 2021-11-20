#/usr/bin/perl -w

# Perl module with main functions needed by at least two of the program components
# WAMM, TIW, and LIE for preparation of the molecular dynamics simulations.
use global;
use prepare_params;
use prepare_top_crd_MD;
use setup_MD;
use check_input;
use ReadBackwards;
use File::Copy;
use strict;


# Subroutine for preparation of leap input files, including charge
# calculation, generation of pdb-files of complex, receptor, and ligand,
# and creation of parameter and library files for the ligand.
sub prepare_leap_input{

	my $struct_path = shift;
	my $struct_files_ref = shift;
	my $struct_names_ref = shift;
	my $struct_b_ref = shift;
	my $c_in_ref = shift;
	
	my @struct_files = @{$struct_files_ref};
	my %struct_names = %{$struct_names_ref};
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$c_in_ref};
	my $leap_dir;	# directory for storing leap-input files

	# Create Leap directory
	if($c_in{'r_1'} != 1){
		$leap_dir = $c_in{'root_path'}."/leap";
       	if(! -e $leap_dir){
			mkdir($leap_dir);
       	}
	}

	# Define charge calculation method
	my $chrg_meth;
	if($c_in{'am1'} == 1){
		$chrg_meth = "am1";
	}
	elsif($c_in{'resp'} == 1){
		$chrg_meth = "resp";
	}
		
	#############################
	#  Total charge definition  #
	#############################
	my $charge;
	my $multi;

	if($c_in{'chg'} == 1){
		($charge, $multi) = define_charge_multi(\%struct_names, \%c_in);
	}

	################################
	#  Calculation of AM1-Charges  #
	################################

	if(($c_in{'am1'} == 1)&&($c_in{'calc_c'}==1)){

		foreach my $s (@struct_files){
			
			# Create ligand specific folders in "leap" directory
			if(! -e "$leap_dir/$struct_b{$s}"){
				mkdir($leap_dir."/".$struct_b{$s});
			}
				
			# Calculate AM1-BCC charges
			my $mol2_out = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_am1.mol2";
			my $am1_log = $c_in{'root_path'}."/antechamber_am1.log";
			my $am1_command = "antechamber -i " . $s . " -fi mol2 ";
			$am1_command .= "-o $mol2_out -fo mol2 ";
			$am1_command .= "-c bcc";
			if(($c_in{'chg'} == 1)&&(($charge != 0)&&($multi != 0))){
				$am1_command = charge_multi_ante($c_in{'chg'}, $am1_command, $struct_names{$s}, $charge, $multi);
			}
			$am1_command .= " > $am1_log 2>&1";

			print "\nCalculating AM1-BCC charges for ligand ".$struct_b{$s}.".\n";
			call_prog($am1_command, $am1_log);
			unlink($am1_log);
			
			check_set_lig_resname_mol2($mol2_out);
		}
	}

	#################################
	#  Calculation of RESP-charges  #
	#################################

	if(($c_in{'resp'} == 1)&&($c_in{'calc_c'} == 1)){

		if($c_in{'r_1'} == 1){
			
			# Create directory for gaussian calculation
			if(! -e $c_in{'root_path'}."/gauss"){
				mkdir($c_in{'root_path'}."/gauss");
			}

			# Set default path for gaussian batch script
			if((! exists $c_in{'pbs_w'})||($c_in{'pbs_w'} eq "")){
				$c_in{'pbs_w'} = $c_in{'root_path'};
			}

			# Generate file for automatic submission of gaussian jobs
			my $gauss_submit = $c_in{'root_path'}."/gauss/submit.sh";
			if(($c_in{'pbs_g'}==1)&&($c_in{'pbs_t'})&&($c_in{'pbs_w'})){
				open(SUBMIT, ">".$gauss_submit);
			}

			foreach my $s (@struct_files){
				
				# Create ligand specific folders in "gauss" directory
				if(! -e $c_in{'root_path'}."/gauss/".$struct_b{$s}){
					mkdir($c_in{'root_path'}."/gauss/".$struct_b{$s});
				}
				
				# Create input file for gaussian calculation
				my $gauss_log = $c_in{'root_path'}."/antechamber_gauss_setup.log";
				my $gauss_command = "antechamber -i " . $s . " -fi mol2 ";
				$gauss_command .= "-o ".$c_in{'root_path'}."/gauss/".$struct_b{$s}."/".$struct_b{$s}.".gcrt -fo gcrt ";
				$gauss_command .= "-ch $struct_b{$s}.chk -gv 0";
				if(($c_in{'chg'} == 1)&&(($charge != 0)&&($multi != 0))){
					$gauss_command = charge_multi_ante($c_in{'chg'}, $gauss_command, $struct_names{$s}, $charge, $multi);
				}
				$gauss_command .= " > $gauss_log 2>&1";

				print "Generating gaussian input files for ligand ".$struct_b{$s}.".\n";
				call_prog($gauss_command, $gauss_log);
				unlink($gauss_log);

				if(($c_in{'pbs_g'}==1)&&($c_in{'pbs_t'})&&($c_in{'pbs_w'})){
					generate_pbs_gauss($struct_b{$s}, \%c_in);

					# Write entry to submission file
					print SUBMIT "qsub " . $c_in{'pbs_w'} . "/gauss/" . $struct_b{$s} . "/" . $struct_b{$s} . ".pbs\n";
				}
			}
			close SUBMIT;
			chmod 0755, $gauss_submit;
		}

		# Processing output of gaussian calculation
		if($c_in{'r_2'} == 1){

			foreach my $s (@struct_files){
				
				# Create ligand specific folders in "leap" directory
				if(! -e "$leap_dir/$struct_b{$s}"){
					mkdir($leap_dir."/".$struct_b{$s});
				}

				# Check gaussian output
				my $gout_file = $c_in{'root_path'}."/gauss/".$struct_b{$s}."/".$struct_b{$s}.".gout";
				
				if(! -e $gout_file){
					print "The gaussian output file $gout_file \n";
					print "required for RESP-charge calculation does not exist. Please ensure that\n";
					print "the file is present in the correct location under ".$gout_file."\n";
					print "and re-start the procedure.\n";
					exit;
				}
				else{
					my $check_normal_term = 0;
					
					my $gout_backwards = File::ReadBackwards->new($gout_file) || "Cannot open file $gout_file for reading.\n";
					my $last_gout_line = $gout_backwards->readline;

					if($last_gout_line =~ m/Normal termination/){
						$check_normal_term = 1;
					}
					
					if($check_normal_term != 1){
						print "Gaussian calculation was not terminated correctly for ligand $struct_b{$s}.\n";
						print "Please make sure that the ESP-calculation with gaussian was normally terminated\n";
						print "and re-start the procedure.\n";
						exit;
					}
				}

				# Check hydrogen bond presence
				check_hbond_gauss($gout_file, $struct_b{$s});

				# Generate mol2-file
				my $r_2_log = $c_in{'root_path'}."/antechamber_r2.log";
				my $r_2_com = "antechamber -i $gout_file -fi gout ";
				$r_2_com .= "-o ".$leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_resp.mol2 -fo mol2 -c resp";
				if(($c_in{'chg'} == 1)&&(($charge != 0)&&($multi != 0))){
					$r_2_com = charge_multi_ante($c_in{'chg'}, $r_2_com, $struct_names{$s}, $charge, $multi);
				}
				$r_2_com .= " > $r_2_log 2>&1";

				print "Calculating RESP charges for ligand ".$struct_b{$s}.".\n";
				call_prog($r_2_com, $r_2_log);
				unlink($r_2_log);

				# Use atom names from original mol2-file to ensure atom name matching
				set_atom_names($struct_path, $leap_dir, $struct_b{$s});
			}
		}
	}

	# Delete temporary antechamber files
	unlink glob "ANTECHAMBER*";
	unlink "ATOMTYPE.INF";
	unlink glob "sqm*";
	

	if(($c_in{'am1'} == 1) || ($c_in{'r_2'} == 1)){

		###############################
		#  Average charges of isomers #
		###############################
		my %avg_charge;
	
		if(-e $c_in{'avg_chg'}){
			my $ref_avg_charge = read_isomer_pairs($c_in{'avg_chg'});
			%avg_charge = %{$ref_avg_charge};			
			average_charges($leap_dir, \%avg_charge, \%c_in);
		}

		##############################
		# Generation of frcmod-files #
		##############################

		foreach my $s (@struct_files){

			if((($c_in{'am1'} == 1)&&(-e $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_am1.mol2"))
			||(($c_in{'resp'} == 1)&&(-e $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_resp.mol2"))){

				my $frc_file = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}.".frcmod";

				if(! -e $frc_file){
					my $frc_log = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_frcmod.log";
					my $frc_command = "parmchk2 -i ".$leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_".$chrg_meth.".mol2 ";
					$frc_command .= "-f mol2 -o ".$frc_file;
					$frc_command .= " 2> ".$frc_log;

					print "\nSetting up parameter file for ligand ".$struct_b{$s}.".\n";
					call_prog($frc_command, $frc_log);
					unlink($frc_log);
					my $check_frcfile = check_frcmod($frc_file);
					if($check_frcfile == 0){
						exit;
					}
				}
			}
		}

		############################################
		# Setup of receptor and complex structures #
		############################################

		foreach my $s (@struct_files){
			
			print "\nGenerating structures for ".$struct_b{$s}.".\n";
			
			# Create ligand specific folders in "leap" directory - In case no charge calculation was requested before
			if(! -e "$leap_dir/$struct_b{$s}"){
				mkdir($leap_dir."/".$struct_b{$s});
			}

			# Copy receptor file to leap input file of the corresponding structure
			my $rec_file = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_rec.pdb";

			# If receptor file in leap directory already exists,
			# provide information to user
			if(-e $rec_file){
				my $count = 1;
				while(-e $rec_file.".setup".$count){
					$count++;
				}
				print "\nWARNING:\nThe existing receptor file $rec_file\n";
				print "will be re-named to ".$rec_file.".setup".$count.".\n"; 
				print "If you want to use this file for simulation setup later,\n";
				print "please re-change the name into $rec_file.\n";
				rename($rec_file, $rec_file.".setup".$count);
			}
			
			if(($c_in{'rec'})&&(-e $c_in{'rec'})){
				copy_receptor($c_in{'rec'}, $rec_file);
			}
			else{
				print "\nERROR: Receptor file was not correctly specified. Please ensure\n";
				print "that at least one receptor file is defined in the command file.\n";
				exit;
			}
			

			# Generate complex
			###################

			# Prepare Ligand-file - Gereate ligand file in PDB format
			my $mol2_lig_file = $struct_path."/".$struct_b{$s}.".mol2"; # Original structure file (mol2 format)
			my $pdb_lig_file = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_lig.pdb"; # PDB-file of ligand

			if(-e $pdb_lig_file){
				if(-e $pdb_lig_file){
					my $count = 1;
					while(-e $pdb_lig_file.".setup".$count){
						$count++;
					}
					print "\nWARNING:\nThe existing ligand file $pdb_lig_file\n";
					print "will be re-named to $pdb_lig_file.setup".$count.".\n";
					print "If you want to use this file for simulation setup later,\n";
					print "please re-change the name into $pdb_lig_file.\n";
					rename($pdb_lig_file, $pdb_lig_file.".setup".$count);
				}
			}
			
			my $pdb_conversion_log = $leap_dir."/".$struct_b{$s}."/leap_pdb_conversion_lig.log";
			my $pdb_command = "antechamber -i ".$mol2_lig_file." -fi mol2 ";
			$pdb_command .= "-o ".$pdb_lig_file." -fo pdb";
			$pdb_command .= " > ".$pdb_conversion_log." 2>&1";
			call_prog($pdb_command, $pdb_conversion_log);
			unlink($pdb_conversion_log);
			
			# Ensure that residue name of ligand is uniform
			check_set_residue_name_ligand($pdb_lig_file);

			# Create complex from Ligand and Receptor files
			my $compl_file = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_com.pdb";

			if(-e $compl_file){
				my $count = 1;
				while(-e $compl_file.".setup".$count){
					$count++;
				}
				print "\nWARNING:\nThe existing complex file $compl_file\n";
				print "will be re-named to ".$compl_file.".setup".$count.".\n";
				print "If you want to use this file for simulation setup later,\n";
				print "please re-change the name into $compl_file.\n\n";
				rename($compl_file, $compl_file.".setup".$count);
			}
			
			if(($c_in{'rec_compl'})&&(-e $c_in{'rec_compl'})){
				$rec_file = $c_in{'rec_compl'};
			}
			
			create_complex($compl_file, $pdb_lig_file, $rec_file, \%c_in);
			
			# If 'membrane' PDB file is present, add membrane information
			# to receptor and complex
			if((($c_in{'membrane'})&&($c_in{'membrane'} ne ""))&&(($c_in{'withMem'})&&($c_in{'withMem'} == 1))){
				if(-e $c_in{'membrane'}){
					add_membrane($rec_file, $compl_file, $c_in{'membrane'}, $c_in{'lig_cutoff'}, $pdb_lig_file);
				}
			}
		}


		#################################
		#  Generation of Library files  #
		#################################

		foreach my $s (@struct_files){

			if((($c_in{'am1'} == 1)&&(-e $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_am1.mol2"))
			||(($c_in{'resp'} == 1)&&($c_in{'r_2'} == 1)&&(-e $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_resp.mol2"))){

				my $lib_leap_in = "$leap_dir/$struct_b{$s}/$struct_b{$s}_".$chrg_meth."_leap.in";
				my $lib_file = "$leap_dir/$struct_b{$s}/$struct_b{$s}_$chrg_meth.lib";

				if(-e $lib_file){
					unlink $lib_file;
				}
				setup_leap_for_lib($lib_leap_in, $lib_file, $leap_dir, $chrg_meth, $s, \%c_in, \%struct_b);

				# Run leap
				my $lib_log = $leap_dir."/".$struct_b{$s}."/".$struct_b{$s}."_".$chrg_meth."_leap.log";
				my $lib_leap_c = "tleap -f $lib_leap_in >> $lib_log 2 >> $lib_log";
				print "\nRunning LEaP for lib-file generation for ".$struct_b{$s}.".\n";
				call_prog($lib_leap_c, $lib_log);

				# Check correct performance
				my $log_file = "$leap_dir/$struct_b{$s}/$struct_b{$s}_".$chrg_meth."_leap.log";
				check_leap($log_file);
				unlink("./leap.log");
			}
		}
	} # End 'skip if r_1'
} # end subroutine prepare_leap_input


# Subroutine for copying the receptor structure and ensuring correct
# handling of blocking groups, since as starting from amber18 the
# representation of the ACE and NME blocking groups has changed, but
# the tutorial still uses a receptor file in the old format.
sub copy_receptor{
	my $org_rec_file = shift;
	my $target_rec_file = shift;

	open(PDB, $org_rec_file) || die "Cannot open pdb-file $org_rec_file.\n";
	open(PDB_MOD, ">$target_rec_file") || die "Cannot open pdb-file $target_rec_file for writing.\n";

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
				print PDB_MOD $pdb_line."\n";
			}
		}
		else{
			print PDB_MOD $pdb_line."\n";
		}
	}
	close PDB;
	close PDB_MOD;
}


# Subroutine for generation of coordinated and topology files for
# MD calculations
sub create_crd_top{
	
	my $procedure = shift;
	my $chrg_meth = shift;
	my $md_dir = shift;
	my $curr_FEW_dir = shift;
	my $struct_files_ref = shift;
	my $struct_b_ref = shift;
	my $c_in_ref = shift;
	my @struct_files = @{$struct_files_ref};
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$c_in_ref};
	my $rec_top_crd = 0; # For receptor setup of 3-trajectory approach
	
	foreach my $s (@struct_files){

		# Create structure specific directory in MD-folder
		if(! -e  $md_dir."/".$struct_b{$s}){
			mkdir($md_dir."/".$struct_b{$s});
		}
			
		# Create 'cryst' folder
		my $cryst_dir = $md_dir."/".$struct_b{$s}."/cryst";
			
		if(! -e $cryst_dir){
			mkdir($cryst_dir);
		}
		else{
			my $count = 1;
			while(-e $md_dir."/".$struct_b{$s}."/cryst_setup".$count){
				$count++;
			}
			print "\nWARNING:\n";
			print "A folder ".$cryst_dir." with\n";
			print "structural information already exists. This is correct, if you ran\n";
			print "a MD setup using the same 'root_path' before.\n";
			print "The toplogy and coordinate files of the previous setup will be kept\n";
			print "in the folder ".$md_dir."/".$struct_b{$s}."/cryst_setup".$count."\n";
			print "from now on.\n";
			move($cryst_dir, $md_dir."/".$struct_b{$s}."/cryst_setup".$count);
			mkdir($cryst_dir);
		}
			
		# Set path to leap directory, since this was not defined globally
		my $curr_leap_dir = $c_in{'root_path'}."/leap/".$struct_b{$s};

		# Re-name CYS to CYX for residues with S-S-bond connectivity
		if(($c_in{'cys_f'} ne "")&&(-e $c_in{'cys_f'})){
			my $complex_file = $curr_leap_dir."/".$struct_b{$s}."_com.pdb";
			my $rec_file = $curr_leap_dir."/".$struct_b{$s}."_rec.pdb";
		
			my $disulf_ref = read_disulf(\%c_in);
			my %disulf = %{$disulf_ref};
			exchange_cys_to_cyx($complex_file, $c_in{'rec_res'}, \%disulf);
			exchange_cys_to_cyx($rec_file, $c_in{'rec_res'}, \%disulf);
		}

		# Create first leap-input script for setup of raw topology and coordinate files
		if((! $c_in{'withMem'})||($c_in{'withMem'} == 0)){
			create_leap_in_top_crd($cryst_dir, $curr_leap_dir, $struct_b{$s}, $chrg_meth, \%c_in);
		}
		else{
			create_mem_leap_in_top_crd($cryst_dir, $curr_leap_dir, $struct_b{$s}, $chrg_meth, $curr_FEW_dir, \%c_in);
		}

		# First leap call
		if(-e $cryst_dir."/leap_1.log"){
			unlink($cryst_dir."/leap_1.log");
		}

		# Status message
		print "\nRunning LEaP to setup raw topology and coordinate files for ".$struct_b{$s}.".\n";

		my $f_leap_log = $cryst_dir."/leap_1.log";
		my $f_leap = "tleap -f ".$cryst_dir."/leap_script_1.in > ".$f_leap_log." 2>&1";
		call_prog($f_leap, $f_leap_log);

		# Check charge and warnings in leap-log file
		my $ref_charge = check_charge_and_warnings($cryst_dir);
		my %charge = %{$ref_charge};
		
		# Check if created topology files exist and are not empty
		foreach my $tag ("com", "rec", "lig"){
			my $topo_vac_file = $cryst_dir."/".$struct_b{$s}."_vac_".$tag.".top";
			if((! -e $topo_vac_file)||(-z $topo_vac_file)){
				print "ERROR: At least one topology file was not correctly created. Please look\n";
				print "in the file $cryst_dir/leap_1.log for a detailed error messarge.\n";
				exit;
			}
		}		

		# Modify library of ligand if charge is marginally smaller than next integer value
		my $int = int($charge{'lig'});
		my $library = $curr_leap_dir."/".$struct_b{$s}."_".$chrg_meth.".lib";
		my $diff = 0;

		# Case charge slightly smaller than next integer
		if(((abs($charge{'lig'}))-$int)>0.9){
			if($charge{'lig'}>0){
				$diff = (($int+1)-$charge{'lig'});
			}
			if($charge{'lig'}<0){
				$diff = (($int-1)+$charge{'lig'});
			}
		}

		# Case charge slightly larger than next integer
		if(((abs($charge{'lig'}))-$int)<0.005){
			if($charge{'lig'}<0){
				$diff = abs($charge{'lig'});
			}
			if($charge{'lig'}>0){
				$diff = (-1)*$diff;
			}
		}

		if((abs($diff)<0.005)&&($diff != 0)){
			$library = change_charge_lib($library, $diff);
		}


		# Prepare input script for second leap run for setup of solvated system
		if((! $c_in{'withMem'})||($c_in{'withMem'} == 0)){
			create_leap_in_setup_solv($cryst_dir, $library, $curr_leap_dir, $struct_b{$s}, \%charge, \%c_in);
		}
		else{
			create_mem_leap_in_setup_solv($cryst_dir, $library, $curr_leap_dir, $struct_b{$s}, $curr_FEW_dir, \%charge, \%c_in);
		}
			

		# Second Leap call
		my $s_leap_log = $cryst_dir."/leap_2.log";
		my $s_leap = "tleap -f $cryst_dir/leap_script_2.in > ".$s_leap_log." 2>&1";
		call_prog($s_leap, $s_leap_log);

		# Check warnings
		print "\nRunning LEaP to setup topology and coordinate files with solvent for $struct_b{$s}.\n";
		check_leap("$cryst_dir/leap_2.log");
		
		# Check if created topology files exist and are not empty
		my @tags;
		if($c_in{'traj'} == 1){
			@tags = ("com");
		}
		if($c_in{'traj'} == 3){
			@tags = ("com", "rec", "lig");
		}
		foreach my $tag (@tags){
			my $topo_vac_file = $cryst_dir."/".$struct_b{$s}."_solv_".$tag.".top";
			if((! -e $topo_vac_file)||(-z $topo_vac_file)){
				print "ERROR: At least one topology file was not correctly created. Please look\n";
				print "in the file $cryst_dir/leap_2.log for a detailed error messarge.\n";
				exit;
			}
		}	
		
		unlink($c_in{'root_path'}."/leap.log");

		# Copy receptor topology and coordinate files to "rec" directory if 3-trajectory approach of WAMM workflow
		if(($c_in{'traj'} == 3)&&($rec_top_crd == 0)&&($procedure eq "wamm")){
			my $rec_cryst_dir = "$md_dir/rec/cryst";
			mkdir($rec_cryst_dir);

			my $source_top = $cryst_dir."/".$struct_b{$s}."_solv_rec.top";
			my $source_crd = $cryst_dir."/".$struct_b{$s}."_solv_rec.crd";
			my $target_top = $rec_cryst_dir."/rec_solv.top";
			my $target_crd = $rec_cryst_dir."/rec_solv.crd";

			copy($source_top, $target_top) || die "Cannot copy $source_top to $target_top\n";
			copy($source_crd, $target_crd) || die "Cannot copy $source_crd to $target_crd\n";

			$rec_top_crd = 1;
		}
	}
}


# Subroutine for the generation of the input files for the equilibration
sub prepare_equilibration{

	my $procedure = shift;
	my $chrg_meth = shift;
	my $md_dir = shift;
	my $struct_files_ref = shift;
	my $struct_b_ref = shift;
	my $c_in_ref = shift;
	
	my @struct_files = @{$struct_files_ref};
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$c_in_ref};

	# Read names of equilibration template files
	my @equi_files = <$c_in{'equi'}/*>;

	# Directory for equilibration
	my $equi_dir;
	
	# Variable for detection of warning messages
	my $warn = 0;
	
	# User status information
	print "\nPreparing input files for equilibration phase of MD simulations.\n";

	# Script for storing pbs-submission commands for equilibration
	if(($c_in{'pbs_e'})&&(-e $c_in{'pbs_e'})){
		open(QSUB_EQUI, ">".$md_dir."/qsub_equi.sh") || die "Cannot open".$md_dir."/qsub_equi.sh\n";
	}

	my @sub_tags = ();
	# One trajectory approach
	if($c_in{'traj'} == 1){
		@sub_tags = ("com");
	}
	# Three trajectory approach
	if($c_in{'traj'} == 3){
		@sub_tags = ("com", "lig");
	}

	foreach my $s (@struct_files){
		foreach my $tag (@sub_tags){
			my $tag_dir = $md_dir."/".$struct_b{$s}."/".$tag;
			if(! -e $tag_dir){
				mkdir($tag_dir);
						
				$equi_dir = $md_dir."/".$struct_b{$s}."/".$tag."/equi";
				if(! -e $equi_dir){
					mkdir($equi_dir);
				}
				
				generate_scripts_equi($tag, $equi_dir, \@equi_files, \%c_in);
				
				if(($c_in{'pbs_e'})&&(-e $c_in{'pbs_e'})){
					generate_pbs_equi($equi_dir, $tag, $struct_b{$s}, $chrg_meth, \%c_in);
					print QSUB_EQUI "qsub ".$c_in{'MD_path'}."/MD_".$chrg_meth."/".$struct_b{$s}."/".$tag."/equi/".$struct_b{$s}.".pbs\n";
				}
			}
			else{
				if($warn == 0){
					my $tag_name;
					my $message;
					if($tag eq "com"){
						$tag_name = "complex";
						$message = "A sub-folder for complex MD simulations already exists.\n";
						$message .= "The files in this folder will not be changed or overwritten. ";
						if(-e $md_dir."/".$struct_b{$s}."/lig"){
							$message = "Sub-folders for complex and ligand MD simulations already exist.\n";
							$message .= "The files in these folders will not be changed or overwritten. ";
							$tag_name .= " and the ligand";
						}
					}
					if($tag eq "lig"){
						$message = "A sub-folder for ligand MD simulations already exists.\n";
						$message .= "The files in this folder will not be changed or overwritten. ";
						$tag_name = "ligand";
					}
			
					print "\nWARNING:\n";
					print $message."This is correct,\n";
					print "if you run a 3-trajectory setup after a 1-trajectory setup. If you want to\n";
					print "generate new MD input files for the $tag_name, please remove or rename the\n";
					print "old MD folder '".$c_in{'root_path'}."/MD_".$chrg_meth."'.\n";
					$warn = 1;
				}
			}
		}
	}

	# Generate the receptor structure input in case of a 3-trajectory approach only once
	if(($c_in{'traj'} == 3)&&($procedure eq "wamm")){
		my $tag = "rec";
		$equi_dir = $md_dir."/".$tag."/equi";
		mkdir($equi_dir);

		# Determine name of file to be copied
		generate_scripts_equi($tag, $equi_dir, \@equi_files, \%c_in);
		
		if(($c_in{'pbs_e'})&&(-e $c_in{'pbs_e'})){
			generate_pbs_equi($equi_dir, $tag, "rec", $chrg_meth, \%c_in);
			print QSUB_EQUI "qsub ".$c_in{'MD_path'}."/MD_".$chrg_meth."/".$tag."/equi/".$tag.".pbs\n";
		}
	}
	if(($c_in{'pbs_e'})&&(-e $c_in{'pbs_e'})){
		close QSUB_EQUI;
		chmod 0755, "$md_dir/qsub_equi.sh";
	}
	return $warn;
}



sub prepare_MD_production{

	my $procedure = shift;
	my $chrg_meth = shift;
	my $md_dir = shift;
	my $warn = shift;
	my $struct_files_ref = shift;
	my $struct_b_ref = shift;
	my $c_in_ref = shift;
	
	my @struct_files = @{$struct_files_ref};
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$c_in_ref};

	# Variables
	my $prod_dir;
	my $total_runs; 	# Total number of MD runs required to reach requested ns
	my @sub_tags = ();
	# One trajectory approach
	if($c_in{'traj'} == 1){
		@sub_tags = ("com");
	}
	# Three trajectory approach
	if($c_in{'traj'} == 3){
		@sub_tags = ("com", "lig");
	}

	# User status information
	print "\nGenerating input files for production phase of MD simulations.\n\n";

	# Generate script for easy qsub-submission
	if(($c_in{'pbs_p'})&&(-e $c_in{'pbs_p'})){
		open(QSUB_MD, ">".$md_dir."/qsub_MD.sh") || die "Cannot open qsub-Script".$md_dir."/qsub_equi.sh for wrinting\n";
	}
	my $rel_MD_dir = "/MD_".$chrg_meth;

	foreach my $s (@struct_files){
		foreach my $tag (@sub_tags){
			my $tag_dir = $md_dir."/".$struct_b{$s}."/".$tag;

			if(! -e $tag_dir){
				print "\nERROR:\n";
				print "In this stage of the workflow a directory $tag_dir should already have been created.\n";
				print "Please ensure that you run the setup of MD equilibration simulations before the setup.\n";
				print "of MD production simulations.\n";
				exit;
			}

			$prod_dir = $md_dir."/".$struct_b{$s}."/".$tag."/prod";
			if(! -e $prod_dir){
				mkdir($prod_dir);

				if(($c_in{'prod'})&&(-e $c_in{'prod'})&&($c_in{'prod_t'})){
					$total_runs = setup_MD_prod_in($prod_dir, \%c_in);
				}

				if(($c_in{'pbs_p'})&&(-e $c_in{'pbs_p'})){
					create_PBS_prod_script($prod_dir, $rel_MD_dir, $struct_b{$s}, $tag, $total_runs, \%c_in);
					print QSUB_MD "qsub ".$c_in{'MD_path'}.$rel_MD_dir."/".$struct_b{$s}."/".$tag."/prod/".$struct_b{$s}.".pbs\n";
				}
			}
			else{
				if($warn == 0){
					my $tag_name;
					if($tag eq "com"){
						$tag_name = "complex";
					}
					if($tag eq "lig"){
						$tag_name = "ligand";
					}
			
					print "\nWARNING:\n";
					print "A sub-folder for $tag_name MD production simulations already exists.\n";
					print "The files in this folder will not be changed or overwritten.\n";
					print "This is correct if you ran a 1-trajectory setup before a \n";
					print "3-trajectory setup.\n";
					print "If you want to generate new MD input files for the $tag_name,\n";
					print "please remove or rename the old MD folder '".$c_in{'root_path'}."/MD_".$chrg_meth."'.\n";
					$warn = 1;
				}
			}
		}
	}

	# Receptor-Preparation - Three trajectory approach
	if(($c_in{'traj'} == 3)&&($procedure eq "wamm")){
		my $tag = "rec";
		$prod_dir = $md_dir."/".$tag."/prod";
		mkdir($prod_dir);

		if(($c_in{'prod'})&&(-e $c_in{'prod'})&&($c_in{'prod_t'})){
			$total_runs = setup_MD_prod_in($prod_dir, \%c_in);
		}

		if(($c_in{'pbs_p'})&&(-e $c_in{'pbs_p'})){
			create_PBS_prod_script($prod_dir, $rel_MD_dir, "rec", $tag, $total_runs, \%c_in);
			print QSUB_MD "qsub ".$c_in{'MD_path'}.$rel_MD_dir."/".$tag."/prod/".$tag.".pbs\n";
		}
	}
	if(($c_in{'pbs_p'})&&(-e $c_in{'pbs_p'})){
		close QSUB_MD;
		chmod 0755, "$md_dir/qsub_MD.sh";
	}
}


# Subroutine for renaming CYS to CYX in PDB file
# at positions indicated for S-S-bonds
sub exchange_cys_to_cyx{
	my $file = shift;
	my $rec_res = shift;
	my $disulf_ref = shift;
	my %disulf = %{$disulf_ref};
	
	my %disulf_cys;
	foreach my $k (keys %disulf){
		$disulf_cys{$k} = 1;
		$disulf_cys{$disulf{$k}} = 1;
	}
	
	my $new_str = "";
	my $current_resno;
	open(PDB_FILE, $file) || die "Cannot open PDB file $file for reading.\n";
	while(my $pdb_l = <PDB_FILE>){
		if(($current_resno < $rec_res)&&(($pdb_l =~ m/ATOM/)||($pdb_l =~ m/HETATM/))){
			my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
			unpack("a6    a5  x   a4     a     a3    x2  a4   x4   a8  a8  a8", $pdb_l);
		
			$resno =~ s/^\s+//g;
			$resno =~ s/\s+$//g;
			$resname =~ s/^\s+//g;
			$resname =~ s/\s+$//g;
			$atom =~ s/^\s+//g;
			$atom =~ s/\s+$//g;
		
			$current_resno = $resno;
		
			if(($disulf_cys{$resno})&&($resname eq "CYS")){
				$pdb_l =~ s/CYS/CYX/;
			}
			if(($disulf_cys{$resno})&&(($resname eq "CYS")||($resname eq "CYX"))){
				# If HG is still present in structure this will be removed
				if($atom eq "HG"){
					print "\nATTENTION: The atom HG of residue $resno present in the provided\n";
					print "receptor structure will be removed to enable a disulfide bond\n";
					print "connection requested for this residue.\n";
					next;
				}
			}
			if(($disulf_cys{$resno})&&(($resname ne "CYS")&&($resname ne "CYX"))){
				print "\nERROR: Residue $resno specified as a residue participating in\n";
				print "a disulfide bond is not a cysteine in the provided PDB file.\n";
				print "Please ensure that the residue numbers given in the file for\n";
				print "the definition of the disulfide bonds are correct.\n";
			}
		}
		$new_str .= $pdb_l;
	}
	close PDB_FILE;
	
	# Write new file with new residue naming
	my $new_pdb_file = $file."_tmp";
	open(NEW_PDB, ">".$new_pdb_file) || die "Cannot open file $new_pdb_file for writing.\n";
	print NEW_PDB $new_str;
	close NEW_PDB;
	
	unlink($file);
	rename($new_pdb_file, $file);
}


1;
