#/usr/bin/perl -w

# Perl module for the preparation of structures for equilibration and for the determination
# of binding energies with the thermodynamic integration approach.

use File::Copy;
use File::Path;
use strict;

use lib "../additional_libs";
use global;
use check_input;
use separate_structures;
use common_prepare_MD;
use prepare_top_crd_MD;
use setup_TI_MD;
use atom_matching;
use ReadBackwards;	# External perl module

sub TIW{

	my $current_dir = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	# Skip all other parts, if only ddG calculation shall be conducted.
	if($c_in{'ti_ddG'}){
		goto TI_DDG;
	}

	#################################
	# Specific check for TIW module #
	#################################
	
	if($c_in{'ti_sim'} == 1){
		my $check_tiw = &check_input_tiw(\%c_in);
		if($check_tiw != 1){
			exit;
		}
	}
	

	#########################
	#  Convert sdf to mol2  #
	#########################

	if($c_in{'sdf'} == 1){
	
		# Check if structure file was provided with extension
		my $ext = substr($c_in{'fi_lig'}, -4);
		my $fi_lig;
		if($ext eq ".sdf"){
			$fi_lig = substr($c_in{'fi_lig'}, 0, -4);
		}
		else{
			$fi_lig = $c_in{'fi_lig'};
		}
		
		# Coversion of sdf to mol2
		my $i_file_base = $c_in{'i_path'}."/".$fi_lig;
		my $conv_log = $c_in{'root_path'}."/babel.log";
		my $conv_c = "babel -isdf ".$i_file_base.".sdf -omol2 ".$i_file_base.".mol2 > ".$conv_log." 2>&1";
		call_prog($conv_c, $conv_log);
		if(-s $conv_log){
			print "\n\nERROR: Problem in Babel conversion of sdf-format to mol2-format. Please consult the file\n";
			print "$conv_log for details.\n\n";
			exit;
		}
		else{
			unlink($conv_log);
		}
	}


	##########################
	#  Structure Separation  #
	##########################

	if($c_in{'sep'} == 1){
		separate(\%c_in);
	}


	#################################
	#  Define paths and structures  #
	#################################

	# Determine folder in which the structures of the ligands (in mol2-format) are located
	my $struct_path;

	# Case 1: Structures automatically deposited in "structs" folder of root_path
	if($c_in{'sep'} == 1){
		$struct_path = $c_in{'root_path'} . "/structs";
	}
	# Case 2: Structures present in input-folder with path "i_path"
	else{
		$struct_path = $c_in{'i_path'};
	}

	# Get names of all files present in struct-directory
	my @struct_files = <$struct_path/*>;    # Complete path to structure files
	my %struct_names;   # Names of structure files
	my %struct_b;       # Basenames of structure files

	foreach my $file (@struct_files){
		my @struct_path_split = split(/\//, $file);
		$struct_names{$file} = $struct_path_split[$#struct_path_split];
		my $struct_b = substr($struct_path_split[$#struct_path_split], 0, -5);
    	$struct_b{$file} = $struct_b;
	}

	# Set MD_path to root_path if no specific MD_path is provided
	if($c_in{'MD_path'} eq ""){
		$c_in{'MD_path'} = $c_in{'root_path'};
	}


	###################
	# Check mol2 file #
	###################
	
	foreach my $s (@struct_files){
		# Check atom and residue naming in provided mol2 file
		my $check_mol = check_mol2($s, "TIW", \%c_in);
		# End program if check failed
		if($check_mol != 1){
			exit;
		}
	}

	###################################
	# Preparation of LEAP input files #
	###################################

	if((exists $c_in{'leap'})&&($c_in{'leap'} == 1)){
		prepare_leap_input($struct_path, \@struct_files, \%struct_names, \%struct_b, \%c_in);
		
		if((! exists $c_in{'r_1'})&&($c_in{'r_1'} != 1)){
			print "\nPlease note:\n";
			print "Parameters have been generated automatically by antechamber without checking\n";
			print "their plausibility. Parameter quality can be analyzed by high level quantum mechanics\n";
			print "calculations, e.g. for checking of torsion parameters and normal modes.\n\n";
		}
	}


	##################################
	#  Setup of MD simulation input  #
	##################################

	if($c_in{'sim'} == 1){
	
		my $warnings = 0;
		my @chrg_meth;

		if($c_in{'P_am1'} == 1){
			push(@chrg_meth, "am1");
		}
		if($c_in{'P_resp'} == 1){
			push(@chrg_meth, "resp");
		}

		foreach my $chrg_meth (@chrg_meth){
		
			#########################################################
			# Check existence of structure specific parameter files #
			#########################################################
	
			my $check_dir_exist = check_file_existence_for_commonMDsetup($chrg_meth, \%c_in, \%struct_b);
			# End program if check failed
			if($check_dir_exist != 1){
				exit;
			}

			# Generate MD-folder
			my $md_dir = $c_in{'root_path'}."/MD_".$chrg_meth;
		
			if(! -e $md_dir){
				mkdir($md_dir);
			}

			############################################
			#  Generate topology and coordinate files  #
			############################################

			create_crd_top("tiw", $chrg_meth, $md_dir, $current_dir, \@struct_files, \%struct_b, \%c_in);

			####################################
			#  Create scripts for equilbration #
			####################################

			# Setup folders, input scripts and pbs script for equilibration
			if(($c_in{'equi'})&&(-e $c_in{'equi'})){
				$warnings = prepare_equilibration("tiw", $chrg_meth, $md_dir, \@struct_files, \%struct_b, \%c_in);
			}
		}
	}

	
	########################
	# Setup TI simulations #
	########################
	
	if($c_in{'ti_sim'} == 1){

		# Create directories
		######################
		
		# TI directory
		my $ti_dir = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'};
		if(! -e $ti_dir){
			mkdir($ti_dir);
		}
		
		# Transformation directory
		my $transform_dir = $ti_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'};
		if(! -e $transform_dir){
			mkdir($transform_dir);
		}
		
		# Create setup directory for storing structures and setup data
		my $setup_dir = $transform_dir."/setup";
		if(! -e $setup_dir){
			mkdir($setup_dir);
		}
		
		# Structure matching and soft core determination #
		##################################################
		
		if($c_in{'match'} == 1){
			my $setup_successful = gen_match_list(\%c_in);
			if($setup_successful == 1){
				print "\n\nAtom matching list is provided in file ".$transform_dir."/match_list.txt\n\n";
			}
			else{
				print "\n\nPlease ensure correctness of matching list, before you proceed with the next setup step.\n";
				exit;
			}
		}
	
		# Generate pdb-files of solvated systems for start- and end-structures
		#######################################################################

		if($c_in{'crd_top'} ==  1){
			
			# Copy coordinated files of start structures to setup directory for consistency
			# and to avoid problems due to paths longer than 80 characters
			my $q_rst_lig = $c_in{'rst_lig'};
			my $t_rst_lig = $c_in{'Lname_start'}."_lig.restrt";
			copy($q_rst_lig, $setup_dir."/".$t_rst_lig) || die "Cannot copy file ".$q_rst_lig." to ".$setup_dir."/".$t_rst_lig."!\n";

			my $q_rst_com = $c_in{'rst_com'};
			my $t_rst_com = $c_in{'Lname_start'}."_com.restrt";
			copy($q_rst_com, $setup_dir."/".$t_rst_com) || die "Cannot copy file ".$q_rst_com." to ".$setup_dir."/".$t_rst_com."!\n";
	
			# Copy toplogy files of start structure to setup directory for consistency
			# and to avoid problems due to paths longer than 80 characters
			my $q_top_lig = $c_in{'top_lig'};
			my $t_top_lig = $c_in{'Lname_start'}."_lig.top";
			copy($q_top_lig, $setup_dir."/".$t_top_lig) || die "Cannot copy file ".$q_top_lig." to ".$setup_dir."/".$t_top_lig."!\n";
			my $q_top_com = $c_in{'top_com'};
			my $t_top_com = $c_in{'Lname_start'}."_com.top";
			copy($q_top_com, $setup_dir."/".$t_top_com) || die "Cannot copy file ".$q_top_com." to ".$setup_dir."/".$t_top_com."!\n";			

			# Change into setup directory
			chdir($setup_dir);
	
			# Convert start structures into PDB-files
			my $pdb_lig_start = "./".$c_in{'Lname_start'}."_lig.pdb";
			my $pdb_convert_lig_log = $setup_dir."/ambpdb_convert_".$c_in{'Lname_start'}."_lig.log";
			my $pdb_convert_lig = "ambpdb -p ".$t_top_lig." -c ".$t_rst_lig." > ".$pdb_lig_start;
			$pdb_convert_lig .= " 2> ".$pdb_convert_lig_log;
			print "Running antechamber for conversion of ligand start-structure.\n";
			call_prog($pdb_convert_lig, $pdb_convert_lig_log);
			unlink($pdb_convert_lig_log);

			my $pdb_com_tmp = "./".$c_in{'Lname_start'}."_com.pdb_tmp";
			my $pdb_convert_com_log = $setup_dir."/ambpdb_convert_".$c_in{'Lname_start'}."_com.log";
			my $pdb_convert_com = "ambpdb -p ".$t_top_com." -c ".$t_rst_com." > ".$pdb_com_tmp;
			$pdb_convert_com .= " 2> ".$pdb_convert_com_log;
			print "Running antechamber for conversion of complex start-structure.\n";
			call_prog($pdb_convert_com, $pdb_convert_com_log);
			unlink($pdb_convert_com_log);
			
			# Insert missing TER if different chains are present in receptor
			# Necessary because ambpdb removes TER breaks between chains in conversion of 
			# coordinate to PDB file
			my $pdb_com_start = $setup_dir."/".$c_in{'Lname_start'}."_com.pdb";
			modify_blocking_groups($pdb_com_tmp);
			insert_ter($pdb_com_tmp, $pdb_com_start, \%c_in);
			unlink($pdb_com_tmp);
	
			# Copy pdb-files of ligand and complex start-structures to end structures
			my $pdb_lig_end = $setup_dir."/".$c_in{'Lname_end'}."_lig.pdb";
			copy($pdb_lig_start, $pdb_lig_end) || die "Cannot copy $pdb_lig_start to $pdb_lig_end.\n";
	
			my $pdb_com_end = $setup_dir."/".$c_in{'Lname_end'}."_com.pdb";
			copy($pdb_com_start, $pdb_com_end) || die "Cannot copy $pdb_com_start to $pdb_com_end.\n";
	
			# Exchange residue name of ligand in start and end-structures
			change_res_name($pdb_lig_start, $setup_dir, $c_in{'Lalias_start'}, \%c_in);
			change_res_name($pdb_com_start, $setup_dir, $c_in{'Lalias_start'}, \%c_in);
			change_res_name($pdb_lig_end, $setup_dir, $c_in{'Lalias_end'}, \%c_in);
			change_res_name($pdb_com_end, $setup_dir, $c_in{'Lalias_end'}, \%c_in);
	
			# Remove superfluous atoms from PDB-files of end-structures
			if(($c_in{'mask0'})&&($c_in{'mask0'} ne "")){
				my @sc1 = split(/\@/, $c_in{'mask0'});
				del_atoms($pdb_lig_end, $setup_dir, $sc1[1], \%c_in);
				del_atoms($pdb_com_end, $setup_dir, $sc1[1], \%c_in);
			}
			
			if($c_in{'sybyl_mol2'} == 1){
				# Generate mol2-file with SYBYL atom types for easy visualization for renaming of end-structure
				my $start_mol2_for_compare = $c_in{'i_path'}."/".$c_in{'Lname_start'}.".mol2";
				my $sybyl_start_output_file = "./".$c_in{'Lname_start'}."_sybyl.mol2";
				my $chg_type_start_log = $setup_dir."/antechamber_prepareSYBYLmol2_v0.log";
				my $chg_type_start_c = "antechamber -i $start_mol2_for_compare -fi mol2 -o $sybyl_start_output_file -fo mol2";
				$chg_type_start_c .= " -at sybyl > ".$chg_type_start_log." 2>&1";
				call_prog($chg_type_start_c, $chg_type_start_log);
				unlink($chg_type_start_log);
	
				my $end_mol2_for_compare = $c_in{'i_path'}."/".$c_in{'Lname_end'}.".mol2"; 
				my $sybyl_end_output_file = "./".$c_in{'Lname_end'}."_sybyl.mol2";
				my $chg_type_end_log = $setup_dir."/antechamber_prepareSYBYLmol2_v1.log";
				my $chg_type_end_c = "antechamber -i $end_mol2_for_compare -fi mol2 -o $sybyl_end_output_file -fo mol2";
				$chg_type_end_c .= " -at sybyl > ".$chg_type_end_log." 2>&1";
				call_prog($chg_type_end_c, $chg_type_end_log);
				unlink($chg_type_end_log);
	
				unlink glob "./ANTECHAMBER*";
				unlink("./ATOMTYPE.INF");
			}
			
			# Set matching-list to default location, if no matching list is provided
			if((! $c_in{'match_list'})||($c_in{'match_list'} eq "")){
				$c_in{'match_list'} = $transform_dir."/match_list.txt";
			}
		
			# Copy mol2 files of start- end end-structures with charges to present folder 
			# to enable direct comparison of atom names
			my $leap_dir = $c_in{'root_path'}."/leap";
			my $start_mol2_org = $leap_dir."/".$c_in{'Lname_start'}."/".$c_in{'Lname_start'}."_".$c_in{'chrg_meth'}.".mol2";
			my $end_mol2_org = $leap_dir."/".$c_in{'Lname_end'}."/".$c_in{'Lname_end'}."_".$c_in{'chrg_meth'}.".mol2";
			my $start_mol2_dest = $setup_dir."/".$c_in{'Lname_start'}.".mol2";
			my $end_mol2_dest = $setup_dir."/".$c_in{'Lname_end'}.".mol2";
	
			copy($start_mol2_org, $start_mol2_dest) || die "Cannot copy $start_mol2_org to $start_mol2_dest.\n";
			copy($end_mol2_org, $end_mol2_dest) || die "Cannot copy $end_mol2_org to $end_mol2_dest.\n";
	
			# Shift position of soft core part to end of mol2-file
			shift_position($start_mol2_dest, "mask0", \%c_in);
			shift_position($end_mol2_dest, "mask1", \%c_in);
	
			# Generation of parameter and library files for end-structure based on mol2-file
			# with charges and atom names corresponding to atom names of start-structure.
					
			# Generate AMBER-conform mol2-file for state V1 including modification of unit name
			my $ante_com1_log = $setup_dir."/antechamber_v1_mol2.log";
			my $ante_com1 = "antechamber -i ".$end_mol2_dest. " -fi mol2";
			$ante_com1 .= " -o ".$end_mol2_dest."_amb -fo mol2";
			$ante_com1 .= " > ".$ante_com1_log." 2>&1";
			print "Generating new mol2-file\n";
			call_prog($ante_com1, $ante_com1_log);
			unlink glob "./ANTECHAMBER*";
			unlink("./ATOMTYPE.INF");
			unlink($ante_com1_log);
		
			open(MOL, $end_mol2_dest."_amb") || die "Cannot open file ".$end_mol2_dest."_amb\n";
			open(TMP, ">tmp.mol2") || die "Cannot open file tmp.mol2 for writing\n";
		
			while(my $mol_l = <MOL>){
				chomp($mol_l);
				$mol_l =~ s/UNN/$c_in{'Lalias_end'}/;
				$mol_l =~ s/<1>/$c_in{'Lalias_end'}/;
				print TMP $mol_l."\n";
			}
			close MOL;
			close TMP;
		
			unlink($end_mol2_dest."_amb");
			copy("tmp.mol2", $end_mol2_dest."_amb");
			unlink("tmp.mol2");
		
			# Generate parameter file for end-structure
			my $frc_log = $setup_dir."/parmchk_".$c_in{'Lname_end'}.".log";
			my $frc_command = "parmchk2 -i ".$end_mol2_dest."_amb -f mol2 -o ./".$c_in{'Lname_end'}.".frcmod";
			$frc_command .= " 2> ".$frc_log; 		
			print "Generating parameters for ".$c_in{'Lname_end'}."\n";
			call_prog($frc_command, $frc_log);
			unlink($frc_log);
		
			# Create library file
			foreach my $k ($c_in{'Lname_start'}, $c_in{'Lname_end'}){
				my $lib_leap_in = $setup_dir."/leap_lib_".$k.".in";
				my $lib_file = $setup_dir."/".$k.".lib";
				my $frcmod_file;
				my $alias;
				if($k eq $c_in{'Lname_start'}){
					$alias = $c_in{'Lalias_start'};
					$frcmod_file = $c_in{'root_path'}."/leap/".$c_in{'Lname_start'}."/".$c_in{'Lname_start'}.".frcmod";
				}
				else{
					$alias = $c_in{'Lalias_end'};
					$frcmod_file = $setup_dir."/".$k.".frcmod";
				}
	
				setup_leap_for_lib_TI($lib_leap_in, $lib_file, $frcmod_file, $alias, $setup_dir, $end_mol2_dest, \%c_in);
	
				# Run leap
				my $lib_leap_log = $setup_dir."/leap_lib_".$k.".log";
				my $lib_leap_c = "tleap -f ".$lib_leap_in." >> ".$lib_leap_log." 2 >> ".$lib_leap_log;
			
				print "Running LEaP for lib-file generation for $k...\n";
				call_prog($lib_leap_c, $lib_leap_log);
			
				# Check correct performance
				check_warn_chrg_leap($lib_leap_log, $lib_file);
			
				unlink("./leap.log");
			}
		
			# Change residue name in library file of the start structure
			my $start_lib = $setup_dir."/".$c_in{'Lname_start'}.".lib";
			change_res_name($start_lib, $setup_dir, $c_in{'Lalias_start'}, \%c_in);
	
			# Create coordinate and topology files
			if((! -e "$setup_dir/$c_in{'Lname_start'}"."_lig_TIin.crd")&&(! -e "$setup_dir/$c_in{'Lname_end'}"."_lig_TIin.crd")
			&&(! -e "$setup_dir/$c_in{'Lname_start'}"."_com_TIin.crd")&&(! -e "$setup_dir/$c_in{'Lname_end'}"."_com_TIin.crd")){
				my $ref_cys_bridges;
				if($c_in{'cys_f'} ne ""){
					$ref_cys_bridges = read_disulf(\%c_in);
				}
	
				gen_leap_crd_top_in($setup_dir, \%c_in, $ref_cys_bridges, "sander");
				print "\nGenerating coordinates and topologies for ".$c_in{'Lname_start'}." and ".$c_in{'Lname_end'}."\n";
	
				my $leap_ct_in = "$setup_dir/leap_top_crd.in";
				my $leap_ct_log = "$setup_dir/leap_top_crd.log";
				my $leap_ct =  "tleap -f ".$leap_ct_in." > ".$leap_ct_log." 2>&1";
				call_prog($leap_ct, $leap_ct_log);
	
				&check_leap($leap_ct_log);
				setBox_info($setup_dir, \%c_in, "sander");
				unlink("./leap.log");
			}

	
			# If setup with pmemd was requested, generate new PDB files with ligands at the top 
			# of the file and create coordinate and topology files out of these PDB files for pmemd
			if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'}==1)){

				my $ref_cys_bridges;
                                if($c_in{'cys_f'} ne ""){
                                        $ref_cys_bridges = read_disulf(\%c_in);
                                }

				# Generate new PDBs for setup with both ligands at the top of the PDB file
				# followed by the complex, ions, and water
				gen_PDBs_pmemd($setup_dir, \%c_in);
				my $leap_pdb_in = "$setup_dir/leap_pdb_pmemd.in";
				my $leap_pdb_log = "$setup_dir/leap_pdb_pmemd.log";
				open(PDB_PMEMD, ">".$leap_pdb_in);

				print PDB_PMEMD "source leaprc.protein.ff14SB\n";
				print PDB_PMEMD "source leaprc.water.tip3p\n";
				print PDB_PMEMD "source leaprc.gaff\n";
				print PDB_PMEMD "loadoff ".$setup_dir."/".$c_in{'Lname_start'}.".lib\n";
				print PDB_PMEMD "loadoff ".$setup_dir."/".$c_in{'Lname_end'}.".lib\n";

				for my $tag ("_com", "_lig"){
					print PDB_PMEMD "PDB".$tag." = loadpdb ".$setup_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}.$tag."_TIin.pdb\n";
					print PDB_PMEMD "savepdb PDB".$tag." ".$setup_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}.$tag."_TIin_leap.pdb\n";
				}
				print PDB_PMEMD "quit\n";	
				my $leap_pdb = "tleap -f ".$leap_pdb_in." > ".$leap_pdb_log." 2>&1";
				call_prog($leap_pdb, $leap_pdb_log);

				# Generate coordinate and topology files
				gen_leap_crd_top_in($setup_dir, \%c_in, $ref_cys_bridges, "pmemd");

				print "\nGenerating coordinates and topologies for ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."\n";
				my $leap_ct_in = "$setup_dir/leap_top_crd_pmemd.in";
				my $leap_ct_log = "$setup_dir/leap_top_crd_pmemd.log";
				my $leap_ct =  "tleap -f ".$leap_ct_in." > ".$leap_ct_log." 2>&1";
				call_prog($leap_ct, $leap_ct_log);

				&check_leap($leap_ct_log);

				# Modify coordinate and topology files with ParmED
				generate_parmed_input($setup_dir, \%c_in);

				print "\nGenerating TI input coordingates and topologies for ".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."\n";
				my $parmed_in = "$setup_dir/parmed.in";
				my $parmed_log = "$setup_dir/parmed.log";
				my $parmed_c = "parmed -i ".$parmed_in." > ".$parmed_log." 2>&1";

				call_prog($parmed_c, $parmed_log);

				# Set new box dimennsions
				setBox_info($setup_dir, \%c_in, "pmemd");
				unlink("./leap.log");
			} 
	
		} # End: Setup of coordinate and topology files
	
	
		######################
		# Prepare Simulation #
		######################
	
		my $com_res_no; # Total number of residues in complex (including water)
		my $lig_res_no; # Total number of residues in ligand (including water)
		my $com_res_lig; # Residue number of ligand in complex

		# Setup scripts for MD equilibration
		######################################

		if($c_in{'ti_equi'} == 1){
	
			# Information for used
			print "\nPreparing files for TI equilibration...\n";

			# Determine residue number
			($com_res_no, $lig_res_no, $com_res_lig) = determine_no_of_res($c_in{'Lname_start'}, $c_in{'Lalias_start'}, $setup_dir);
	
			# Determine new softcore of V1 after automatic renaming by leap
			$c_in{'mask1'} = determine_new_sc_str($setup_dir, \%c_in);
	
			if(($c_in{'equi_templ'})&&(-e $c_in{'equi_templ'})){
	
				# Create equilibration directory
 				my $equi_dir = $transform_dir."/equi";
				if(! -e $equi_dir){
					mkdir($equi_dir);
				}
				chdir($equi_dir);
		
				my @tags = ("lig", "com");

				foreach my $tag (@tags){
				
					# Create equilibration "tag" directory
					my $equi_tag_dir = $equi_dir."/".$tag;
					if(! -e $equi_tag_dir){
						mkdir($equi_tag_dir);
					}

					# Copy topology and coordinate files to TI equilibration directory
					if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
						for my $file_type("top", "crd"){
							my $q_file = $setup_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.".$file_type;
							my $t_file = $equi_tag_dir."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."_".$tag."_TIin.".$file_type;
							copy($q_file, $t_file) || die "Cannot copy file $q_file to $t_file!\n";
						}
					}
					if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
						for my $name ($c_in{'Lalias_start'}, $c_in{'Lalias_end'}){
							for my $file_type ("top", "crd"){
								my $q_file = "";
						
								if($name eq $c_in{'Lalias_start'}){
									$q_file = $setup_dir."/".$c_in{'Lname_start'}."_".$tag."_TIin.".$file_type;
								}
								else{
									$q_file = $setup_dir."/".$c_in{'Lname_end'}."_".$tag."_TIin.".$file_type;
								}	
						
								my $t_file = $equi_tag_dir."/".$name."_".$tag."_TIin.".$file_type;
								copy($q_file, $t_file) || die "Cannot copy file $q_file to $t_file!\n";
							}
						}
					}
			
					# Setup equilibration
					my $prod_no = equi_setup_TI($equi_tag_dir, $tag, $lig_res_no, $com_res_no, \%c_in);
				
					my $pbs_work_dir = "";
					if((exists $c_in{'pbs_t_path'})&&($c_in{'pbs_t_path'} ne "")){
						$pbs_work_dir = $c_in{'pbs_t_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."/equi/".$tag;
					}
					else{
						$pbs_work_dir = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."/equi/".$tag;
					}
					generate_pbs_equi_TI($equi_tag_dir, $pbs_work_dir, $tag, $prod_no, \%c_in);
				}
			}
		} # End: Prepare equilibration
	
	
		# Setup directory structure and scripts for MD_production
		###########################################################

		if($c_in{'ti_prod'} == 1){
			if(($c_in{'prod_templ'})&&(-e $c_in{'prod_templ'})){
	
				# Information for user
				print "\nPreparing files for TI production...\n";
	
				# Determine new softcore of V1 after automatic renaming by leap
				if($c_in{'ti_equi'} == 0){
					$c_in{'mask1'} = determine_new_sc_str($setup_dir, \%c_in);
				}
	
				# Determine number of lambda steps
				my @lambda = split(/\,/, $c_in{'lambda'});
				my @prod_lambda = split(/\,/, $c_in{'prod_lambda'});
		
				# Assign production time to lambda
				my @prod_time = split(/\,/, $c_in{'prod_time'});	
				my %prod_time;
				for my $l (0..@prod_lambda){
					$prod_time{$prod_lambda[$l]} = $prod_time[$l];
				}
	
				# Determine residue number
				($com_res_no, $lig_res_no, $com_res_lig) = determine_no_of_res($c_in{'Lname_start'}, $c_in{'Lalias_start'}, $setup_dir);
	
				# Set equilibration directory
				my $equi_tag_dir = $transform_dir."/equi";
			
				if(! -e $equi_tag_dir){
					print "\nERROR:\nTI-equilibration must be performed before production simulations are setup.\n";
					exit;
				}
				else{
				
					# Check completion of equilibration phase
					my $equi_check = check_TIequi($transform_dir, \%c_in);
				
					if($equi_check != 1){
						print "\n\n ATTENTION:\n";
						print "The equilibration check was not successful. Nevertheless input files for the production\n";
						print "calculations were prepared. These can be used to run longer equilibration simulations\n";
						print "at individual lambda steps. It is strongly recommended to carefully inspect the equilibration\n";
						print "simulations to detect potential problems. In case the equilibration runs are extended with\n";
						print "the prepared production input files, the recorded dV/dL values should be plotted against the\n";
						print "simulation time to find out from what point onward the systems can be regarded as equilibrated.\n";
					}
				
					# Create production directory
					my $prod_dir = $transform_dir."/prod";
					if(! -e $prod_dir){
						mkdir($prod_dir);
					}
					chdir($prod_dir);
		
					# Determine mask for re-imaging of atoms common in V0 and V1
					# using cpptraj during production run
					my $mask = identify_mask($setup_dir, \%c_in);
					my @tags = ("com", "lig");

					# Prepare qsub-file for submission of batch jobs for TI production
					open(QSUB_TI_PROD, ">".$prod_dir."/qsub_TI_prod.sh");

					foreach my $tag (@tags){
			
						# Create "Tag" directory
						my $MD_tag_dir = $prod_dir."/".$tag;
						if(! -e $MD_tag_dir){
							mkdir($MD_tag_dir);
						}
		
						# Determine basename of coordinate files from last step of equilibration
						my $equi_path_and_base = "";
						if((! exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 0)){
							$equi_path_and_base = $equi_tag_dir."/".$tag."/".$c_in{'Lalias_end'};
						}
						if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
							$equi_path_and_base = $equi_tag_dir."/".$tag."/md_";
						}
						my $last_equi_basename = detect_last_equi_basename($equi_path_and_base, \@prod_lambda, \%c_in);
		
						# Copy topology and coordinate files to TI MD-directory
						my @prefix;
						if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
							@prefix = ($c_in{'Lalias_start'}, $c_in{'Lalias_end'});
						}
						if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
							@prefix = ($c_in{'Lname_start'}."_".$c_in{'Lname_end'});
						}

						for my $pre (@prefix){
					
							#Copy topology
							my $q_topo = $equi_tag_dir."/".$tag."/".$pre."_".$tag."_TIin.top";
							my $t_topo = $MD_tag_dir."/".$pre."_".$tag."_TIin.top";
						
							copy($q_topo, $t_topo) || die "Cannot copy $q_topo to $t_topo.\n";
							
							# Copy coordinate files from last step of equilibration									
							foreach my $l (@prod_lambda){
								my $q_crd = "";
								if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'}) == 0){
									$q_crd = $equi_tag_dir."/".$tag."/".$pre."_".$last_equi_basename."_";
								}
								if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'}) == 1){
									$q_crd = $equi_tag_dir."/".$tag."/md_".$last_equi_basename."_";
								}

								my $t_crd = "$MD_tag_dir/equi_";
								if($pre eq $c_in{'Lalias_start'}){
									$q_crd = $q_crd."v0_l".$l.".rst";
									$t_crd = $t_crd."v0_l".$l.".rst";
								}
								elsif($pre eq $c_in{'Lalias_end'}){
									$q_crd = $q_crd."v1_l".$l.".rst";
									$t_crd = $t_crd."v1_l".$l.".rst";
								}
								else{
									$q_crd = $q_crd."l".$l.".rst";
									$t_crd = $t_crd."l".$l.".rst";
								}
							
								# Copy all except time statement in first line
								my $first_line = 0;
								my $crd_str_out = "";
								
								open(QCRD, $q_crd) || die "Cannot open file $q_crd for reading.\n";
								$crd_str_out = <QCRD>; # Comment line
								while(my $qcrd_l = <QCRD>){
									my $str_chk = substr($qcrd_l, 0, 6);
									if($first_line == 0){
										$qcrd_l =~ s/^\s+//;
										my @qcrd_l = split(/\s+/, $qcrd_l);
										$qcrd_l = sprintf("%5s", $qcrd_l[0]);
										$qcrd_l .= "  ";
										$qcrd_l .= sprintf("%.7e\n", 0);
										$first_line = 1;
									}
									$crd_str_out .= $qcrd_l;
								}
								close QCRD;
							
								open(TCRD, ">$t_crd") || die "Cannot open file $t_crd for writing.\n";	
								print TCRD $crd_str_out;
								close TCRD;
							}
						}
						
					
						# Create MD-setup					
						generate_cpptraj_img($MD_tag_dir, $tag, $com_res_lig, $mask, \%c_in);
					
						my $requested_prod_time_templ = determine_requested_time("production", $c_in{'prod_templ'}, \%c_in);
						foreach my $lambda (@prod_lambda){
							my $prod_no = $prod_time{$lambda}/$requested_prod_time_templ;
							my $restrt = "";

							if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
								group_file_prod_setup($MD_tag_dir, $restrt, $tag, "1", $lambda, \%c_in);
							}

							prod_setup_TI($MD_tag_dir, $tag, $lig_res_no, $com_res_no, "0", "1", $lambda, \%c_in);
							
							my $pbs_work_dir;
							if((exists $c_in{'pbs_t_path'})&&($c_in{'pbs_t_path'} ne "")){
								$pbs_work_dir = $c_in{'pbs_t_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."/prod/".$tag;
							}
							else{
								$pbs_work_dir = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."/prod/".$tag;
							}
								
							if(($c_in{'pbs_prod_t'} ne "")&&(-e $c_in{'pbs_prod_t'})){
								generate_pbs_prod_TI($MD_tag_dir, $pbs_work_dir, $current_dir, $tag, $lambda, $prod_no, $requested_prod_time_templ, \%c_in);	 
								print QSUB_TI_PROD "qsub ".$pbs_work_dir."/run_prod.pbs.l".$lambda."\n";
							}
						}
					
					} # End: tags
					close QSUB_TI_PROD;
					chmod 0755, "$prod_dir/qsub_TI_prod.sh";
				} # End: Ensure presence of equilibration data		
			} # End: Ensure existence of production template
		} # End: ti_prod flag
	} # End: TI simulation setup
	
	
	###############################################
	# Calculate difference in binding free energy #
	###############################################
	
	TI_DDG:
	
	if(($c_in{'ti_ddG'})&&($c_in{'ti_ddG'} == 1)){
	
		my $check_dG_params = check_dG_params(\%c_in);
		if($check_dG_params != 1){
			exit;
		}
		
		print "\nCalculationg difference in binding free energy between ".$c_in{'Lname_start'}." and ".$c_in{'Lname_end'}.".\n";
	
		# Check existence of required folders / files
		my $transform_dir = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'};
		my $prod_dir = $transform_dir."/prod";

		# Create analysis directory
		my $calc_dir = $transform_dir."/TI_results";
		if(! -e $calc_dir){
			mkdir($calc_dir);
		}
		
		# Create log-file
		open(LOG, ">".$calc_dir."/TI_calculation.log") || die "Cannot open log-file for writing.\n";
		print LOG "# Log-file for calculation of ddG energy with TIW module of FEW\n";
		print LOG "#######################################################################\n";
		
		# Create result-file
		open(OUT, ">".$calc_dir."/TI_dG.out") || die "Cannot open result-file for writing.\n";
		print OUT "# Difference in binding free energy determined by thermodynamic integration #\n";
		print OUT "#############################################################################\n";
		
		# Identify dV/dL source
		my @dVdL_src = ();
		if(defined $c_in{'dVdL_src'}){
			if($c_in{'dVdL_src'} == 0){
				push(@dVdL_src, "1");
				push(@dVdL_src, "all");
			}
			elsif($c_in{'dVdL_src'} =~ m/(\d+)\-0/){
				push(@dVdL_src, $1);
				push(@dVdL_src, "all");
			}
			elsif($c_in{'dVdL_src'} =~ m/(\d+)\-(\d+)/){
				push(@dVdL_src, $1);
				push(@dVdL_src, $2);
			}
			else{
				print "\nERROR: Wrong pattern used for definition of the files that shall be\n";
				print "considered in the calculation of ddG.\n";
				exit;
			}
		}
		
		# Calculate ddG
		my $ddG;
		my $ddG_err;
		my %dG;
		my %dG_err;
		foreach my $tag ("com", "lig"){
		
			my $sys;
			if($tag eq "com"){
				$sys = "Complex";
			}
			else{
				$sys = "Ligand";
			}
			print OUT "\n".$sys.":\n";
		
			# Detect available files and sort them by lambda
			my @sander_out_files = <$prod_dir/$tag/*.out>;
			my %files_per_lambda;
			my $end_structure = $c_in{'Lalias_end'};
			my %no_files_per_lambda;
			my %lambdas;
			my $init = 0;
			my $pmemd_ti = 0;
			
			foreach my $file (@sander_out_files){
				my @filename_components = split(/\//, $file);

				if($filename_components[$#filename_components] =~ m/md\_prod(\d+)\_l(\d+)/){
					$pmemd_ti = 1;
				}	

				if(($filename_components[$#filename_components] =~ m/$end_structure\_prod(\d+)\_v1\_l(\d+)/)
				or($filename_components[$#filename_components] =~ m/md\_prod(\d+)\_l(\d+)/)){
					my $prod_no = $1;
					my $lambda = $2;
					$files_per_lambda{$lambda}->{$prod_no} = $file;
					
					if($init == 0){
						$no_files_per_lambda{$lambda} = $prod_no;
						$lambdas{$lambda} = $lambda;
					}
					else{
						if(! defined $lambdas{$lambda}){
							$lambdas{$lambda} = $lambda;
						}
						if($no_files_per_lambda{$lambda} < $prod_no){
							$no_files_per_lambda{$lambda} = $prod_no;
						}
					}
				}
			}
			
			my @sorted_lambdas = sort {$lambdas{$a} <=> $lambdas{$b}} keys %lambdas;
		
			# Values for lamba need to be always larger 1
			if(@sorted_lambdas == 0){
				print "\nERROR: No Lambda values were found. Please ensure that the file location and format\n";
				print "is consistent with the one used by FEW.\n";
				exit;
			}
			if(@sorted_lambdas == 1){
				print "\nWARNING: Only one lambda value was found for analysis. Calculation will be conducted\n";
				print "by conventional numerical integration.\n";
				if($c_in{'calc_meth'} != 1){
					$c_in{'calc_meth'} = 1;
				}
			}
			
			# Determine files that shall be regarded
			foreach my $lambda (@sorted_lambdas){
				for(my $prod_no=1; $prod_no<=$no_files_per_lambda{$lambda}; $prod_no++){
					if($prod_no < 10){
						$prod_no = "0".$prod_no;
					}
					if($prod_no < $dVdL_src[0]){
						delete($files_per_lambda{$lambda}->{$prod_no});
					}
					if(($dVdL_src[1] ne "all")&&($prod_no > $dVdL_src[1])){
						delete($files_per_lambda{$lambda}->{$prod_no});
					}
				}
			}
			
			# Print information about analysis to log-file
			print LOG "\n";
			print LOG "$sys files considered in calculation:\n";
			foreach my $lambda (@sorted_lambdas){
				print LOG "\nLambda: 0.".$lambda."\n";
				my @sorted_prods = sort {$files_per_lambda{$lambda}->{$a} cmp $files_per_lambda{$lambda}->{$b}} keys %{$files_per_lambda{$lambda}};
				foreach my $prod_no (@sorted_prods){
					if(defined $files_per_lambda{$lambda}->{$prod_no}){
						print LOG $files_per_lambda{$lambda}->{$prod_no}."\n";
					}
				}
			}
			
			# Read dV/dL values and calculate average and standard deviation
			my $sum_dVdL = 0;
			my @dVdL;
			my %avg_dVdL;
			my %std_dVdL;
			my %std_err;
			my $value_no = 0;
			my $read_dVdL = 0;
			my $read_dVdL_two = 0;
			my $curr_time = 0;
			my $time_inter = 0;

			foreach my $lambda (@sorted_lambdas){
				my @sorted_prods = sort {$files_per_lambda{$lambda}->{$a} cmp $files_per_lambda{$lambda}->{$b}} keys %{$files_per_lambda{$lambda}};
				my $dVdL_out = $calc_dir."/dVdL_".$tag."_l".$lambda;
				my $save_time = 0;
				$time_inter = 0;
				my $time_read = 0;
				
				open(DVDL_OUT, ">".$dVdL_out) || die "Cannot open file $dVdL_out for writing.\n";
				print DVDL_OUT "NSTEP\tdVdL\n";
				
				foreach my $prod_no (@sorted_prods){
					open(DVDL_SRC, $files_per_lambda{$lambda}->{$prod_no}) || die "Cannot open file ".$files_per_lambda{$lambda}->{$prod_no}." for reading.\n";

					while(my $dVdL_line = <DVDL_SRC>){
						chomp($dVdL_line);
			
						if($dVdL_line =~ m/4\.  RESULTS/){
							$read_dVdL = 1;
						}

						if($dVdL_line =~ m/TI region  1/){
							$read_dVdL_two = 1;
						}

						if($dVdL_line =~ m/TI region  2/){
							$read_dVdL_two = 0;
						}

						if($dVdL_line =~ m/A V E R A G E S/){
							$read_dVdL = 0;
							close DVDL_SRC;
							last;
						}
			
						if((($pmemd_ti == 0)&&($read_dVdL == 1)&&($dVdL_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))
						||(($pmemd_ti == 1)&&($read_dVdL_two == 1)&&($dVdL_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))){
							$curr_time = $2.".".$3;
						}
						
						if((($pmemd_ti == 0)&&($read_dVdL == 1)&&($time_read == 0)&&($dVdL_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))
						||(($pmemd_ti == 1)&&($read_dVdL_two == 1)&&($time_read == 0)&&($dVdL_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/))){
							my $time = $2.".".$3;
			
							if($save_time != 0){
								$time_inter = $time - $save_time;
								$time_read = 1;
							}
							$save_time = $time;
						}
						
						if((($pmemd_ti == 0)&&($read_dVdL == 1)&&($dVdL_line =~ m/^ DV\/DL/))
						||(($pmemd_ti == 1)&&($read_dVdL_two == 1)&&($dVdL_line =~ m/^ DV\/DL/))){
							my @value = split(/\s+/, $dVdL_line);
							$sum_dVdL = $sum_dVdL + $value[3];
							$dVdL[$value_no] = $value[3];
							print DVDL_OUT $curr_time."\t".$value[3]."\n";
							$value_no++;
						}
					}
				}
				close DVDL_OUT;
				
				$avg_dVdL{$lambda} = $sum_dVdL / $value_no;
				$avg_dVdL{$lambda} = sprintf("%.4f", $avg_dVdL{$lambda});
				print LOG "Mean dVdL for Lambda 0.".$lambda.": ".$avg_dVdL{$lambda}." ";
				
				# Determine standard deviation
				my $std_sum = 0;
				$std_dVdL{$lambda} = 0;
				$std_err{$lambda} = 0;

				for(my $v=0; $v<$value_no; $v++){
					$std_sum = ($std_sum + (($dVdL[$v] - $avg_dVdL{$lambda})**2));
				}
				my $factor = (1 / ($value_no - 1));
				$std_dVdL{$lambda} = $factor * $std_sum;
				$std_dVdL{$lambda} = sqrt($std_dVdL{$lambda});
				
				print LOG "+/- ".$std_dVdL{$lambda}."\t";
				
				# Determine total simulation time
				my $total_time = $time_inter * $value_no;
				
				# Read autocorrelation time
				my $autocorr = 0;
				my $autocorr_file = $prod_dir."/".$tag."/stderr_".$lambda;
				open(AUTOCORR, $autocorr_file)|| die "Cannot open file $autocorr_file.\n";
				while(my $line = <AUTOCORR>){
					if($line =~ m/Autocorrelation: (\d+)/){
						$autocorr = $1;
					}
				}
				
				# Determine standard error
				my $factor_autocorr = $total_time / (2*$autocorr*$time_inter);
				$std_err{$lambda} = 1/sqrt($factor_autocorr);
				$std_err{$lambda} = $std_err{$lambda} * $std_dVdL{$lambda};
				$std_err{$lambda} = sprintf("%.4f", $std_err{$lambda});
				
				print LOG "STDERR: ".$std_err{$lambda}."\n";
				
				$value_no = 0;
				$sum_dVdL = 0;
			}
			
			# Determine real lambda_values
			my @lambda_values;
			foreach my $l (@sorted_lambdas){
				my $l_value = '0.'.$l;
				push(@lambda_values, $l_value);
			}
			
			# Compute dG
			my %weights;
			if(($c_in{'calc_meth'})&&($c_in{'calc_meth'} == 1)){
			
				# Interpolate to lambda=0.0
				if($sorted_lambdas[0] != 0) {
					$avg_dVdL{'0'} = ($lambda_values[0]*$avg_dVdL{$sorted_lambdas[1]} 
									- $lambda_values[1]*$avg_dVdL{$sorted_lambdas[0]}) 
									/ ($lambda_values[0]-$lambda_values[1]);
					$avg_dVdL{'0'} = sprintf("%.4f", $avg_dVdL{'0'});
					$std_err{'0'} = $std_err{$sorted_lambdas[0]};
					unshift(@sorted_lambdas, "0");
					unshift(@lambda_values, "0.0");
				}
				
				# Interpolate to lambda=1.0
				# Lambda of 1 cannot be handled and will therefore never occur -> Values of lambda from 1 to 99
				# are accepted -> Use "100" as alias for lamba=1.0
				my $n = @sorted_lambdas;
				$avg_dVdL{'100'} = (($lambda_values[$n-2] - 1)*$avg_dVdL{$sorted_lambdas[$n-1]}
								  + (1 - $lambda_values[$n-1])*$avg_dVdL{$sorted_lambdas[$n-2]}) 
								  / ($lambda_values[$n-2] - $lambda_values[$n-1]);
				$avg_dVdL{'100'} = sprintf("%.4f", $avg_dVdL{'100'});
				$std_err{'100'} = $std_err{$sorted_lambdas[$n-1]};
				push(@sorted_lambdas, 100);
				push(@lambda_values, "1.0");
				
				# Determine interval weights
				for(my $n=0; $n<@sorted_lambdas; $n++) {
					if(($n == 0)&&($sorted_lambdas[$n] == 0)){
						$weights{$sorted_lambdas[$n]} = 0.5 * ( $lambda_values[$n] + $lambda_values[$n+1] );
					}
					elsif(($n == (@sorted_lambdas-1))&&($sorted_lambdas[$n] == 100)){
						$weights{$sorted_lambdas[$n]} = 1 - 0.5 * ( $lambda_values[$n] + $lambda_values[$n-1] );
					}
					else {
						$weights{$sorted_lambdas[$n]} = ( $lambda_values[$n] - $lambda_values[$n-1] );
					}
				}
			}
			
			if(($c_in{'calc_meth'})&&($c_in{'calc_meth'} == 2)){
				
				# Determine interval weights
				# It was ensured before that at least two values are present.
				for(my $n=0; $n<$#sorted_lambdas; $n++) {
					$weights{$sorted_lambdas[$n]} = ( $lambda_values[$n+1] - $lambda_values[$n] );
				}
				$weights{$sorted_lambdas[$#sorted_lambdas]} = $lambda_values[$#sorted_lambdas] - $lambda_values[$#sorted_lambdas-1];
			}
			
			# Determine dG and error
			$dG{$tag} = 0;
			$dG_err{$tag} = 0;
			for(my $n=0; $n<@sorted_lambdas; $n++) {
				$dG{$tag} = $dG{$tag} + ($weights{$sorted_lambdas[$n]} * $avg_dVdL{$sorted_lambdas[$n]});
				$dG_err{$tag} = $dG_err{$tag} + (($weights{$sorted_lambdas[$n]}**2) * ($std_err{$sorted_lambdas[$n]}**2));
				print OUT "Lambda: ".$lambda_values[$n]."; Weight: ".$weights{$sorted_lambdas[$n]}.
				"; dVdL: ".$avg_dVdL{$sorted_lambdas[$n]}." +/- ".$std_err{$sorted_lambdas[$n]}." kcal/mol\n";
			}
			
			$dG_err{$tag} = sqrt($dG_err{$tag});	
			print OUT "________________________________________________________________________________\n";
			$dG{$tag} = sprintf("%.4f", $dG{$tag});
			$dG_err{$tag} = sprintf("%.4f", $dG_err{$tag});
			print OUT "dG-".$sys.": ".$dG{$tag}." +/- ".$dG_err{$tag}." kcal/mol\n";
		}
		$ddG = $dG{'com'} - $dG{'lig'};
		$ddG = sprintf("%.4f", $ddG);
		$ddG_err = ($dG_err{'com'}**2) + ($dG_err{'lig'}**2);
		$ddG_err = sqrt($ddG_err);
		$ddG_err = sprintf("%.4f", $ddG_err);
		print OUT "\n================================================================================\n";
		print OUT "ddG-binding: ".$ddG." +/- ".$ddG_err." kcal/mol\n\n";
	}
} # End subroutine

1;
