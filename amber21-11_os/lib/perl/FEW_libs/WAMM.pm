#/usr/bin/perl -w

# Perl module for the preparation of structures for MD input and for
# analysis of snapshots from the simulations with the MM-PBSA and MM-GBSA methods.

use File::Copy;
use strict;

use global;
use check_input;
use separate_structures;
use common_prepare_MD;
use setup_PB_GB;
use read_data;



sub WAMM{

	my $current_dir = shift;
	my $n_traj_inter = shift;
	my $traj_files_ref = shift;
	my $traj_inter_ref = shift;
	my $c_in_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my %traj_inter = %{$traj_inter_ref};
	my %c_in = %{$c_in_ref};
	
	
	##################################
	# Specific check for WAMM module #
	##################################
	
	if(($c_in{'analyze'} == 1)||($c_in{'sim'} == 1)){
		my $check_wamm = &check_input_wamm($n_traj_inter, \@traj_files, \%traj_inter, \%c_in);
		if($check_wamm != 1){
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
	#START LW_03_2018
	my $prot_pdb;
	#END LW_03_2018

	# Case 1: Structures automatically deposited in "structs" folder of root_path
	if($c_in{'sep'} == 1){
		$struct_path = $c_in{'root_path'} . "/structs";
	}
	#START LW_03_2018
	# Case 2: Protein-Protein complex pdb as only structure
	if($c_in{'prot_prot'} == 1){
		$prot_pdb = (split '/', $c_in{'prot_com'} )[-1];
		$struct_path = $c_in{'prot_com'};
		$struct_path =~ s/\/$prot_pdb//g;
	}
	# Case 3: Structures present in input-folder with path "i_path"
	#END LW_03_2018
	else{
		$struct_path = $c_in{'i_path'};
	}

	# Get names of all files present in struct-directory
	my @struct_files = <$struct_path/*>;    # Complete path to structure files
	my %struct_names;   # Names of structure files
	my %struct_b;       # Basenames of structure files
	my $mmpbsa_pl;		# Path to mm_pbsa.pl script, ensures correct handling if defined for AMBERHOME variable
	#START LW_03_2018
	my $substr_del;

	if($c_in{'prot_prot'} != 1){
		foreach my $file (@struct_files){
			my @struct_path_split = split(/\//, $file);
			$struct_names{$file} = $struct_path_split[$#struct_path_split];
			$substr_del = 5;
			my $struct_b = substr($struct_path_split[$#struct_path_split], 0, -$substr_del);
			$struct_b{$file} = $struct_b;
		}
	}
	# For prot_prot = 1, only use one structure
	else{
		my @struct_path_split = split(/\//, $c_in{'prot_com'});
		# For protein ligands suffix is .pdb
		$substr_del = 4;
		my $struct_b = substr($struct_path_split[$#struct_path_split], 0, -$substr_del);
		$struct_b{$c_in{'prot_com'}} = $struct_b;
		# Only one structure is allowd to be analized with prot_prot = 1
		@struct_files = ();
		push (@struct_files, $struct_path."/".$prot_pdb);
	}
	#END LW_03_2018
	
	# Set MD_path to root_path if no specific MD_path is provided
	if($c_in{'MD_path'} eq ""){
		$c_in{'MD_path'} = $c_in{'root_path'};
	}

	my $amberhome =  $ENV{"AMBERHOME"};
	$mmpbsa_pl = $c_in{'mmpbsa_pl'};
	
	if(($c_in{'mmpbsa_pl'} =~ m/AMBERHOME/)&&($amberhome ne "")){
		$mmpbsa_pl =~ s/\$AMBERHOME//;
		$mmpbsa_pl = $amberhome.$mmpbsa_pl;
	}
	if(((! exists $c_in{'mmpbsa_pl'})||($c_in{'mmpbsa_pl'} eq ""))&&($amberhome ne "")){
		$mmpbsa_pl = $amberhome."/bin/mm_pbsa.pl";
		$c_in{'mmpbsa_pl'} = $amberhome."/bin/mm_pbsa.pl";
	}
	
	if(! -e $mmpbsa_pl){
		print "The path to the mm_pbsa.pl file you specified is not correct.\n";
		exit;
	}


	###################
	# Check mol2 file #
	###################
	
	if($c_in{'prot_prot'} != 1){ # LW_03_2018
		foreach my $s (@struct_files){
			# Check atom and residue naming in provided mol2 file
			my $check_mol = check_mol2($s, "WAMM", \%c_in);
			# End program if check failed
			if($check_mol != 1){
				exit;
			}
		}
	} # LW_03_2018

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

			# 3-Trajectory approach - receptor
			if($c_in{'traj'} == 3){
				mkdir("$md_dir/rec");
			}

			############################################
			#  Generate topology and coordinate files  #
			############################################

			create_crd_top("wamm", $chrg_meth, $md_dir, $current_dir, \@struct_files, \%struct_b, \%c_in);

			####################################
			#  Create scripts for equilbration #
			####################################

			# Setup folders, input scripts and pbs script for equilibration
			if(($c_in{'equi'})&&(-e $c_in{'equi'})){
				$warnings = prepare_equilibration("wamm", $chrg_meth, $md_dir, \@struct_files, \%struct_b, \%c_in);
			}


			#####################################
			#  Create scripts for MD-production #
			#####################################

			if(($c_in{'prod'})&&(-e $c_in{'prod'})&&($c_in{'prod_t'})){
				prepare_MD_production("wamm", $chrg_meth, $md_dir, $warnings, \@struct_files, \%struct_b, \%c_in);
			}
		}
	}


	####################################
	#    MM-PBSA / MM-GBSA analysis    #
	####################################

	if($c_in{'analyze'} == 1){

		# Read trajectories into array 'traj_files' if 'all' was specified instead of trajectories.
		if(($traj_files[0] eq "all")||($c_in{'trj_file'} eq "all")){
			@traj_files = ();
			my $traj_files_ref = read_trajectory_files($struct_b{$struct_files[0]}, \%c_in);
			@traj_files = @{$traj_files_ref};
		}

		# Check if specified trajectories needed for analysis exist
		my $check_traj = check_traj_files("wamm", \@traj_files, \%struct_b, \%c_in);
		if($check_traj != 1){
			exit;
		}
		
		# When imaging is requested, check if imaged trajectories are already present
		if($c_in{'image'} == 1){
			check_imaged_traj_presence("wamm", \@traj_files, \%struct_b, \%c_in);
		}
		
		# Check provided residue number of receptor structure using receptor structure prepared for first ligand
		# This check cannot be performed for MM-PBSA calculations based on explicit membrane simulations, where the no. of
		# residues in the receptor structure is equivalent to rec_res + membrane lipids + ions + water after initial setup.
		if((! $c_in{'ImplMem'})||($c_in{'ImplMem'} == 0)){
			my $rec_file = $c_in{'root_path'}."/leap/".$struct_b{$struct_files[0]}."/".$struct_b{$struct_files[0]}."_rec.pdb";
			if(! -e $rec_file){
				print "\nWARNING: Could not find receptor file $rec_file\n";
				print "for checking consistency of receptor residue information.\n";
			}
			else{
				my ($check_rec_resno, $ca_no, $resno) = check_resno($rec_file, $c_in{'rec_res'});
				if($check_rec_resno == 0){
					print "\nWARNING: $ca_no C-alpha atoms and $resno residues were detected in the\n";
					print "receptor structure ".$rec_file.".\n";
					print "The number of residues differs from the one specified in the command-file.\n"; 
					print "Please check thoroughly whether there is an error in the no. of residues\n";
					print "specified for the receptor structure.\n";
				}
			}
		}
	
		
		# If PB=2 or PB=3 requested set GB explicitly to 0
		# Decomposition is currently only possible with mbondi-radii (PB=4 and GB=1).
		if((($c_in{'PB'} == 2) || ($c_in{'PB'} == 3))&&($c_in{'GB'} > 0)){
			print "WARNING : PB ".$c_in{'PB'}." cannot be run in combination with GB.\n";
			print "GB will be set to zero in calculation setup to avoid inconsistencies.\n";
			$c_in{'GB'} = 0;
		}
		
		# Set relative path for batch submission to root_path if not batch-specific
		# path is specified.
		if((!$c_in{'pbs_m_path'})||($c_in{'pbs_m_path'} eq "")){
			$c_in{'pbs_m_path'} = $c_in{'root_path'};
		}				

		# Generate directory for storing analysis results
		my $analyze_dir;
		my $chrg_alias;
		my %qsub;
		my $rec_img = 0;

		#START LW_03_2018
		# for prot_prot, no charge set is required
		# use "p" instead of charge method for directories
		if($c_in{'prot_prot'} == 1){
			$c_in{'chrg_meth'} = "p";
			$analyze_dir = $c_in{'root_path'}."/calc_p_".$c_in{'traj'}."t";
			$chrg_alias = "p";
		}
		
		if(($c_in{'chrg_meth'} eq "resp")&&($c_in{'prot_prot'} != 1)){
		#END LW_03_2018
			$analyze_dir = $c_in{'root_path'}."/calc_r_".$c_in{'traj'}."t";
			$chrg_alias = "r";
		}
		#START LW_03_2018
		if(($c_in{'chrg_meth'} eq "am1")&&($c_in{'prot_prot'} != 1)){
		#END LW_03_2018
			$analyze_dir = $c_in{'root_path'}."/calc_a_".$c_in{'traj'}."t";
			$chrg_alias = "a";
		}

		if(! -e $analyze_dir){
			mkdir($analyze_dir);
		}

		foreach my $s (@struct_files){

			if(! -e "$analyze_dir/".$struct_b{$s}){
				mkdir("$analyze_dir/".$struct_b{$s});
			}
			
			##########################
			# Prepare topology files #
			##########################
			
			my $topo_dir = "$analyze_dir/".$struct_b{$s}."/topo";
			if(! -e $topo_dir){
				mkdir($topo_dir);
			}

			# If normal MD simulation with crystal water
			if(($c_in{'wat'} == 1)&&($c_in{'ImplMem'} == 0)){
				delete_cryst_wat($topo_dir, $struct_b{$s}, \%c_in);
			}
			
			# If membrane shall be considered implicitly in MM-PBSA analysis
			#START LW_03_2018
			if($c_in{'prot_prot'} != 1){
				if(($c_in{'ImplMem'} == 1)&&(exists $c_in{'memResNo'})&&($c_in{'memResNo'} != 0)){
					delete_lipids_ions_wat($topo_dir, $struct_b{$s}, \%c_in);
				}
			}
			else{
				prot_delete_lipids_ions_wat($struct_path, $topo_dir, $struct_b{$s}, \%c_in);
			}
			#END LW_03_2018
			

			# Generate leap input script for topology setup for calculations without water
			#START LW_03_2018
			if($c_in{'prot_prot'} != 1){
				gen_leap_in_top_setup($topo_dir, $struct_b{$s}, \%c_in);
			}
			else{
				prot_gen_leap_in_top_setup($topo_dir, $struct_b{$s}, \%c_in);
			}
			#END LW_03_2018

			my $topo_gen_log = $topo_dir."/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".log";
			my $topo_gen_c = "tleap -f $topo_dir/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".in > "
			                 .$topo_gen_log." 2>&1";
			call_prog($topo_gen_c, $topo_gen_log);

			# Check warnings
			print "\nRunning Leap to setup topology files for analysis for $struct_b{$s}.\n";
			check_leap($topo_dir."/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".log");
			unlink("./leap.log");
			
			# Check if created topology files exist and are not empty
			foreach my $tag ("com", "rec", "lig"){
				my $topo_file = $topo_dir."/".$struct_b{$s}."_gb".$c_in{'GB'}."_pb".$c_in{'PB'}."_".$tag.".top";
				if((! -e $topo_file)||(-z $topo_file)){
					print "\nERROR: At least one topology file was not correctly created. Please look\n";
					print "in the file $topo_dir/leap_topo_gb".$c_in{'GB'}."_pb".$c_in{'PB'}.".log for a detailed error messarge.\n";
					exit;
				}
			}

			# Copy topology files of solvated systems to topology folder (for PB=2)
			if($c_in{'PB'} == 2){
				copy_topWAT($topo_dir, $struct_b{$s}, \%c_in);
			}

			######################################
			# Prepare MM-PBSA / MM-GBSA analysis #
			######################################
			
			if((!$c_in{'coord_templ'})||($c_in{'coord_templ'} eq "")){
				$c_in{'coord_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/extract_snaps.in";
			}

			if((-e $c_in{'coord_templ'})&&($c_in{'get_snaps'} == 1)){

				# Coordinate extraction
				#########################

				# Generate directory for storing snapshots - Whole trajectory
				my $snaps_dir = $analyze_dir."/".$struct_b{$s}."/snapshots";

				if(! -e $snaps_dir){
					mkdir($snaps_dir);
				}

				# Create separate directory for storing snapshots with water if PB=2
				my $solv_snaps_dir = "";

				if($c_in{'PB'} == 2){
					$solv_snaps_dir = $analyze_dir."/".$struct_b{$s}."/s_Hyd";

					if(! -e $solv_snaps_dir){
						mkdir($solv_snaps_dir);
						mkdir("$solv_snaps_dir/com");
						mkdir("$solv_snaps_dir/rec");
						mkdir("$solv_snaps_dir/lig");
					}
				}

				# Retrieve atom numbers from complex pdb file
				#START LW_03_2018
				my ($lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, $total_atoms);

				if($c_in{'prot_prot'} != 1){
					($lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, $total_atoms) = identify_atom_no($struct_b{$s}, \%c_in);
				}
				else{
					($lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, $total_atoms) = prot_identify_atom_no($struct_path, $struct_b{$s}, \%c_in);
				}
				#END LW_03_2018

				# Regard difference in total atoms of 1- and 3-trajectory approach
				my %total_atoms;
				$total_atoms{'com'} = $total_atoms;
				if($c_in{'traj'} == 3){
					$total_atoms{'rec'} = determine_total_atom_no("rec", $struct_b{$s}, \%c_in);
					$total_atoms{'lig'} = determine_total_atom_no("lig", $struct_b{$s}, \%c_in);
				}

				# Set path to trajectory files
	 			my @tags;
				if($c_in{'traj'} == 1){
					push(@tags, "com");
				}
				if($c_in{'traj'} == 3){
					@tags = ("com", "lig", "rec");
				}

				my %traj_path;
				foreach my $tag (@tags){
					if($tag eq "rec"){
						$traj_path{$tag} = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$tag."/prod";
					}
					else{
						$traj_path{$tag} = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod";
					}
	 			}

				# Determine topology files - Same location for 1- and 3-trajectory approach
				my %top;
				for my $tag ("lig", "com", "rec"){
					if($tag ne "rec"){
						$top{$tag} = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_".$tag.".top";
					}
					if($tag eq "rec"){
						$top{$tag} = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$tag."/cryst/".$tag."_solv.top";
					}
				}

				# Perform imaging of solute if requested
				my $img_res;
				my $tag_name;
				
				if($c_in{image} == 1){
					print "\nImaging trajectories for ".$struct_b{$s}."\n";
				}

				foreach my $tag (@tags){

					if(($tag eq "rec")&&($rec_img == 1)){
						last;
					}

					# Determine residues to image
					if($tag eq "com"){
						if($c_in{'ImplMem'} == 1){
							$img_res = $c_in{'rec_res'};
						}
						else{
							$img_res = $c_in{'rec_res'}+1;
						}
						$tag_name = "complex";
					}
					elsif($tag eq "rec"){
						$img_res = $c_in{'rec_res'};
						$tag_name = "receptor";
					}
					else{
						$img_res = 1;
						$tag_name = "ligand";
					}
					
					for(my $f=1; $f<=@traj_files; $f++){
						# Check if imaged file is not already present
						my @traj_file_s = split(/\./, $traj_files[$f-1]);
						my $imaged_file = $traj_path{$tag}."/".$traj_file_s[0]."_img.mdcrd.gz";

						# If imaging shall be performed although imaged file already exists, the imaged file will be
						# renamed to avoid overwriting.
						if(($c_in{'image'} == 1)&&($c_in{'use_imaged'} != 1)&&((-e $imaged_file)||(-l $imaged_file))){
							my $sub_imaged = substr($imaged_file, 0, -9);
							my $count_imaged = 1;
							my $new_imaged_traj = $sub_imaged."_".$count_imaged.".mdcrd.gz";
							while(-e $new_imaged_traj){
								$count_imaged++;
								$new_imaged_traj = $sub_imaged."_".$count_imaged.".mdcrd.gz";
							}
							rename($imaged_file, $new_imaged_traj);
							print "Imaged trajectory\n";
							print "$imaged_file\n";
							print "already exists. It will be renamed to\n";
							print "$new_imaged_traj\n";
							print "to avoid overwriting. If you wish to consider the old imaged\n";
							print "trajectory for snapshot extraction, please rename the imaged\n";
							print "trajectory and set the flag 'use_imaged_trajectories' to 1.\n";
						}
						# Perform imaging
						if(($c_in{'image'} == 1)&&(! ((-e $imaged_file)||(-l $imaged_file)))){
							my ($image_c, $image_log, $cpptraj_file, $imaged_trj) = image_traj($traj_files[$f-1], $f, $img_res, $traj_path{$tag}, $top{$tag}, $c_in{'image_backwards'});
							print "Imaging trajectory $f of $tag_name.\n";
							call_prog($image_c, $image_log);
							#unlink($cpptraj_file);
							unlink($image_log);
							my $gzip_log = $c_in{'root_path'}."/gzip.log";
 							my $gzip_c = "gzip $imaged_trj 2> $gzip_log";
 							call_prog($gzip_c, $gzip_log);
							unlink($gzip_log);
						}
						elsif(($c_in{'image'} == 1)&&(((-e $imaged_file)||(-l $imaged_file))&&($c_in{'use_imaged'} == 1))){
							print "Existing imaged trajectory\n";
							print $imaged_file."\n";
							print "will be used for analysis.\n";
						}
						else{
							next;
						}
					}
					if($tag eq "rec"){
						$rec_img = 1;
					}
				}
				
				# Change command-file for snapshot extraction with mm_pbsa.pl
				#START LW_03_2018
				if($c_in{'prot_prot'} != 1){
					modify_extract_snaps($snaps_dir, $lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, \%total_atoms, $struct_b{$s}, \@traj_files, \%c_in);
				}
				else{
					prot_modify_extract_snaps($snaps_dir, $lig_atom_start, $lig_atom_end, $rec_atom_start, $rec_atom_end, \%total_atoms, $struct_b{$s}, \@traj_files, \%c_in);
				}
				#END LW_03_2018

				# Execute snapshot extraction
				if(($c_in{'get_snaps'} == 1)&&(-e $mmpbsa_pl)){
					chdir($snaps_dir);
					print "\nPerforming snapshot extraction for ".$struct_b{$s}.".\n";

					if($c_in{'traj'} == 1){
						my $extract_snaps_log = "$snaps_dir/extract_coordinates_com.log";
						my $extract_snaps_c = "perl ".$mmpbsa_pl." $snaps_dir/extract_coordinates_com.in"; 
						$extract_snaps_c .= " 1>".$extract_snaps_log." 2>".$extract_snaps_log;
						call_prog($extract_snaps_c, $extract_snaps_log);
					}
					if($c_in{'traj'} == 3){
					
						# Generate complex structures
						for my $tag ("com"){
							my $extract_snaps_log = "$snaps_dir/extract_coordinates_$tag.log";
							my $extract_snaps_c = "perl ".$mmpbsa_pl." $snaps_dir/extract_coordinates_$tag.in";
							$extract_snaps_c .= " 1>".$extract_snaps_log." 2>".$extract_snaps_log;
							call_prog($extract_snaps_c, $extract_snaps_log);
						}
						# Rename complex *rec* in complex *com* structures
						my $snaps_no = ($c_in{'Nstop'} - $c_in{'Nstart'}) + 1;
						for(my $n=1; $n<=$snaps_no; $n++){
							my $q = $snaps_dir."/".$struct_b{$s}."_rec.crd.".$n;
							my $t = $snaps_dir."/".$struct_b{$s}."_com.crd.".$n;
							copy($q, $t) || die "Cannot copy $q to $t\n";
							unlink $q;
						}
						# Generate ligand structures
						for my $tag ("lig"){
							my $extract_snaps_log = "$snaps_dir/extract_coordinates_$tag.log";
							my $extract_snaps_c = "perl ".$mmpbsa_pl." $snaps_dir/extract_coordinates_$tag.in";
							$extract_snaps_c .= " 1>".$extract_snaps_log." 2>".$extract_snaps_log;
							call_prog($extract_snaps_c, $extract_snaps_log);
						}
						# Rename ligand *rec* in Ligand *lig* structures
						my $snaps_no = ($c_in{'Nstop'} - $c_in{'Nstart'}) + 1;
						for(my $n=1; $n<=$snaps_no; $n++){
							my $q = $snaps_dir."/".$struct_b{$s}."_rec.crd.".$n;
							my $t = $snaps_dir."/".$struct_b{$s}."_lig.crd.".$n;
							copy($q, $t) || die "Cannot copy $q to $t\n";
							unlink $q;
						}
						# Generate receptor structures
						for my $tag ("rec"){
							my $extract_snaps_log = "$snaps_dir/extract_coordinates_$tag.log";
							my $extract_snaps_c = "perl ".$mmpbsa_pl." $snaps_dir/extract_coordinates_$tag.in";
							$extract_snaps_c .= " 1>".$extract_snaps_log." 2>".$extract_snaps_log;
							call_prog($extract_snaps_c, $extract_snaps_log);
						}
					}
				}

				# Extraction of snapshots with water for PB=2
				# if snapshots were not generated before
				if(($c_in{'PB'} == 2)&&(-e $solv_snaps_dir)){
				
					print "Extracting snapshots with water for PB-hybrid model calculation for ".$struct_b{$s}.".\n";
				
					# Change command-file for snapshot extraction
					modify_extract_snaps_solv($solv_snaps_dir, \%total_atoms, $struct_b{$s}, \@tags, \%traj_path, \@traj_files, \%c_in);

					# Execute extraction
					foreach my $tag (@tags){
						if($mmpbsa_pl){
							chdir("$solv_snaps_dir/$tag");
							my $extract_solv_snaps_log = "$solv_snaps_dir/$tag/extract_coordinates_$tag.log";
							my $extract_solv_snaps_c = "perl ".$mmpbsa_pl." ".$solv_snaps_dir."/".$tag."/extract_coordinates_".$tag.".in"; 
							$extract_solv_snaps_c .= " 1>".$extract_solv_snaps_log." 2>".$extract_solv_snaps_log;
							call_prog($extract_solv_snaps_c, $extract_solv_snaps_log);
						}
					}
				}
				
				# Generate PQR-files for implicit membrane analysis with APBS
				if($c_in{'ImplMem'} == 1){
					my $pqr_dir = $analyze_dir."/".$struct_b{$s}."/pqr_snaps";
					if(! -e $pqr_dir){
						mkdir($pqr_dir);
					}
					chdir($pqr_dir);
					
					# Extract pqr-files for each trajectory
					my $snaps_count = 0;
					my $traj_count = 0;
					my $total_snaps_no = 0;
					mkdir($pqr_dir."/tmp");
					
					#START LW_03_2018
					#Case 1: Default 
					if($c_in{'prot_prot'} != 1){
						foreach my $traj (@traj_files){
							$traj_count++;
							my $internal_count = 0;
							my $cpptraj_pqr;
							
							$cpptraj_pqr = generate_cpptraj_rec_img_pqr($pqr_dir, $traj, $traj_count, $struct_b{$s}, \%c_in);
							
							my $cpptraj_pqr_log = $pqr_dir."/pqr_".$traj_count.".log";
							my $pqr_c = "cpptraj -p ".$c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_com.top";
							$pqr_c .= " -i $cpptraj_pqr > $cpptraj_pqr_log 2> $cpptraj_pqr_log";
							call_prog($pqr_c, $cpptraj_pqr_log);
							
							# Change name of PQR files
							my @pqr_files = <tmp/snap*>;
							my $tmp_files_count = 0;
							while(@pqr_files){
								$total_snaps_no++;
								$tmp_files_count++;
								
								if($total_snaps_no <= $c_in{'Nstop'}){
								
									if(($total_snaps_no >= $c_in{'Nstart'})&&(($total_snaps_no % $c_in{'Nfreq'}) == 0)){
										$snaps_count = $snaps_count + 1;
									
										# Copy complex
										my $q_pqr = "tmp/snap.pqr.".$tmp_files_count;
										my $t_pqr = $struct_b{$s}."_com.pqr.".$snaps_count;
										copy($q_pqr, $t_pqr);
										
										# Generate receptor and ligand
										gen_rec_lig_pqr($pqr_dir,  $struct_b{$s}, $t_pqr, $c_in{'rec_res'}, $snaps_count);
									}
								}
								else{
									last;
								}
								shift(@pqr_files);
							}
							unlink glob "tmp/snap*";
						}
					}
					#Case 2: for prot_prot = 1, use crd to create pqr files
					else{
						my $cpptraj_pqr;
						my $snaps_no = ($c_in{'Nstop'} - $c_in{'Nstart'}) + 1;
						$cpptraj_pqr = prot_generate_cpptraj_rec_img_pqr($snaps_dir, $topo_dir, $pqr_dir, $struct_b{$s}, \%c_in);
						my $cpptraj_pqr_log = $pqr_dir."/pqr_".$traj_count.".log";
						my $pqr_c = "cpptraj -i $cpptraj_pqr > $cpptraj_pqr_log 2> $cpptraj_pqr_log";
						call_prog($pqr_c, $cpptraj_pqr_log);
						
						# Change name of PQR files
						my @pqr_files = <tmp/snap*>;
						my $tmp_files_count = 0;
						while(@pqr_files){
							$total_snaps_no++;
							$tmp_files_count++;
							
							if($total_snaps_no <= $c_in{'Nstop'}){
							
								if(($total_snaps_no >= $c_in{'Nstart'})&&(($total_snaps_no % $c_in{'Nfreq'}) == 0)){
									$snaps_count = $snaps_count + 1;
							
									# Copy complex
									my $q_pqr = "tmp/snap.pqr.".$tmp_files_count;
									my $t_pqr = $struct_b{$s}."_com.pqr.".$snaps_count;
									copy($q_pqr, $t_pqr);
							
									# Generate receptor and ligand
									my ($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
									prot_gen_rec_lig_pqr($pqr_dir,  $struct_b{$s}, $t_pqr, $rec_end , $snaps_count);
								}
							}
							else{
								last;
							}
							shift(@pqr_files);
						}
						unlink glob "tmp/snap*";
					}
					#END LW_03_2018
					rmtree("$pqr_dir/tmp");
				}
				chdir($c_in{'root_path'});
			}


			#  Setup analysis
			###################

			if((! $c_in{'mmpbsa_templ'})||($c_in{'mmpbsa_templ'} eq "")){
				$c_in{'mmpbsa_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/mmpbsa.in";
			}

    		if(-e $c_in{'mmpbsa_templ'}){
      			for(my $n=1; $n<=$n_traj_inter; $n++){

					# Generate directory for storing mm_pbsa_results
					my $ana_dir = $analyze_dir."/".$struct_b{$s}."/s".$traj_inter{$n}->{'Start'}."_".$traj_inter{$n}->{'Stop'}."_".$traj_inter{$n}->{'Offset'};

					if(! -e $ana_dir){
						mkdir($ana_dir);
					}

					# Generate directory for mmpbsa_analysis
					my $calc_dir;
					my $rel_calc;
					if((($c_in{'PB'} == 1)||($c_in{'PB'} == 2)||(($c_in{'PB'} == 3)&&($c_in{'ImplMem'} == 0)))&&($c_in{'Decomp'} > 0)){
						print "WARNING: Currently no decomposition is implemented for PB=".$c_in{'PB'}."\n";
						print "Setup will be done without decomposition.\n";
						$c_in{'Decomp'} = 0;
					}
					if((($c_in{'GB'} == 2)||($c_in{'GB'} == 5))&&($c_in{'Decomp'} > 0)){
						print "WARNING: Currently no decomposition is implemented for GB=".$c_in{'GB'}."\n";
						print "Setup will be done without decomposition.\n";
						$c_in{'Decomp'} = 0;
					}
					if($c_in{'Decomp'} == 0){
						$calc_dir = $ana_dir."/pb".$c_in{'PB'}."_gb".$c_in{'GB'};
						$rel_calc = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/s".$traj_inter{$n}->{'Start'}."_"
						            .$traj_inter{$n}->{'Stop'}."_".$traj_inter{$n}->{'Offset'}."/pb".$c_in{'PB'}."_gb".$c_in{'GB'};
					}
					if($c_in{'Decomp'} > 0){
						$calc_dir = $ana_dir."/pb".$c_in{'PB'}."_gb".$c_in{'GB'}."_dec";
						$rel_calc = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/s".$traj_inter{$n}->{'Start'}."_"
									.$traj_inter{$n}->{'Stop'}."_".$traj_inter{$n}->{'Offset'}."/pb".$c_in{'PB'}."_gb".$c_in{'GB'}."_dec";
					}
					if(! -e $calc_dir){
						mkdir($calc_dir);
					}
					
					my $rel_topo = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/topo";
					
					if($c_in{'PB'} != 2){

						# Create calculation command-file in present directory
			  			my $mmpbsa_in = "$rel_calc/mmpbsa.in";
						my $rel_snaps = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/snapshots";
						my $rel_pqr_dir = "";
						if($c_in{'ImplMem'} == 1){
							$rel_pqr_dir = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/pqr_snaps";
						}

						# Change command-file for mmpbsa calculation with mm_pbsa.pl
						modify_mmpbsa_in($calc_dir, $mmpbsa_in, $struct_b{$s}, $n, \%traj_inter, \@traj_files, \%c_in);

						if(-e $c_in{'pbs_mmpbsa'}){
							gen_pbs_script($calc_dir, $mmpbsa_in, $rel_calc, $struct_b{$s}, $struct_b{$s}, $rel_snaps, $rel_pqr_dir, $rel_topo, $current_dir, \%c_in);
						}

						# Store information for submission
						if((($c_in{'pbs_mmpbsa'})&&(-e $c_in{'pbs_mmpbsa'}))&&(-e $mmpbsa_pl)){
							$qsub{$n}->{$struct_b{$s}}->{'calc'}->{'calc'} = "qsub ".$c_in{'pbs_m_path'}."/".$rel_calc."/run_mmpbsa.pbs\n";
						}
					}

					# Peform seperate analysis if PB=2
					if($c_in{'PB'} == 2){

						for my $pb2_tag ("gas", "hyd"){
							my $pb2_calc_dir = $ana_dir."/pb".$c_in{'PB'}."_gb".$c_in{'GB'}."/".$pb2_tag;

							if(! -e $pb2_calc_dir){
								mkdir($pb2_calc_dir);
								mkdir("$pb2_calc_dir/com");
								mkdir("$pb2_calc_dir/rec");
								mkdir("$pb2_calc_dir/lig");
							}

							my $pb2_rel_snaps;
							if($pb2_tag eq "hyd"){
								$pb2_rel_snaps = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/s_Hyd";
							}
							else{
								$pb2_rel_snaps = "calc_".$chrg_alias."_".$c_in{'traj'}."t/".$struct_b{$s}."/snapshots";
							}

							modify_mmpbsa_in_pb2($pb2_calc_dir, $pb2_tag, $struct_b{$s}, $n, \%traj_inter, \@traj_files, \%c_in);

							if(-e $c_in{'pbs_mmpbsa'}){
								for my $tag ("com", "rec", "lig"){
									my $pb2_mmpbsa_in = "$rel_calc/$pb2_tag/$tag/mmpbsa.in";
									my $job_name = $struct_b{$s}."_".$pb2_tag."_".$tag;
									gen_pbs_script("$pb2_calc_dir/$tag", $pb2_mmpbsa_in, "$rel_calc/$pb2_tag/$tag", $struct_b{$s}, $job_name, $pb2_rel_snaps, "", $rel_topo, $current_dir, \%c_in);
									$qsub{$n}->{$struct_b{$s}}->{$pb2_tag}->{$tag} = "qsub ".$c_in{'pbs_m_path'}."/".$rel_calc."/".$pb2_tag."/".$tag."/run_mmpbsa.pbs\n";
								}
							}
						}
					}

					# Generate input script for implicit membrane calculation with APBS
					if($c_in{'ImplMem'} == 1){
						gen_mmpbsa_few_mem_input($calc_dir, $current_dir, $struct_b{$s}, $traj_inter{$n}->{'Start'}, $traj_inter{$n}->{'Stop'}, $traj_inter{$n}->{'Offset'}, \%c_in);
					}
				}
			}
		}

		# Write qsub_scripts for submission based on selected intervall, offset and PB/GB procedure
		for(my $n=1; $n<=$n_traj_inter; $n++){
			my $qsub_file = $analyze_dir."/qsub_s".$traj_inter{$n}->{'Start'}."_".$traj_inter{$n}->{'Stop'}."_".$traj_inter{$n}->{'Offset'}."_pb".$c_in{'PB'}."_gb".$c_in{'GB'}.".sh";
			open(QSUB, ">$qsub_file") || die "Cannot open $qsub_file for writing.\n";

			foreach my $s (@struct_files){
				foreach my $pb2_tag (keys %{$qsub{$n}->{$struct_b{$s}}}){
					foreach my $tag (keys %{$qsub{$n}->{$struct_b{$s}}->{$pb2_tag}}){
						print QSUB $qsub{$n}->{$struct_b{$s}}->{$pb2_tag}->{$tag};
					}
				}
			}
			close QSUB;
			chmod 0755, $qsub_file;
		}
	}
}

1;
