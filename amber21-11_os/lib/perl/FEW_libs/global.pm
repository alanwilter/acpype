# Association of more intuitive keywords with old keywords used by FEW
sub keyword_association{
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my %new_to_old_keywords;
	
	# Association list
	%new_to_old_keywords = (
		"lig_struct_path"  =>  "i_path",
		"output_path"  =>  "root_path",
		"lig_format_sdf"  =>  "sdf",
		"lig_format_mol2"  =>  "mol2",
		"rec_structure"  =>  "rec",
		"multi_structure_lig_file"  =>  "fi_lig",
		"bound_rec_structure"  =>  "rec_compl",
		"membrane_file" => "membrane",
		"water_in_rec"  =>  "wat",
		
		"structure_separation"  =>  "sep",
		
		"prepare_leap_input"  =>  "leap",
		"non_neutral_ligands"  =>  "chg",
		"lig_charge_file"  =>  "c_file",
		"am1_lig_charges"  =>  "am1",
		"resp_lig_charges"  =>  "resp",
		"calc_charges"  =>  "calc_c",
		"resp_setup_step1"  =>  "r_1",
		"resp_setup_step2"  =>  "r_2",
		"prepare_gauss_batch_file"  =>  "pbs_g",
		"gauss_batch_template"  =>  "pbs_t",
		"gauss_batch_path"  =>  "pbs_w",
		"average_charges"  =>  "avg_chrg",
		"prepare_membrane" => "withMem",
		"ligand_water_cutoff" => "lig_cutoff",
		
		"additional_library"  =>  "add_lib",
		"additional_frcmod"  =>  "add_frcmod",
		"use_lipid14_ff" => "Lipid_ff",
		"use_gaff_lipid_ff" => "Lipid_gaff",
		"SSbond_file"  =>  "cys_f",
		"no_of_rec_residues"  =>  "rec_res",
		"protein1_res_range" =>  "rec_range", # Added by LW_03_2018
		"protein2_res_range" =>  "lig_range", # Added by LW_03_2018
		"charge_method"  =>  "chrg_meth",
		"image_trajectories"  =>  "image",
		"trajectory_files"  =>  "trj_file",
		"protein_protein_com" => "prot_com", # Added by LW_03_2018
		
		"setup_MDsimulations"  =>  "sim",
		"traj_setup_method"  =>  "traj",
		"MD_am1"  =>  "P_am1",
		"MD_resp"  =>  "P_resp",
		"MD_batch_path"  =>  "MD_path",
		"MDequil_template_folder"  =>  "equi",
		"total_MDequil_time"  =>  "equi_t",
		"MDequil_batch_template"  =>  "pbs_e",
		"restrain_membrane_residues" => "rstMem",
		"MDprod_template"  =>  "prod",
		"total_MDprod_time"  =>  "prod_t",
		"MDprod_batch_template"  =>  "pbs_p",
		"restart_file_for_MDprod"  =>  "MD_rst",
		
		"mmpbsa_calc"  =>  "analyze",
		"1_or_3_traj"  =>  "traj",
		"protein_protein"  =>  "prot_prot", # Addeed by LW_03_2018
		"extract_snapshots"  =>  "get_snaps",
		"use_imaged_trajectories" => "use_imaged",
		"image_mass_origin" => "image_backwards",
		"snap_extract_template"  =>  "coord_templ",
		"first_snapshot"  =>  "Nstart",
		"last_snapshot"  =>  "Nstop",
		"offset_snapshots"  =>  "Nfreq",
		"mmpbsa_template"  =>  "mmpbsa_templ",
		"decomposition"  =>  "Decomp",
		"total_no_of_intervals"  =>  "tot_inter",
		"first_PB_snapshot"  =>  "Start",
		"last_PB_snapshot"  =>  "Stop",
		"offset_PB_snapshots"  =>  "Offset",
		"mmpbsa_batch_template"  =>  "pbs_mmpbsa",
		"mmpbsa_batch_path"  =>  "pbs_m_path",
		"mmpbsa_sander_exe"  =>  "sander_exe",
		"parallel_mmpbsa_calc"  =>  "parallel",
		
		"membrane_residue_no" => "memResNo",
		"implicit_membrane" => "ImplMem",
		"apbs_executable" => "apbs_exe",
		"epsilon_solute" => "indi",
		"bottom_membrane_boundary" => "zmem",
		"membrane_thickness" => "tmem",
		"membrane_dielc" => "dielc_mem",
		"second_slab_thickness" => "t_sec_slab",
		"second_slab_dielc" => "dielc_sec_slab",
		"third_slab_thickness" => "t_third_slab",
		"third_slab_dielc" => "dielc_third_slab",
		"ion_concentration" => "ionconc",
		"upper_exclusion_radius" => "r_top",
		"lower_exclusion_radius" => "r_bottom",
		"do_focussing" => "focus",
		"size_large_grid" => "glen_l",
		"size_medium_grid" => "glen_m",
		"size_small_grid" => "glen_s",
		"grid_dimensions" => "dime",
		"nonpolar_solv" => "npsolv", # Added by LW_03_2018
		"lie_calc"  =>  "analyze",
		"lie_executable"  =>  "lie_exe",
		"lie_batch_template"  =>  "pbs_lie",
		"lie_batch_path"  =>  "pbs_l_path",
		"snaps_per_trajectory"  =>  "trj_snaps",
		"first_lie_snapshot"  =>  "start",
		"last_lie_snapshot"  =>  "stop",
		"offset_lie_snapshots"  =>  "offset",
		"sander_executable"  =>  "sander_exe",
		"parallel_lie_call"  =>  "parallel_call",
		"delete_lie_trajectories"  =>  "del_trajs",
		
		"ti_simulation_setup"  =>  "ti_sim",
		"lig_name_v0_struct"  =>  "Lname_start",
		"lig_name_v1_struct"  =>  "Lname_end",
		"lig_alias_v0"  =>  "Lalias_start",
		"lig_alias_v1"  =>  "Lalias_end",
		"softcore_mask_v0"  =>  "mask0",
		"softcore_mask_v1"  =>  "mask1",
		"prepare_match_list"  =>  "match",
		"prepare_inpcrd_prmtop"  =>  "crd_top",
		"lig_inpcrd_v0"  =>  "rst_lig",
		"com_inpcrd_v0"  =>  "rst_com",
		"lig_prmtop_v0"  =>  "top_lig",
		"com_prmtop_v0"  =>  "top_com",
		"match_list_file"  =>  "match_list",
		"create_sybyl_mol2"  =>  "sybyl_mol2",
		"ti_batch_path"  =>  "pbs_t_path",
		"ti_prod_template"  =>  "prod_templ",
		"no_shake"  =>  "noshake",
		"ti_equil"  =>  "ti_equi",
		"ti_equil_batch_template"  =>  "pbs_equi_t",
		"ti_equil_lambda"  =>  "equi_lambda",
		"ti_equil_template"  =>  "equi_templ",
		"ti_production"  =>  "ti_prod",
		"ti_prod_lambda"  =>  "prod_lambda",
		"total_ti_prod_time"  =>  "prod_time",
		"ti_prod_batch_template"  =>  "pbs_prod_t",
		"converge_check_script"  =>  "conv_exe",
		"converge_check_method"  =>  "conv_meth",
		"converge_error_limit"  =>  "error_limit",
		"dVdL_calc_source"  =>  "dVdL_src",
		"ddG_calc_method"  =>  "calc_meth",
	);
	
	# Reset keywords
	foreach my $c_in_key (keys %c_in){
		if(exists $new_to_old_keywords{$c_in_key}){
			$c_in{$new_to_old_keywords{$c_in_key}} = $c_in{$c_in_key};
			delete($c_in{$c_in_key});
		}
	}
	
	return \%c_in;
}


# Set default values
sub set_defaults{
	my $current_dir = shift;
	my $procedure = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	# Set AMBERHOME path
	$c_in{'amberhome'} =  $ENV{"AMBERHOME"};
	
	# If the respective keywords are not provided in the command file:
	
	# Set 'image' to 1
	if((($c_in{'get_snaps'})&&($c_in{'get_snaps'}==1))&&(! exists $c_in{'image'})){
		$c_in{'image'} = 1;
	}
	# Set 'mol2' to 1 and 'sdf' to 0
	if(! exists $c_in{'mol2'}){
		$c_in{'mol2'} = 1;
	}
	if(! exists $c_in{'sdf'}){
		$c_in{'sdf'} = 0;
	}
	
	# Set 'sep' to 0
	if(! exists $c_in{'sep'}){
		$c_in{'sep'} = 0;
	}
	
	# Set 'calc_c' to 1
	if(! exists $c_in{'calc_c'}){
		$c_in{'calc_c'} = 1;
	}
	
	# Set 'pbs_g' to 0
	# No gaussian batch script will be setup per default
	if(! exists $c_in{'pbs_g'}){
		$c_in{'pbs_g'} = 0;
	}

	# Set 'withMem' to 0
	if(! exists $c_in{'withMem'}){
		$c_in{'withMem'} = 0;
	}

	# Set 'withMem' to 1.0 A
	if(($c_in{'withMem'} == 1)&&((! exists $c_in{'lig_cutoff'})||($c_in{'lig_cutoff'}))){
		$c_in{'lig_cutoff'} = 1.0;
	}
	
	# Set 'Lipid_ff' to 1
	# Lipid14 force field will be used by default for membrane
	# simulation setup
	if(! exists $c_in{'Lipid_ff'}){
		$c_in{'Lipid_ff'} = 1;
	}
	# Set 'Lipid_gaff' to 1
	# GaffLipid force field will by default not be used for membrane
	# simulation setup
	if(! exists $c_in{'Lipid_gaff'}){
		$c_in{'Lipid_gaff'} = 0;
	}

        # Set 'water_model' to TIP3P
        # TIP3P water model will be default, if no water model is specified
        # by the user.
        if(! exists $c_in{'water_model'}){
                $c_in{'water_model'} = "TIP3P";
        }
	
	# Set 'equi' folder to example 'equi' folder provided with FEW
	# and corresponding equi_t to 400 ps, if solvated protein shall be
	# simulated 
	# Not in case of "TI", because here only structure preparation 
	# should be possible.
	if(((! exists $c_in{'equi'})||($c_in{'equi'} eq ""))&&($procedure ne "TI")
	&&((! exists $c_in{'withMem'})||($c_in{'withMem'} == 0))){
		$c_in{'equi'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi";
		if((! exists $c_in{'equi_t'})||($c_in{'equi_t'} eq "")){
			$c_in{'equi_t'} = 400;
		}
	}
	
	# Set 'equi' folder to example 'equi_with_membrane' folder provided with FEW
	# and corresponding equi_t to 1120 ps, if a protein embedded in an explicit
	# membrane shall be simulated.
	# Not in case of "TI", because here only structure preparation 
	# should be possible.
	if(((! exists $c_in{'equi'})||($c_in{'equi'} eq ""))&&($procedure ne "TI")
	&&((exists $c_in{'withMem'})&&($c_in{'withMem'} == 1))){
		$c_in{'equi'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi_with_membrane";
		if((! exists $c_in{'equi_t'})||($c_in{'equi_t'} eq "")){
			$c_in{'equi_t'} = 1120;
		}
	}
	
	# Set template for MD production to example input script provided with FEW
	if(((! exists $c_in{'prod'})||($c_in{'prod'} eq ""))&&((! exists $c_in{'withMem'})||($c_in{'withMem'} == 0))){
		$c_in{'prod'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/MD_prod.in";
	}

	if(((! exists $c_in{'prod'})||($c_in{'prod'} eq ""))&&((exists $c_in{'withMem'})&&($c_in{'withMem'} == 1))){
		$c_in{'prod'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/MD_prod_membrane.in";
	}
	
	# Set 'use_imaged' to 1
	if((($c_in{'get_snaps'})&&($c_in{'get_snaps'}==1))&&(! exists $c_in{'use_imaged'})){
		$c_in{'use_imaged'} = 1;
	}

	# If 'image_origin_mass' is not specified, set it to 0
	if(! exists $c_in{'image_backwards'}){
		$c_in{'image_backwards'} = 0;
	}

	# Set 'trajectory_files' to consider per default to 'all'
	if(! exists $c_in{'trj_file'}){
		$c_in{'trj_file'} = "all";
	}
	
	# Set 'ImplMem=0' to not consider implicit membrane if not explicitly requested
	if(! exists $c_in{'ImplMem'}){
		$c_in{'ImplMem'} = 0;
	}
	
	# In case of implicit membrane MM-PBSA calculation set the following defaults
	if($c_in{'ImplMem'} == 1){
		if(! exists $c_in{'indi'}){
			$c_in{'indi'} = 1;
		}
		if(! exists $c_in{'zmem'}){
			$c_in{'zmem'} = -18;
		}
		if(! exists $c_in{'tmem'}){
			$c_in{'tmem'} = 36;
		}
		if(! exists $c_in{'ionconc'}){
			$c_in{'ionconc'} = 0.15;
		}
		if(! exists $c_in{'focus'}){
			$c_in{'focus'} = 0,
		}
		if((! exists $c_in{'dielc_mem'})||($c_in{'dielc_mem'} eq "")){
			$c_in{'dielc_mem'} = 2;
		}
		if(! exists $c_in{'t_sec_slab'}){
			$c_in{'t_sec_slab'} = 0;
		}
		if(! exists $c_in{'dielc_sec_slab'}){
			$c_in{'dielc_sec_slab'} = 0;
		}
		if(! exists $c_in{'t_third_slab'}){
			$c_in{'t_third_slab'} = 0;
		}
		if(! exists $c_in{'dielc_third_slab'}){
			$c_in{'dielc_third_slab'} = 0;
		}
	}
	# Set 'calc_sasa' to 0
	# Do not perform SASA calculation in LIE analysis per default
	if(! exists $c_in{'calc_sasa'}){
		$c_in{'calc_sasa'} = 0;
	}
	
	# Set 'sybyl_mol2' to 0
	# Do not perform setup of mol2 file with sybyl atom types per default
	if(! exists $c_in{'sybyl_mol2'}){
		$c_in{'sybyl_mol2'} = 0;
	}
	
	# Set 'noshake' to 0
	# It is assumed per default that shake will be performed in TI simulations
	if(! exists $c_in{'noshake'}){
		$c_in{'noshake'} = 0;
	}
	
	# Set template for TI equilibration to example input script provided with FEW
	if((! exists $c_in{'equi_templ'})||($c_in{'equi_templ'} eq "")){
		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			$c_in{'equi_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi_TI.in";
		}
		if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
			$c_in{'equi_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi_TI_pmemd.in";
		}
	}
	
	# Set template for TI production to example input script provided with FEW
	if((! exists $c_in{'prod_templ'})||($c_in{'prod_templ'} eq "")){
		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			$c_in{'prod_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/MD_prod_TI.in";
		}
		if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
			$c_in{'prod_templ'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/MD_prod_TI_pmemd.in";
		}
	}

	# Set template for TI equilibration batch script to example batch script provided with FEW
	if((! exists $c_in{'pbs_equi_t'})||($c_in{'pbs_equi_t'} eq "")){
		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			$c_in{'pbs_equi_t'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi_TI.pbs";
		}
		if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
			$c_in{'pbs_equi_t'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/equi_TI_pmemd.pbs";
		}
	}

	# Set template for TI production batch script to example batch script provided with FEW
	if((! exists $c_in{'pbs_prod_t'})||($c_in{'pbs_prod_t'} eq "")){
		if((! exists $c_in{'use_pmemd'})||($c_in{'use_pmemd'} == 0)){
			$c_in{'pbs_prod_t'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/prod_TI.pbs";
		}
		if((exists $c_in{'use_pmemd'})&&($c_in{'use_pmemd'} == 1)){
			$c_in{'pbs_prod_t'} = $ENV{"AMBERHOME"}."/lib/perl/FEW_libs/input_info/prod_TI_pmemd.pbs";
		}
	}
	
	# Set convergence check method to 1
	if((! exists $c_in{'conv_meth'})||($c_in{'conv_meth'} eq "")){
		$c_in{'conv_meth'} = 1;
	}
	
	# Set error limit for convergence check to 0.01 kcal/mol
	if((! exists $c_in{'error_limit'})||($c_in{'error_limit'} eq "")){
		$c_in{'error_limit'} = 0.01;
	}
	
	# Set dV/dL source to 0
	# All files generated in TI production will be considered per default
	if(! exists $c_in{'dVdL_src'}){
		$c_in{'dVdL_src'} = 0;
	}
	
	# Set ddG calculation method to 1
	if((! exists $c_in{'calc_meth'})||($c_in{'calc_meth'} eq "")){
		$c_in{'calc_meth'} = 1;
	}
	
	# Set 'backwards' to 0
	# Flag for backwards compartibility. If set to 1, ff99SB will be used
	# instead of ff12SB
	if(! exists $c_in{'backwards'}){
		$c_in{'backwards'} = 0;
	}

	return \%c_in;
}


# Calling programs
sub call_prog{
	my $command = shift;
	my $log = shift;

	system($command) == 0
				or die "\nERROR:
Execution of external program failed. Please consult the log-file
$log 
for more detailed information about the reasons for the program failure.";
	if($? != 0){
		exit;
	}
}


# Check if leap run ended without warnings
sub check_leap{
	my $log_file = shift;
	my $consider_charges = shift;
	my $warn = 0;

	open(LOG, $log_file) || die "Cannot open log-file $log_file of leap.\n";

	my $line_count_after_warning = 0;
	my $detected_warning = 0;
	while(my $log_line = <LOG>){
		chomp($log_line);
		if($detected_warning == 1){
			$line_count_after_warning++;
		}
		if($line_count_after_warning == 2){
			$line_count_after_warning = 0;
			$detected_warning = 0;
		}
		if(($detected_warning == 1)&&($line_count_after_warning == 1)){
			if(($log_line =~ m/charge/)||($log_line =~ m/check/)){
				next;
			}
			elsif($log_line =~ m/residue name to PDB format/){
				next;
			}
			elsif($log_line =~ m/name change in pdb file residue/){
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
			if(($log_line =~ m/charge/)||($log_line =~ m/check/)){
				next;
			}
			elsif($log_line =~ m/residue name to PDB format/){
				next;
			}
			elsif($log_line =~ m/name change in pdb file residue/){
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
		
		if($log_line =~ m/missing parameters/){
			print $log_line . "\n";
			$warn = 1;
		}
		if($log_line =~/usage:  bond <atom1> <atom2> \[order\]/){
			print "\nERROR: Disulfide bridge definition is not correct. Please ensure that the\n";
			print "provided S-S connectivity information is consistent with the residue numbering\n";
			print "generated by LEaP. In order to be able to directly use the residue numbers\n";
			print "present in the provided receptor structure for disulfide bridge definition,\n";
			print "it is required to load the receptor PDB structure into LEaP and re-write it\n";
			print "employing the 'savepdb' command of LEaP.\n";
			exit;
		}
	}

	if($warn == 1){
		print "LEaP run terminated with warnings. Please inspect the LEaP log file at ".$log_file."\n";
	}
	if($warn != 1){
		print "Terminated without warnings.\n";
	}
}

1;
