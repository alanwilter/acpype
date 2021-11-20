# Perl module for checking the provided input data

use strict;

# Check common command file input
sub check_input_common{
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $check = 0;

	my $file_extension;
	if($c_in{'mol2'} == 1){
		$file_extension = "mol2";
	}
	if($c_in{'sdf'} == 1){
		$file_extension = "sdf";
	}

	print "\nChecking input parameters and files.\n";
	
	if($c_in{'amberhome'} eq ""){
		print "\nERROR: The AMBERHOME variable is not set on your system.\n";
		print "Please set the AMBERHOME variable according to the instructions\n";
		print "in the AMBER manual before calling FEW.\n";
	} 
	#Next line modified LW_03_2018
	elsif(((! $c_in{'fi_lig'})||($c_in{'fi_lig'} eq ""))&&($c_in{'prot_prot'} != 1)&&(($c_in{'i_path'} eq "")||(! -e $c_in{'i_path'}))){
		print "\nERROR: Path to folder with input structures not correctly specified.\n\n";
	}
	elsif((!$c_in{'root_path'})||($c_in{'root_path'} eq "")){
		print "\nERROR: No root path for output specified.\n\n";
	}
	elsif(($c_in{'sep'} == 1)&&(!$c_in{'fi_lig'})){
		print "\nERROR: No file containing the structures for requested separation\n";
		print "was specified. Please make sure that you give the base name of the\n";
 	   	print "ligand structure file in the command file.\n\n";
	}
	elsif(($c_in{'sep'} == 1) && ((! -e $c_in{'i_path'}."/".$c_in{'fi_lig'}.".".$file_extension)
	                          &&  (! -e $c_in{'i_path'}."/".$c_in{'fi_lig'}))){
		print "\nERROR: Specified file containing the input structures\n";
		print $c_in{'fi_lig'}." does not exist.\n";
		print "Please make sure you define the correct file in the command file.\n\n";
	}
	elsif(($c_in{'sdf'} == 0)&&($c_in{'mol2'} == 0)){
		print "\nERROR: Format of ligand input file(s) not specified. Please make sure you\n";
		print "use files in either sdf- or mol2-format and give the corresponding format in\n";
		print "the command file.\n\n";
	}
	elsif(($c_in{'sdf'} == 1)&&($c_in{'mol2'} == 1)){
		print "\nERROR: Flags set to 1 for both mol2- and sdf-format. The ligand structures\n";
		print "can be provided only in either mol2- or sdf-format. Please correct your settings\n";
		print "in the command file.\n";
	}
	elsif(($c_in{'leap'} == 1)&&((($c_in{'pbs_g'} == 1)&&($c_in{'resp'} != 1))||(($c_in{'pbs_t'} != "")&&($c_in{'resp'} != 1)))){
		print "\nERROR: You requested preparation of PBS-Input scripts for gaussian runs\n";
	   	print "although the selected charge determination does not require gaussian\n";
		print "calculations. Please make sure that your input is correct!\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&(($c_in{'pbs_g'} == 1)&&(!$c_in{'pbs_t'}))){
		print "\nERROR: The preparation of a pbs-Input script for gaussian was requested.\n";
		print "Rhis procedure requires a template file. Please specify such a file in the\n";
		print "command file\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&(($c_in{'pbs_g'} == 1)&&($c_in{'pbs_t'})&&(! -e $c_in{'pbs_t'}))){
		print "\nERROR: The preparation of a pbs-Input script for gaussian was requested.\n";
		print "This procedure requires a template file. The template file you specified\n";
		print $c_in{'pbs_t'}."\n";
		print "does not exist. Please make sure you give the correct path and name of the\n";
		print "template pbs-script for the gaussian calulation.\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&(($c_in{'r_2'} == 1)&&(! -e $c_in{'root_path'}."/gauss"))){
		print "\nERROR: Setup of files with RESP charges was requested. A prerequisite for this\n";
		print "procedure is the existens of ESP data in form of gaussian output files. Please make\n";
		print "sure you calculate the ESP (resp_setup_step1 keyword) before running the setup.\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&(($c_in{'r_1'} == 1)&&($c_in{'r_2'} == 1))){
		print "\nERROR: Steps 1 & 2 of the RESP charge calculation procedure cannot be performed\n";
		print "in the same step. Please run step 1 'ESP calculation with gaussian' before you\n";
		print "start the second step of the calculation.\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&(($c_in{'am1'} == 1)&&($c_in{'resp'} == 1))){
		print "\nERROR: Setup for AM1-BCC charges and RESP chargers cannot be performed\n";
		print "simultaneously.\n\n";
	}
	elsif(($c_in{'leap'} == 1)&&((($c_in{'am1'} == 1)&&($c_in{'r_1'} == 1))
						 	   ||(($c_in{'am1'} == 1)&&($c_in{'r_2'} == 1)))){  
		print "\nERROR: Calculation of AM1-BCC charges and RESP charges cannot be performed\n";
		print "simultaneously. Please request either AM1-BCC charge calculation or RESP\n";
		print "charge calculation.\n\n";
	}
	elsif((($c_in{'leap'} == 1)&&($c_in{'chg'} == 1))&&((! $c_in{'c_file'})||(! -e $c_in{'c_file'}))){
		print "\nERROR: File in which total charges for non-neutral ligands are defined\n";
		print "not correctly specified.\n\n";
	}
	elsif(($c_in{'avg_chg'} ne "")&&(! -e $c_in{'avg_chg'})){
		print "\nERROR: The specified file containing the isomer pairs for charge averaging\n";
		print $c_in{'avg_chg'}." does not exist.\n\n";
	}
	elsif((($c_in{'rec_compl'})&&($c_in{'rec_compl'} ne ""))&&((! $c_in{'rec'})&&($c_in{'rec'} eq ""))){
		print "\nERROR: If a seperate structure of the bound conformation of the receptor is provided\n";
		print "(keyword 'bound_rec_structure') an apo structure must also be given as 'rec_structure'.\n";
		print "Please ensure that you define receptor PDB structures of both the unbound and the bound\n";
		print "state in the command file or specify only a PDB structure of the bound state after the\n";
		print "keyword 'rec_structure'.\n\n";
	}
	elsif(($c_in{'membrane'})&&($c_in{'membrane'} ne "")&&(! -e $c_in{'membrane'})){
		print "\nERROR: The file with the membrane, ion, and water information you specified in the\n";
		print "command file does not exist. Please provide the correct path and name of the membrane file.\n";
	}
	elsif(($c_in{'leap'} == 1)&&((!$c_in{'rec'})||(! -e $c_in{'rec'}))){
		print "\nERROR: Receptor PDB file not correctly specified.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&($c_in{'equi'})&&($c_in{'equi'} ne "")&&(! -e $c_in{'equi'})){
		print "\nERROR: The folder with template input files for the equilibration phase does not\n";
		print "exist. Please ensure you specified the correct path.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&($c_in{'pbs_e'})&&($c_in{'pbs_e'} ne "")&&(! -e $c_in{'pbs_e'})){
		print "\nERROR: The specified batch template file for the equilibration phase of the\n";
		print "MD simulations does not exist. Please ensure you specified the file correctly.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&($c_in{'lib'})&&(! -e $c_in{'lib'})){
		print "\nERROR: The additional library file you specified\n";
		print $c_in{'lib'}." does not exist.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&((! $c_in{'rec_res'})||($c_in{'rec_res'} eq ""))){
		print "\nERROR: Number of residues in receptor structure not given in command file.\n";
		print "Please specify the number of residues in the receptor PDB structure.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&($c_in{'P_am1'} == 0)&&($c_in{'P_resp'} == 0)){
		print "\nERROR: Types of charges that shall be considered in the simulation setup\n\n";
		print "not specified. Please set either 'MD_am1' or 'MD_resp' to one.\n\n";
	}
	elsif(($c_in{'sim'} == 1)&&(exists $c_in{'cys_f'})&&($c_in{'cys_f'} ne "")&&(! -e $c_in{'cys_f'})){
		print "\nERROR: The 'SSbond_file' with disulfide bridge information that is specified\n";
		print "in the command file does not exist. Please check that the absolute path and name\n";
		print "of the file are correct.\n";
	}
	elsif(($c_in{'sim'} == 1)&&($c_in{'withMem'} == 1)&&($c_in{'traj'} != 1)){
		print "\nERROR: Setup of MD simulations with membrane is currently only possible according\n";
		print "to the 1-trajectory approach.\n";
	}
	elsif(($c_in{'add_lib'})&&($c_in{'add_lib'} ne "")&&(! -e $c_in{'add_lib'})){
		print "\nERROR: The additional library you specified in the command file\n";
		print $c_in{'add_lib'}." does not exist.\n";
		print "Please correct the specification of the additional library in the command file.\n\n";
	}
	elsif(($c_in{'add_frcmod'})&&($c_in{'add_frcmod'} ne "")&&(! -e $c_in{'add_frcmod'})){
		print "\nERROR: The additional parameter (frcmod) file you specified in the command file\n";
		print "does not exist. Please correct the specification of the additional parameter\n";
		print "file in the command file.\n\n";
	}
	elsif((($c_in{'Lipid_gaff'})&&($c_in{'Lipid_gaff'} == 1))&&((! $c_in{'add_lib'})||($c_in{'add_lib'} eq ""))){
		print "\nERROR: Usage of the GAFFlipid force field was requested for the setup of the MD simulations.\n";
		print "For using this force field an additional library file with the specifications for the respective\n";
		print "lipids needs to be provided. It can be obtained from http://lipidbook.bioch.ox.ac.uk\n";
	}
	elsif((($c_in{'Lipid_gaff'})&&($c_in{'Lipid_gaff'} == 1))&&((! $c_in{'add_frcmod'})||($c_in{'add_frcmod'} eq ""))){
		print "\nERROR: Usage of the GAFFlipid force field was requested for the setup of the MD simulations.\n";
		print "For using this force field an additional parameter file containing the parameters for the respective\n";
		print "lipids needs to be provided. It can be obtained from http://lipidbook.bioch.ox.ac.uk\n";
	}
	elsif(((($c_in{'leap'})&&($c_in{'leap'} == 1))||(($c_in{'sep'})||($c_in{'sep'} == 1)))&&
		 (($c_in{'sdf'} == 1)&&((! $c_in{'fi_lig'})||($c_in{'fi_lig'} eq "")))){
		print "\nERROR: If ligands are in sdf-format a multi-structure file must be provided.\n";
		print "Base name of sdf-file with provided structures is not given.\n\n";
	}
	elsif(((($c_in{'leap'})&&($c_in{'leap'} == 1))||(($c_in{'sep'})||($c_in{'sep'} == 1)))&&
			(($c_in{'fi_lig'})&&($c_in{'fi_lig'} ne ""))&&((! $c_in{'sep'})||( $c_in{'sep'} == 0))){
		print "\nERROR: If a multi-structure ligand file is provided separation of ligand structures\n";
		print "must be requested by 'structure_separation=1'. Please correct the settings in the\n";
		print "command file.\n\n";
	}
	else{
		$check = 1;
	}
	
	# Check receptor
	if($check == 1){
		if(($c_in{'sim'} == 1)||($c_in{'leap'} == 1)){
			if(! -e $c_in{'rec'}){
					print "\nERROR: The receptor file ".$c_in{'rec'}." does not exist.\n";
					$check = 0;
			}
			else{
				# 1. Structure
				my $rec_struct_passed = check_rec_struct($c_in{'rec'});
				if($rec_struct_passed == 0){
					$check = 0;
				}
		
				if($c_in{'sim'} == 1){
				
					# 2. Number or receptor residues
					my ($check_rec_resno, $ca_no, $resno) = check_resno($c_in{'rec'}, $c_in{'rec_res'});
					if($check_rec_resno == 0){
						print "\nWARNING: $ca_no C-alpha atoms and $resno residues were detected in the\n";
						print "receptor structure ".$c_in{'rec'}.".\n"; 
						print "The number of residues differs from the one specified in the command file.\n"; 
						print "Please ensure that the number of receptor residues specified in the\n";
						print "command file represents the actual number of residues in the receptor\n";
						print "structure.\n";
					}
			
					if(($c_in{'rec_compl'})&&(-e $c_in{'rec_compl'})){
						my ($check_rec_resno, $ca_no, $resno) = check_resno($c_in{'rec_compl'}, $c_in{'rec_res'});
						if($check_rec_resno == 0){
							print "\nWARNING: $ca_no C-alpha atoms and $resno residues were detected in the\n";
							print "receptor structure ".$c_in{'rec_compl'}.".\n";
							print "The number of residues differs from the one specified in the command file.\n"; 
							print "Please ensure that the number of receptor residues specified in the\n";
							print "command file represents the actual number of residues in the receptor\n";
							print "structure.\n";
						}
					}
				}
			}
		}
	}
	
	# Check existence of file with ligand structures in case of multi-structure file
	if($check == 1){
		if((exists $c_in{'fi_lig'})&&($c_in{'fi_lig'} ne "")){
			my $fi_lig_name = "";
			my @fi_lig_components = split(/\./, $c_in{'fi_lig'});
			if(($fi_lig_components[1] eq "mol2")||($fi_lig_components[1] eq "sdf")){
				$fi_lig_name = $fi_lig_components[0].".".$file_extension;
			}
			if($fi_lig_components[1] eq ""){
				$fi_lig_name = $c_in{'fi_lig'}.".".$file_extension;
			}
			
			if(! -e $fi_lig_name){
				print "\nERROR: The multi-structure ligand file could not be found. Based on the information\n";
				print "provided in the command file it is expected that the multi-structure ligand file can\n";
				print "be found under $fi_lig_name\n\n";
				$check = 0;
			}
		}
	}
	
	# Check file containing membrane, ions, and water
	if($check == 1){
		if((($c_in{'sim'} == 1)||($c_in{'leap'} == 1))&&((exists $c_in{'membrane'})&&($c_in{'membrane'} ne ""))){
			if(-e $c_in{'membrane'}){
				# Check format of provided file
				my $check_mem_file_passed = check_mem_file($c_in{'membrane'}, \%c_in);
				if($check_mem_file_passed == 0){
					$check = 0;
				}
			}
			else{
				print "\nERROR: The file with membrane, ion, and water information specified\n";
				print "in the command file could not be found. Based on the information\n";
				print "provided in the command file it is expected that the file with membrane\n";
				print "information can be found under ".$c_in{'membrane'}."\n\n";
			}
		}
	}
	
	
	# Ensure that not two lipid force fields are used at the same time
	if($check == 1){
		if($c_in{'sim'} == 1){
			if(((exists $c_in{'Lipid_ff'})&&(exists $c_in{'Lipid_gaff'}))&&(($c_in{'Lipid_ff'} == 1)&&($c_in{'Lipid_gaff'} ==1))){
				print "WARNING: Setup of MD simulations using the GAFFlipid force field was requested, therefore\n";
				print "the Lipid14 force field will not be used.\n";
				$c_in{'Lipid_ff'} = 0;
			}
			if((($c_in{'withMem'} == 1)&&((! exists $c_in{'Lipid_ff'})||($c_in{'Lipid_ff'} == 0))
										&&((! exists $c_in{'Lipid_gaff'})||($c_in{'Lipid_gaff'} == 0)))){
				print "WARNING: Setup of MD simulations with explicit membrane was requested, but no lipid force field\n";
				print "was specified. The Lipid14 force field will be used per default.\n";
				$c_in{'Lipid_ff'} = 1;
			}
		}
	}
	
	# Check existence of PBS template file for production simulation
	if($check == 1){
		if(($c_in{'sim'} == 1)&&($c_in{'pbs_p'})&&($c_in{'pbs_p'} ne "")&&(! -e $c_in{'pbs_p'})){
			print "\nWARNING: The batch template for production specified in the command file:\n";
			print $c_in{'pbs_p'}."\n";
			print "does not exist. Therefore no batch script for the production will be generated.\n\n";
		}
	}
	
	# Check if the program Babel is available if sdf-format was provided
	if($check == 1){
		if($c_in{'sdf'} == 1){
			my $call = "which babel >log.txt 2>error.txt";
			system($call);
			
			if(-z "log.txt"){
				print "\nERROR: The program 'babel' required for conversion of sdf- to mol2-format\n";
				print "cannot be invoked. Please ensure that 'babel' is available on your system and\n";
				print "can be called by its name\n";
				unlink("log.txt");
				unlink("error.txt");
				exit;
			}
			else{
				unlink("log.txt");
				unlink("error.txt");
			}
		}
	}
	
	return $check;
}


# Check existence of structure specific folder in "leap" directory
sub check_file_existence_for_commonMDsetup{
	my $chrg_meth = shift;
	my $ref_c_in = shift;
	my $ref_struct_b = shift;
	my %c_in = %{$ref_c_in};
	my %struct_b = %{$ref_struct_b};
	my $check = 1;
	
	# Check existence of leap directory and library file
	foreach my $s (keys %struct_b){
		my $leap_dir = $c_in{'root_path'}."/leap/".$struct_b{$s};
		if((! -e $leap_dir)&&($c_in{'leap'}== 0)){
			print "ERROR: Parameter files missing for structure ".$struct_b{$s}.". Please ensure\n";
			print "that for every structure for which a MD simulation shall be prepared library\n";
			print "and parameter files are available in the \'leap\' directory.\n";
			$check = 0;
			last;
		}
		# Check if library file exists that is required for setup of coordinate and topology files
		my $lib_file = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}."_".$chrg_meth.".lib";
		if((! -e $lib_file)||(-z $lib_file)){
			print "ERROR: The library file needed for setup of coordinate and topology files for\n";
			print "MD input, that should have been created in the previous step, cannot be found\n";
			print "under the default location: $lib_file.\n";
			$check = 0;
			last;
		}
		# Check if parameter file exists that is required for setup of coordinate and topology files
		my $parm_file = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}.".frcmod";
		if((! -e $parm_file)||(-z $parm_file)){
			print "ERROR: The parameter file for the structure ".$struct_b{$s}." that should have been\n";
			print "created in the previous step, cannot be found under the default location: $parm_file\n";
			$check = 0;
			last;
		} 
	}
	return $check;
}


# Specific check of input parameters for WAMM module
sub check_input_wamm{

	my $n_traj_inter = shift;
	my $ref_traj_files = shift;
	my @traj_files = @{$ref_traj_files};
	my $ref_traj_inter = shift;
	my %traj_inter = %{$ref_traj_inter};
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $check = 0;
	
	if(($c_in{'sim'} == 1)&&(($c_in{'traj'} != 1)&&($c_in{'traj'} != 3))){
		print "\nERROR: Please select either setup according to the 1-trajectory approach\n";
		print "or according the 3-trajectory approach.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'PB'} == 2)&&($c_in{'traj'} != 3)){
    	print "\nERROR: PB procedure 2 requires 3-trajectory input. Please make sure that\n";
		print "traj=3 if PB=2 is used.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&(($n_traj_inter == 0)&&(! $traj_inter{'1'}->{'Start'})
								 &&(! $traj_inter{'1'}->{'Stop'})&&(! $traj_inter{'1'}->{'Offset'}))){
		print "\nERROR: Please specify the total number of intervals (tot_inter) to regard\n";
		print "in the MM-PBSA or MM-GBSA analysis together with the corresponding 'first_PB_snapshot',\n";
		print "'last_PB_snapshot', and 'offset_PB_snapshots' entries.\n";
	}
	#Next line modified LW_03_2018
	elsif(($c_in{'analyze'} == 1)&&(($c_in{'chrg_meth'} ne "am1")&&($c_in{'chrg_meth'} ne "resp"))&&($c_in{'prot_prot'} != 1)){
		print "\nERROR: Charge method that shall be used for MM-PBSA or MM-GBSA analyses not\n";
		print "correctly specified.\n";
	}	
	elsif(($c_in{'analyze'} == 1)&&($c_in{'Decomp'} == 1)&&(($c_in{'PB'} != 4)||($c_in{'GB'} != 1))&&($c_in{'ImplMem'} == 0)){
		print "\nERROR: Decomposition only works with options PB=4 and GB=1.\n"; 
	}
	elsif(($c_in{'analyze'} == 1)&&(($c_in{'mol2'} == 0)||($c_in{'sdf'} == 1))){
		print "\nERROR: Structures in mol2-format are expected as input structures. At this stage\n";
		print "of analysis structures in mol2-format should be available from the MD setup step.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'get_snaps'}==1)&&(($c_in{'coord_templ'} ne "")&&(! -e $c_in{'coord_templ'}))){
		print "\nERROR: The specified template file for coordinate extraction\n";
		print $c_in{'coord_templ'}." does not exist.\n";
		print "Please ensure that the correct path and name are given for the template file.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'mmpbsa_templ'})&&($c_in{'mmpbsa_templ'} ne "")&&(! -e $c_in{'mmpbsa_templ'})){
		print "\nERROR: The input file for mm_pbsa.pl specified in the command file\n";
		print $c_in{'mmpbsa_templ'}." does not exist.\n";
		print "Please make sure that the correct path and name of the template file are given in\n";
		print "the command file.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'pbs_mmpbsa'})&&($c_in{'pbs_mmpbsa'} ne "")&&(! -e $c_in{'pbs_mmpbsa'})){
		print "\nERROR: The pbs template file for the analysis specified in the command file\n";
		print $c_in{'pbs_mmpbsa'}." does not exist.\n";
		print "Please make sure that the correct path and name of the template file are given in the command file\n";
	}
	#Next LW_03_2018
	elsif(($c_in{'analyze'} == 1)&&((! $c_in{'rec_res'})||($c_in{'rec_res'} eq ""))&&($c_in{'prot_prot'} != 1)){
		print "\nERROR: Information about the number of residues in the receptor is missing. Please specify in\n";
		print "the command file the actual number of residues that shall be regarded as part of the receptor.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'PB'} == 1)&&(($c_in{'GB'} == 2) || ($c_in{'GB'} == 5) || ($c_in{'GB'} == 6))){
		print "\nERROR: PB method 1 can only be run in combination with GB method 1 to avoid radii set inconsistencies.\n";
		print "Please set GB=1 or GB=0, when you use PB=1.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'PB'} == 4)&&(($c_in{'GB'} == 2) || ($c_in{'GB'} == 5) || ($c_in{'GB'} == 6))){
		print "\nERROR: PB method 4 can only be run in combination with GB method 1 to avoid radii set inconsistencies.\n";
		print "Please set GB=1 or GB=0, when you use PB=4.\n";
	}
	# MM-PB/GBSA calculations considering the membrane implicitly
	elsif(($c_in{'analyze'} == 1)&&($c_in{'ImplMem'} == 1)&&(($c_in{'PB'} != 3)||($c_in{'GB'} != 0))){
		print "\nERROR: Implicit membrane MM-PBSA calculations are currently only supported with PB=3 and GB=0.\n";
	}

	#START LW_03_2018
	# MM-PB/GBSA calculations with two proteins (prot_prot) and implicit membrane 
	elsif(($c_in{'prot_prot'} == 1)&&($c_in{'sim'} == 1)){
		print "\nERROR: Setup of Protein-Protein MD simulations is currently not supported.\n";
	}
	elsif(($c_in{'prot_prot'} == 1)&&($c_in{'leap'} == 1)){
		print "\nERROR: Setup of Protein-Protein topology files is currently not supported.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'prot_prot'} == 1)&&($c_in{'ImplMem'} != 1)){
		print "\nERROR: Protein-Protein MM-PBSA calculations are currently only supported with implicit membrane.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'prot_prot'} == 1)&&($c_in{'npsolv'} == 1)&&(($c_in{'Decomp'} == 0) || ($c_in{'Decomp'} eq ""))){
		print "\nERROR: Protein-Protein MM-PBSA calculations with \"nonpolar_solv\" option\n";
		print "are currently only supported for per-residue decomposition of energies\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'prot_prot'} == 1)&&(($c_in{'memResNo'} != 0) || ($c_in{'memResNo'} ne ""))){
		print "\nERROR: Protein-Protein MM-PBSA calculations do not need an input of the membrane residue number.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'prot_prot'} == 1)&&(((! $c_in{'rec_range'}) || ($c_in{'rec_range'} eq ""))||((! $c_in{'lig_range'})||($c_in{'lig_range'} eq "")))){
		print "\nERROR: Information about the number of residues for protein-protein calculation is missing. \n";
		print "Please specify in the command file the actual range of residues that shall be regarded.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&(($c_in{'prot_prot'} != 1) || ($c_in{'ImplMem'} != 1) )&&(($c_in{'npsolv'} != 0)||($c_in{'npsolv'} ne ""))){
		print "\nERROR: \"nonpolar_solv\" option can only be used with \"protein_protein\" mode \n";
		print "and implicit membrane.\n";
	}
	#END LW_03_2018
	else{
		$check = 1;
	}
	
	# Check consistency of provided interval information
	if(($c_in{'analyze'} == 1)&&($n_traj_inter)&&($n_traj_inter != 0)){
		for(my $n=1; $n<=$n_traj_inter; $n++){
			if((! $traj_inter{$n}->{'Start'})||(! $traj_inter{$n}->{'Stop'})||(! $traj_inter{$n}->{'Offset'})){
				print "\nERROR: The number of 'first_PB_snapshot', 'last_PB_snapshot', and 'offset_PB_snapshots'\n";
				print "definitions is not consistent with the number of intervals (tot_inter) that shall be analyzed.\n";
				print "Please ensure that the settings in the command file are consistent.\n";
				$check = 0;
				last;
			}
		}
	}
	
	# Check correct specification of membrane residues
	if($check == 1){
		if((($c_in{'analyze'} == 1)&&($c_in{'ImplMem'} == 1))&&((exists $c_in{'membrane'})&&($c_in{'membrane'} ne ""))){
			if(-e $c_in{'membrane'}){
				# Check format of provided file
				my $check_mem_file_passed = check_mem_file($c_in{'membrane'}, \%c_in);
				if($check_mem_file_passed == 0){
					$check = 0;
				}
			}
		}
	}
	
	# Check correct input for membrane slab regions
	if($check == 1){
		my @dime = ("97", "129", "161"); # Allowed grid dimensions
		#Next line modified LW_03_2018
		if(($c_in{'analyze'} == 1)&&($c_in{'ImplMem'} == 1)&&($c_in{'prot_prot'} != 1)){
			if((! exists $c_in{'r_top'})||($c_in{'r_top'} == 0)){
				print "\nERROR: MM-PBSA calculation with implicit membrane was requested, but\n";
				print "no upper exclusion radius was specified for the region in which the\n";
				print "receptor is inserted into the membrane. Please specify the upper exclusion\n";
				print "radius in the command file.\n";
			}
			if((! exists $c_in{'r_bottom'})||($c_in{'r_bottom'} == 0)){
				print "\nERROR: MM-PBSA calculation with implicit membrane was requested, but\n";
				print "no lower exclusion radius was specified for the region in which the\n";
				print "receptor is inserted into the membrane Please specify the lower exclusion\n";
				print "radius in the command file.\n";
			}
			if((exists $c_in{'dime'})&&(!(grep {$c_in{'dime'} == $_} @dime))){
				print "\nERROR: The specified grid dimensions are not allowed.\n";
				print "Possible values for 'grid_dimensions' are 97, 129, and 161\n";
			}
			if(($c_in{'t_sec_slab'} != 0)&&($c_in{'dielc_sec_slab'} == 0)){
				print "ERROR: Setup of a membrane with a second membrane slab was requested,\n";
				print "but no dielectric constant was provided for this membrane slab region.\n";
				$check = 0;
			}
			elsif(($c_in{'t_third_slab'} != 0)&&($c_in{'t_sec_slab'} == 0)){
				print "ERROR: Setup of a third membrane slab was requested, but no parameters\n";
				print "for the second membrane slab were provided. A third slab can only be\n";
				print "setup if a second slab is also specified.\n";
				$check = 0;
			}
			elsif(($c_in{'t_third_slab'} != 0)&&($c_in{'dielc_third_slab'} == 0)){
				print "ERROR: Setup of a membrane with a third membrane slab was requested,\n";
				print "but no dielectric constant was provided for this membrane slab region.\n";
				$check = 0;
			}
			elsif(($c_in{'t_sec_slab'} + $c_in{'t_third_slab'}) > ($c_in{'tmem'} / 2)){
				print "ERROR: The specified thickness of the membrane slab(s) is larger than\n";
				print "the specified thickness of the membrane. Please ensure that the sum\n";
				print "of the thicknesses of the membrane slabs is not larger than half the\n";
				print "thickness of the membrane bilayer.\n";
				$check = 0;
			}
			else{
				$check = 1;
			}
		}
	}
	#START LW_03_2018
	# Check correct input for protein-protein calculations
	if($check == 1){
		if(($c_in{'analyze'} == 1)&&($c_in{'prot_prot'} == 1)){
			my ($rec_start, $rec_end) = split(/-/, $c_in{'rec_range'});
			my ($lig_start, $lig_end) = split(/-/, $c_in{'lig_range'});
			if(($rec_start eq "")||($rec_end eq "")){
				print "\nERROR: The specified residue range for protein 1 cannot be read.\n";
				print "Please specifiy a range in the style of \"RESNUM1-RESNUM2\"\n";
				$check = 0;
			}
			elsif(($lig_start eq "")||($lig_end eq "")){
				print "\nERROR: The specified residue range for protein 2 cannot be read.\n";
				print "Please specifiy a range in the style of \"RESNUM1-RESNUM2\"\n";
				$check = 0;
			}
			elsif(($rec_end < $rec_start) || (($lig_end < $lig_start))){
				print "\nERROR: The specified residue range cannot be read.\n";
				print "Please check if the range is given from small to big residue number\n";
				$check = 0;
			}
			elsif(($rec_start > $lig_start) || (($lig_end < $lig_end))){
				print "\nERROR: The specified residue range cannot be read.\n";
				print "Please define the range of protein 1 to cover the lower residue number range\n";
				$check = 0;
			}
		}
	}
	#END LW_03_2018

	return $check;
}	


# Specific check of input parameters for LIEW module
sub check_input_liew{

	my $current_FEW_dir = shift;
	my $ref_c_in = shift;
	my $ref_traj_files = shift;
	my %c_in = %{$ref_c_in};
	my @traj_files = @{$ref_traj_files};
	my $check = 0;
	my $traj_no = @traj_files;
	
	if(($c_in{'sim'} == 1)&&($c_in{'traj'} != 3)){
		print "\nERROR: You requested a ".$c_in{'traj'}."-trajectory setup in LIE calculation modus.\n";
		print "Only 3-trajectory approach works for LIE calculations.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&((! $c_in{'chrg_meth'})||($c_in{'chrg_meth'} eq ""))){
		print "\nERROR: Charge method used for determination of ligand charges for molecular dynamics simulations\n";
		print "not correctly specified. Please use either \"resp\" or \"am1\".\n";
	}
	elsif(($c_in{'analyze'} == 1)&&((! $c_in{'trj_snaps'})||($c_in{'trj_snaps'} == 0))){
		print "\nERROR: LIE calculation was requested, but the number of snapshots\n";
		print "per trajectory (trj_snaps) was not specified. This information is required\n";
		print "for the preparation of the LIE calculations. Please provide the number of\n";
		print "snapshots per trajectory in the command file.\n"; 
	}
	elsif(($c_in{'analyze'} == 1) && ((! $c_in{'start'})||($c_in{'start'} == 0)
									||(! $c_in{'stop'})||($c_in{'stop'} == 0)
									||(! $c_in{'offset'})||($c_in{'offset'} == 0))){
		print "\nERROR: Energy calculation was requested, but not all parameters required\n";
		print "for snapshot extraction are specified in the command file.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&(($c_in{'mol2'} == 0)||($c_in{'sdf'} == 1))){
		print "\nERROR: Structures in mol2-format are expected as input structures. At this stage\n";
		print "of analysis structures in mol2-format should be available from the MD setup step.\n";
	}
	elsif(($c_in{'analyze'} == 1)&&($c_in{'pbs_lie'})&&($c_in{'pbs_lie'} ne "")&&(! -e $c_in{'pbs_lie'})){
		print "\nERROR: The batch file you specified\n";
		print $c_in{'pbs_lie'}." does not exist.\n";
		print "Please correct the name and path of the batch file.\n";
	}
	else{
		$check = 1;
	}
	
	# Check if enough structures are present in the specified trajectories
	if($check == 1){
		if((exists $c_in{'analyze'})&&($c_in{'analyze'} == 1)){
			my $total_snap_no = $c_in{'stop'} - $c_in{'start'} + 1;
			my $total_struct_no = $c_in{'trj_snaps'} * @traj_files;
			if($total_struct_no < $total_snap_no){
				print "\nERROR: The number of requested snapshots exceeds the number of structures present\n";
				print "in the given trajectories.\n";
				$check = 0;
			}

			if($check == 1){
				if((! $c_in{'lie_exe'})||($c_in{'lie_exe'} eq "")){
					print "\nWARNING: LIE executable specification missing in command file. It is assumed\n";
					print "that the program is located under ".$current_FEW_dir."/miscellaneous/LIE.pl\n";
				}
				my $total_snaps_no = (($c_in{'stop'} - $c_in{'start'}) / $c_in{'offset'});
				my $modulo_total = $total_snaps_no % $c_in{'offset'};
				if($modulo_total == 0){
					$total_snaps_no = $total_snaps_no + 1;
				}
				if($total_snaps_no < 3){
					print "\nERROR: LIE analysis can only be run for three snapshots or more to allow\n";
					print "a statistical analysis.\n";
					$check = 0;
				}
				if(($check == 1)&&($c_in{'stop'} > ($c_in{'trj_snaps'} * $traj_no))){
					print "\nERROR: LIE calculation was requested for a larger number of structures, than\n";
					print "are present in the specified trajectories. Please check the settings for\n";
					print "'snaps_per_trajectory', 'first_lie_snapshot', 'last_lie_snapshot', and\n";
					print "'offset_lie_snapshots' in the command file.\n";
					$check = 0;
				}
			}
		}
	}
	
	return $check;
}


# Specific check of input parameters for TIW module
sub check_input_tiw{

	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $check = 0;
	
	if(($c_in{'sim'} == 1)&&($c_in{'traj'} != 3)){
		print "\nERROR: You requested a ".$c_in{'traj'}."-trajectory setup in TI calculation modus.\n";
		print "Only 3-trajectory approach works for Thermodynamic integration calculations.\n";
	}
	elsif(($c_in{'sim'} == 1)&&((($c_in{'prod'})&&($c_in{'prod'} ne ""))||(($c_in{'prod_t'})&&($c_in{'prod_t'} ne ""))
	                         ||(($c_in{'pbs_p'})&&($c_in{'pbs_p'} ne "")))){
		print "\nWARNING: You requested setup of molecular dynamics minimization and equilibration\n";
		print "for later use in thermodynamic integration calculations. Please note that the parameters\n";
		print "you specified for MD production will not be regarded.\n";
		$check = 1;
	}
	
	elsif(($c_in{'ti_sim'} == 1)&&($c_in{'crd_top'} == 1)&&($c_in{'match_list'} ne "")&&(! -e $c_in{'match_list'})){
		print "\nERROR:\nThe matching list you provided could not be found. Please ensure that the correct\n";
		print "path and name of the matching file are given in the command file\n";
	}
	
	elsif(($c_in{'ti_sim'} == 1)&&(($c_in{'Lname_start'} eq "")||($c_in{'Lname_end'} eq "")
								||($c_in{'Lalias_start'} eq "")||($c_in{'Lalias_end'} eq ""))){
		print "\nERROR:\nNot all names and alias for start- and end-structures specified in command file.\n";
		print "Please make sure you specified names and alias for start- and end-structures correctly.\n";
	}

	elsif(($c_in{'ti_sim'} == 1)&&(($c_in{'chrg_meth'} ne "resp")&&($c_in{'chrg_meth'} ne "am1"))){
		print "\nERROR: Charge method must be either 'resp' or 'am1'\n";
		print "Please ensure you specified a correct charge method.\n";
	}

	elsif(($c_in{'ti_sim'} == 1)&&(! -e $c_in{'rst_lig'})){
		print "\nERROR: Specified RESTRT-file of ligand start-structure (state V0)\n";
		print $c_in{'rst_lig'}." does not exist.\n";
	}
	
	elsif(($c_in{'ti_sim'} == 1)&&(! -e $c_in{'rst_com'})){
		print "\nERROR: Specified RESTRT-file of complex start-structure (state V0)\n";
		print $c_in{'rst_com'}." does not exist.\n";
	}

	elsif(($c_in{'ti_sim'} == 1)&&(! -e $c_in{'top_lig'})){
		print "\nERROR: Specified topology of ligand start-structure (state V0)\n";
		print $c_in{'top_lig'}." does not exist.\n";
	}

	elsif(($c_in{'ti_sim'} == 1)&&(! -e $c_in{'top_com'})){
		print "\nERROR: Specified topology of complex start-structure (state V0)\n";
		print $c_in{'top_com'}." does not exist.\n";
	}
	
	elsif(($c_in{'ti_sim'} == 1)&&(($c_in{'mol2'} == 0)||($c_in{'sdf'} == 1))){
		print "\nERROR: Structures in mol2-format are expected as input structures. At this stage\n";
		print "of analysis structures in mol2-format should be available from the MD setup step.\n";
	}
	elsif(($c_in{'ti_sim'} == 1)&&((length $c_in{'Lalias_start'} > 3)||(length $c_in{'Lalias_end'} > 3))){
		print "\nERROR:\nOne of the aliasses you specified for state V0 or V1 comprises too many characters.\n";
		print "Please ensure that the alias comprises exactly 3 characters.\n";
	}
	elsif(($c_in{'ti_sim'} == 1)&&((length $c_in{'Lalias_start'} < 3)||(length $c_in{'Lalias_end'} < 3))){
		print "\nERROR:\nOne of the aliasses you specified for state V0 or V1 comprises not enough characters.\n";
		print "Please ensure that the alias comprises exactly 3 characters.\n";
	}
	else{
		$check = 1;
	}
	
	# Check mask pattern
	if($check == 1){
		if(! ($c_in{'mask0'} =~ m/\w@[\w\,]+\w/)){
			print "\nERROR: The 'softcore_mask_v0' information in the command file has not the required format.\n";
			$check = 0;
		}
		if(! ($c_in{'mask1'} =~ m/\w@[\w\,]+\w/)){
			print "\nERROR: The 'softcore_mask_v1' information in the command file has not the required format.\n";
			$check = 0;
		}
	}
	
	# Check existence of mol2 files
	if(($check == 1)&&(($c_in{'ti_sim'} == 1)&&($c_in{'sybyl_mol2'} == 1)&&($c_in{'Lname_start'})&&($c_in{'Lname_end'})&&($c_in{'i_path'}))){
		my $mol2_file_start = $c_in{'i_path'}."/".$c_in{'Lname_start'}.".mol2";
		my $mol2_file_end = $c_in{'i_path'}."/".$c_in{'Lname_start'}.".mol2";
		if(! -e $mol2_file_start){
			print "\nERROR:\nMol2-file for state V0 cannot be found in specified input structure directory ".$c_in{'i_path'}."\n";
			print "Please ensure you specified the correct input directory.\n";
			$check = 0;
		}
		if(! -e $mol2_file_end){
			print "\nERROR:\nMol2-file for state V1 cannot be found in specified input structure directory ".$c_in{'i_path'}."\n";
			print "Please ensure you specified the correct input directory.\n";
			$check = 0;
		}
	}
	
	# Check number of time values provided
	if(($check == 1)&&(($c_in{'ti_sim'} == 1)&&($c_in{'ti_prod'} == 1))){
		my @prod_time = split(/\,/, $c_in{'prod_time'});
		my @lambda = split(/\,/, $c_in{'equi_lambda'});
		my @prod_lambda = split(/\,/, $c_in{'prod_lambda'});
		
		if(@prod_time != @prod_lambda){
			print "\nERROR:\nNumber for time values given for MD-production not equal to\n";
			print "the number of lambda values. Please ensure that a production\n";
			print "time value is given for every lambda value.\n";
			$check = 0;
		}
		if(@prod_lambda > @lambda){
			print "\nERROR:\nInvalid number of lambda steps specified for production run.\n";
			print "The number of lambda steps used in the production must be smaller\n";
			print "or equal to the number of lambda steps for which the equilibration\n";
			print "was performed. Please correct your input script.\n";
			$check = 0;
		}
		
		foreach my $prod_l (@prod_lambda){
			if(grep {$prod_l eq $_} @lambda){
				next;
			}
			else{
				print "ERROR: For lambda step $prod_l no equilibration was conducted according to your\n";
				print "command file settings.\n";
				$check = 0;
				last;
			}
		}
	}  
	
	# Check format of lambda values
	if(($check == 1)&&(($c_in{'ti_sim'} == 1)&&($c_in{'ti_equi'}))){
		my @lambda = split(/\,/, $c_in{'equi_lambda'});
		my $ref_lambda = length($lambda[0]);
		for(my $l=1; $l<@lambda; $l++){
			my $length_l = length($lambda[$l]);
			if($length_l != $ref_lambda){
				print "\nERROR:\nProvided lambda values for TI equilibration have different number of digits.\n";
				print "Please ensure that all lambda values have the same number of digits, e.g.:\n";
				print "if requesting lambda=0.05, 0.1, 0.15, 0.2, ..., specify lambda=05,10,15,20,...\n";
				$check = 0;
			}
		}
	}
	
	# Check format of lambda values
	if(($check == 1)&&(($c_in{'ti_sim'} == 1)&&($c_in{'ti_prod'}))){
		my @lambda = split(/\,/, $c_in{'prod_lambda'});
		my $ref_lambda = length($lambda[0]);
		for(my $l=1; $l<@lambda; $l++){
			my $length_l = length($lambda[$l]);
			if($length_l != $ref_lambda){
				print "\nERROR:\nProvided lambda values for TI production have different number of digits.\n";
				print "Please ensure that all lambda values have the same number of digits, e.g.:\n";
				print "if requesting lambda=0.05, 0.1, 0.15, 0.2, ..., specify lambda=05,10,15,20,...\n";
				$check = 0;
			}
		}
	}
	
	return $check;	
}


# Checking existence of files required for preparation of free energy calculations
sub check_file_existence{
	my $procedure = shift;
	my $ref_c_in = shift;
	my $ref_struct_b = shift;
	my %c_in = %{$ref_c_in};
	my %struct_b = %{$ref_struct_b};
	my $check = 1;
		
	foreach my $s (keys %struct_b){
		if($procedure eq "liew"){
			# Ligand library file needed for topology setup
			my $ligand_lib = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}."_".$c_in{'chrg_meth'}."_mod.lib";
			my $alternative_lib = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}."_".$c_in{'chrg_meth'}.".lib";
			if(! -e $ligand_lib){
				if(! -e $alternative_lib){
					print "The library file needed for setup of the topology could not be found for structure $struct_b{$s}.\n";
					print "Please note that the same root path as for the setup of the MD simulations is required for\n";
					print "conducting the analysis.\n"; 
					$check = 0;
				} 
			}
		
			# Ligand parameter file
			my $frcmod = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}.".frcmod";
			if(! -e $frcmod){
				print "The parameter file required for setup of the topology files could not be found for structure $struct_b{$s}.\n";
				print "Please note that the same root path as for the setup of the MD simulations is required for\n";
				print "conducting the analysis.\n";
				$check = 0;
			} 
		
			for my $tag ("com", "lig"){
				# PDB structure used for setup of topology files
				my $pdb_template = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_".$tag.".pdb";
				if(! -e $pdb_template){
					print "PDB template file required for topology setup cannot be found for structure $struct_b{$s}.\n";
					print "Please ensure that the specified root-path is correct and the file $pdb_template exists.\n";
					$check = 0;
				}
			}
		}
	}
	return $check;
}


# Checking existence of trajectory files
sub check_traj_files{

	my $procedure = shift;
	my $ref_traj_files = shift;
	my $struct_b_ref = shift;
	my $ref_c_in = shift;
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$ref_c_in}; 
	my @traj_files = @{$ref_traj_files};
	my $check = 1;
	my $traj_file = "";
	
	foreach my $s (keys %struct_b){
		foreach my $traj (@traj_files){
			if($procedure eq "wamm"){
				if($c_in{'traj'} == 1){			
					#START LW_03_2018
					# allow custom path for protein-protein calculations #
					if($c_in{'prot_prot'} == 1){
						$traj_file = $traj;
						if(! -e $traj_file){
							$check = 0;
							last;
						}
					}
					else{
						$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/com/prod/".$traj;
						if(! -e $traj_file){
							$check = 0;
							last;
						}
					}
					#END LW_03_2018
				}
				if($c_in{'traj'} == 3){
					foreach my $tag ("lig", "com"){
						$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj;
						if(! -e $traj_file){
							$check = 0;
							last;
						}
					}
					if($check == 1){
						$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/rec/prod/".$traj;
						if(! -e $traj_file){
							$check = 0;
							last;
						}
					}
					else{
						last;
					}
				}
			}
			if($procedure eq "liew"){
				foreach my $tag ("lig", "com"){
					$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj;
					if(! -e $traj_file){
						$check = 0;
						last;
					}
				}
			}
			if($check == 0){
				last;
			}
		}
		if($check == 0){
			last;
		}
	}
	
	if($check == 0){
		print "\nERROR: At least one of the expected trajectory files does not exist.\n";
		print "Please check existence of file $traj_file.\n\n";
	}
	
	return $check;
}


# Checking existence of imaged trajectories in case imaging was requested.
sub check_imaged_traj_presence{

	my $procedure = shift;
	my $ref_traj_files = shift;
	my $struct_b_ref = shift;
	my $ref_c_in = shift;
	my %struct_b = %{$struct_b_ref};
	my %c_in = %{$ref_c_in}; 
	my @traj_files = @{$ref_traj_files};
	my $check = 1;
	my $traj_file = "";
	my $unzipped_traj_file = "";
	
	foreach my $s (keys %struct_b){
		foreach my $traj (@traj_files){
		
			my @traj_file_s = split(/\./, $traj);
		
			if(($procedure eq "wamm")&&($c_in{'traj'} == 1)){		
				$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/com/prod/".$traj_file_s[0]."_img.mdcrd.gz";
				$unzipped_traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/com/prod/".$traj_file_s[0]."_img.mdcrd";
				if((-e $traj_file)||(-e $unzipped_traj_file)){
					$check = 0;
					last;
				}
			}
			if((($procedure eq "wamm")&&($c_in{'traj'} == 3))||($procedure eq "liew")){
				foreach my $tag ("lig", "com"){
					$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj_file_s[0]."_img.mdcrd.gz";
					$unzipped_traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj_file_s[0]."_img.mdcrd";
					if((-e $traj_file)||(-e $unzipped_traj_file)){
						$check = 0;
						last;
					}
				}
				if(($check == 1)&&($procedure eq "wamm")){
					$traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/rec/prod/".$traj_file_s[0]."_img.mdcrd.gz";
					$unzipped_traj_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/rec/prod/".$traj_file_s[0]."_img.mdcrd";
					if((-e $traj_file)||(-e $unzipped_traj_file)){
						$check = 0;
						last;
					}
				}
				else{
					last;
				}
			}
			if($check == 0){
				last;
			}
		}
		if($check == 0){
			last;
		}
	}
	
	if($check == 0){
		my $detected_img_traj = "";
		if(-e $traj_file){
			$detected_img_traj = $traj_file;
		}
		if(-e $unzipped_traj_file){
			$detected_img_traj = $unzipped_traj_file;
		}
	
		print "\nWARNING:\n";
		print "You requested imaging of trajectories before snapshot extraction. For at least\n";
		print "one of the specified trajectories there exists already an imaged version of the\n";
		print "trajectory, see e.g. ".$detected_img_traj.".\n";
	}
}



# Subroutine for checking the provided receptor structure.
# Currently: Check for connectivity of backbone atoms.
sub check_rec_struct{
	my $rec_struct = shift;
	my $res_count = 0;
	my %back_atoms;
	my $prev_resno = 0;
	my $prev_resname = "";
	my $prev_rescount;
	my %ters;
	my $struct_check = 1;
	my %nme_block_groups;
	
	# Read backbone information from file
	open(REC, $rec_struct) || die "Cannot open file $rec_struct for reading.\n";
	while(my $pdb_l = <REC>){
		chomp($pdb_l);
		$pdb_l =~ s/^\s+//g;
		$pdb_l =~ s/\s+$//g;

		if($pdb_l =~ m/^ATOM/){

			my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
		unpack("a6    a5  x   a4     a     a3    x2  a4   x4   a8  a8  a8", $pdb_l);
			
			$atom =~ s/^\s+//g;
			$atom =~ s/\s+$//g;
			$resno =~ s/^\s+//g;
			$resno =~ s/\s+$//g;
			$resname =~ s/^\s+//g;
			$resname =~ s/\s+$//g;
		
			if((($prev_resno != 0)&&($prev_resno != $resno))||(($prev_resname ne "")&&($prev_resname ne $resname))){
				$res_count++;
				$back_atoms{$res_count}->{'org_res'} = $resno;
			}
			if(($prev_resno == 0)&&($prev_resname eq "")){
				$res_count++;
				$back_atoms{$res_count}->{'org_res'} = $resno;
			}
			
			if(($atom eq "N")||($atom eq "CA")||($atom eq "C")){
				$x =~ s/^\s+//g;
				$x =~ s/\s+$//g;
				$y =~ s/^\s+//g;
				$y =~ s/\s+$//g;
				$z =~ s/^\s+//g;
				$z =~ s/\s+$//g;
				$back_atoms{$res_count}->{$atom}->{'x'} = $x;
				$back_atoms{$res_count}->{$atom}->{'y'} = $y;
				$back_atoms{$res_count}->{$atom}->{'z'} = $z;
			}
			
			if($resname eq "NME"){
				$nme_block_groups{$res_count} = 1;
			}
			
			$prev_resno = $resno;
			$prev_resname = $resname;
			$prev_rescount = $res_count;
		}
		if($pdb_l =~ m/^TER/){
			$ters{$prev_rescount} = 1;
		}
	}
	close REC;
	
	# Check presence of TER card after ACE group
	foreach my $k (keys %nme_block_groups){
		if(! (exists $ters{$k})){
			print "ERROR:\n";
			print "Missing TER card in provided receptor PDB structure after residue\n";
			print "NME $k. Please ensure that all residues that are not directly bonded\n";
			print "to the next residue in the sequence are followed by a TER card in\n";
			print "the provided receptor PDB structure.\n";
			$struct_check = 0;
			last;
		}
	}
	
	# Check distances between residues - C of residue i and N of residue i+1
	for(my $i=1; $i<$res_count; $i++){
		if(exists $back_atoms{$i}->{'C'}){
			my $diff_x = ($back_atoms{$i}->{'C'}->{'x'} - $back_atoms{$i+1}->{'N'}->{'x'})**2;
			my $diff_y = ($back_atoms{$i}->{'C'}->{'y'} - $back_atoms{$i+1}->{'N'}->{'y'})**2;
			my $diff_z = ($back_atoms{$i}->{'C'}->{'z'} - $back_atoms{$i+1}->{'N'}->{'z'})**2;
			my $distance = sqrt($diff_x + $diff_y + $diff_z);
		
			if(($distance > 3)&&(! exists $ters{$i})){
				print "ERROR:\n";
				print "The distance between the C-alpha atom of residue ".$back_atoms{$i}->{'org_res'}." and the amide nitrogen\n";
				print "of residue ".$back_atoms{$i+1}->{'org_res'}." is larger than 3 Angstroem. Please check if there are missing\n";
				print "residues in your receptor structure or if the receptor consists of different chains\n";
				print "that are not separated by TER cards in the provided receptor PDB structure. In the\n";
				print "latter case, please ensure that all parts of the receptor that are not directly\n";
				print "connected are separated by TER cards in the receptor PDB structure.\n";
				$struct_check = 0;
				last;
			}
		}
	}
	return $struct_check;
}

# Subroutine for checking the format of the PDB file with membrane, ions, and
# water molecules. 
sub check_mem_file{
	my $mem_file = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $check_mem_file = 1;
	my @lipid_resnames = ("AR", "CHL", "DHA", "LA", "LEN", "LEO", "MY", "OL", "P2-",
					 	  "PA", "PC", "PE", "PGR", "PGS", "PH-", "PI", "PS", "ST",
						  "DMP", "PPE", "DLP", "DOP", "DPP", "DOP");
	
	# Read atom and residue information from file
	my $prev_resname = 0;
	my $prev_resno = 0;
	my $mem_residue_count = 0;
	my $lipids_found = 0;
	my $ions_found = 0;
	my $water_found = 0;
	
	open(MEM, $mem_file) || die "Cannot open file $mem_file for reading.\n";
	while(my $pdb_l = <MEM>){
		chomp($pdb_l);
		$pdb_l =~ s/^\s+//g;
		$pdb_l =~ s/\s+$//g;

		if(($pdb_l =~ m/^ATOM/)||($pdb_l =~ $pdb_l =~ m/^HETATM/)){

			my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
		unpack("a6    a5  x   a4     a     a3    x2  a4   x4   a8  a8  a8", $pdb_l);
			
			$atom =~ s/^\s+//g;
			$atom =~ s/\s+$//g;
			$resno =~ s/^\s+//g;
			$resno =~ s/\s+$//g;
			$resname =~ s/^\s+//g;
			$resname =~ s/\s+$//g;
			
			if((($resname ne $prev_resname)||($resno ne $prev_resno))&&(grep {$resname eq $_} @lipid_resnames)){
				$mem_residue_count++;
			}
			
			if(($mem_residue_count <= 3)&&(! (grep {$resname eq $_} @lipid_resnames))){
				print "\nERROR: The provided $mem_file does not contain any known type\n";
				print "of lipids from the GAFFlipid or the Lipid14 force field at the\n";
				print "beginning of the file. Please note that the PDB file with membrane\n";
				print "information has to have the following order: membrane lipids, (ions),\n";
				print "water molecules.\n";
				$check_mem_file = 0;
				last;
			}
			
			if($mem_residue_count > 3){
				$lipids_found = 1;
			}
			
			if(($lipids_found == 1)&&($water_found == 0)&&(($resname eq "Cl-")||($resname eq "Na+")||($resname eq "K+"))){
				$ions_found = 1;
			}
			
			if(($lipids_found == 1)&&(($resname eq "WAT")||($resname eq "HOH"))){
				$water_found = 1;
			}
			
			$prev_resname = $resname;
			$prev_resno = $resno;
		}
	}
	
	if(! (($lipids_found == 1)&&($water_found == 1))){
		print "\nERROR: The PDB file with membrane information specified in the command file\n";
		print "does not have the expected format. Please note that the PDB file with membrane\n";
		print "information has to have the following order: membrane lipids, (ions),\n";
		print "water molecules.\n";
		$check_mem_file = 0;
	}
	
	if(((exists $c_in{'rstMem'})&&($c_in{'rstMem'} != 0)&&($c_in{'rstMem'} != $mem_residue_count))
	 ||((exists $c_in{'memResNo'})&&($c_in{'memResNo'} != 0)&&($c_in{'memResNo'} != $mem_residue_count))){
		print "WARNING: It was specified in the command file that the membrane comprises\n";
		print $c_in{'rstMem'} ." residues, but $mem_residue_count number of residues were\n";
		print "found. Please ensure that the number of provided residues is correct.\n";
	}

	return $check_mem_file;
}


# Subroutine for checking if the number of residues in a receptor PDB structure
# is consistent with the residue number
sub check_resno{
	my $rec_file = shift;
	my $rec_resno = shift;
	my $prev_resno = 0;
	my $prev_resname = "";
	my $res_count = 0;
	my $ca_count = 1;
	
	open(REC, $rec_file) || die "Cannot open file $rec_file for reading.\n";
	while(my $pdb_l = <REC>){
		chomp($pdb_l);
		$pdb_l =~ s/^\s+//g;
		$pdb_l =~ s/\s+$//g;

		if(($pdb_l =~ m/^HETATM/)||($pdb_l =~ m/^ATOM/)){

			my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
		unpack("a6    a5  x   a4     a     a3    x2  a4   x4   a8  a8  a8", $pdb_l);
		
			if(($resname ne "HOH")&&($resname ne "WAT")){
				$resno =~ s/^\s+//g;
				$resno =~ s/\s+$//g;
				$resname =~ s/^\s+//g;
				$resname =~ s/\s+$//g;

				if((($prev_resno != 0)&&($prev_resno != $resno))||(($prev_resname ne "")&&($prev_resname ne $resname))){
					$res_count++;
				}
				if(($prev_resno == 0)&&($prev_resname eq "")){
					$res_count++;
				}
				
				$atom =~ s/^\s+//g;
				$atom =~ s/\s+$//g;
				if(($atom eq "CA")&&($resname ne "CA")){
					$ca_count++;
				}
				$prev_resno = $resno;
				$prev_resname = $resname;
			}
		}
	}
	close REC;
	
	my $check_passed = 0;
	if($res_count == $rec_resno){
		$check_passed = 1;
	}	

	return $check_passed, $ca_count, $res_count;
}


# Checking atom name length and residue name composition in provided mol2-files
sub check_mol2{
	my $s = shift;
	my $calling_module = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $start_reading = 0;
	my $check_mol = 1;
	my $consider_struct = 1;
	my @struct = split(/\//, $s);


	# Check if extension of ligand file is ".mol2"
	my $file_extension = substr($s, -4);
	if($file_extension ne "mol2"){
		print "\nERROR:\nThe file extension of the ligand file ".$s."\n";
		print "is not '.mol2'. Please ensure that the corresponding file is\n";
		print "a mol2-file and has the correct file extension.\n";
		exit;
	}

	# Check ligand structure name
	my $struct_name = substr($struct[$#struct], 0, -5);
	
	# Case A: First character is a digit
	my $first_char = substr($struct_name, 0, 1);
	if($first_char =~ m/(\d+)/){
		print "\nERROR:\nFirst character of ligand file ".$struct[$#struct]." is a digit.\n";
		print "Ligands that start with a number cannot be handled by FEW.\n";
		print "Please provide a ligand name that starts with a character.\n";
		exit;
	}
	
	# Case B: Residue name is equivalent to residue in AMBER force field
	my @ff_res_names = ();
	
	if($c_in{'backwards'} == 0){
	
		#ff14SB and ions13.lib
		@ff_res_names = ("A", "A3", "A5", "ACE", "ALA", "AN", "ARG", "ASH", 
						"ASN", "ASP", "Ag2+", "Ba2+", "Be2+", "Br-", "C", "C3",
						"C5", "CALA", "CARG", "CASN", "CASP", "CCYS", "CCYX",
						"CGLN", "CGLU", "CGLY", "CHCL3BOX", "CHID", "CHIE",
						"CHIP", "CHIS", "CHYP", "CILE", "CLEU", "CLYS", "CMET",
						"CN", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR",
						"CVAL", "CYM", "CYS", "CYX", "Ca2+", "Cd2+", "Cl-", "Co2+",
						"Cr2+", "Cs+", "Cu2+", "DA", "DA3", "DA5", "DAN", "DC",
						"DC3", "DC4", "DC5", "DCN", "DG", "DG3", "DG5", "DGN",
						"DT", "DT3", "DT5", "DTN", "Eu2+", "F-", "Fe2+", "G", "G3",
						"G5", "GLH", "GLN", "GLU", "GLY", "GN", "HID", "HIE",
						"HIP", "HIS", "HOH", "HYP", "Hg2+", "I-", "ILE", "K+",
						"LEU", "LYN", "LYS", "Li+", "MEOHBOX", "MET", "Mg+", "Mg2+",
						"Mn2+", "NALA", "NARG", "NASN", "NASP", "NCYS", "NCYX",
						"NGLN", "NGLU", "NGLY", "NHE", "NHID", "NHIE", "NHIP",
						"NHIS", "NILE", "NLEU", "NLYS", "NMABOX", "NME", "NMET",
						"NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
						"Na+", "Ni2+", "OHE", "PHE", "PL3", "POL3BOX", "PRO",
						"Pb2+", "Pd2+", "Pt2+", "QSPCFWBOX", "Ra2+", "Rb+", "SER",
						"SPC", "SPCBOX", "SPCFWBOX", "SPF", "SPG", "Sm2+", "Sn2+",
						"Sr2+", "T4E", "THR", "TIP3PBOX", "TIP3PFBOX", "TIP4PBOX",
						"TIP4PEWBOX", "TIP5PBOX", "TP3", "TP4", "TP5", "TPF", "TRP",
						"TYR", "U", "U3", "U5", "UN", "V2+", "VAL", "WAT", "Yb2+", "Zn2+");
	}
	
	# ff99SB
	if($c_in{'backwards'} == 1){
		@ff_res_names = ("ACE", "ALA", "ARG", "ASH", "ASN", "ASP", "CALA", "CARG",
						"CASN", "CASP", "CCYS", "CCYX", "CGLN", "CGLU", "CGLY", "CHCL3BOX",
						"CHID", "CHIE", "CHIP", "CHIS", "CILE", "CIO", "CLEU", "CLYS",
						"CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL",
						"CYM", "CYS",  "CYX", "Cl-", "Cs+", "DA", "DA3", "DA5",
						"DAN", "DC", "DC3", "DC4", "DC5", "DCN", "DG", "DG3",
						"DG5", "DGN", "DT", "DT3", "DT5", "DTN", "GLH", "GLN",
						"GLU", "GLY", "HID", "HIE", "HIP", "HIS", "HOH", "IB",
						"ILE", "K+", "LEU", "LYN", "LYS", "Li+", "MEOHBOX", "MET",
						"MG2", "NALA", "NARG", "NASN", "NASP", "NCYS", "NCYX", "NGLN",
						"NGLU", "NGLY", "NHE", "NHID", "NHIE", "NHIP", "NHIS", "NILE",
						"NLEU", "NLYS", "NMABOX", "NME", "NMET", "NPHE", "NPRO", "NSER",
						"NTHR", "NTRP", "NTYR", "NVAL", "Na+", "PHE", "PL3", "POL3BOX",
						"PRO", "QSPCFWBOX", "RA", "RA3", "RA5", "RAN", "RC", "RC3",
						"RC5", "RCN", "RG", "RG3", "RG5", "RGN", "RU", "RU3",
						"RU5", "RUN", "Rb+", "SER", "SPC", "SPCBOX", "SPCFWBOX", "SPF",
						"SPG", "T4E", "THR", "TIP3PBOX", "TIP3PFBOX", "TIP4PBOX", "TIP4PEWBOX",
						"TIP5PBOX", "TP3", "TP4", "TP5", "TPF", "TRP", "TYR", "VAL", "WAT");
	}
	

	if(grep {$struct_name eq $_} @ff_res_names){
		if($c_in{'backwards'} == 0){
			print "ERROR:\n";
			print "Name of ligand ".$struct_name." is equivalent to the name of a residue pre-defined in the\n";
			print "ff14SB force field or ions13.lib file of AMBER. Such a ligand name cannot be used for\n";
			print "calculation setup with FEW.\n";
			print "Please change the ligand name to a name not pre-defined in the ff14SB force field\n";
			print "of AMBER.\n";
		}
		else{
			print "ERROR:\n";
			print "Name of ligand ".$struct_name." is equivalent to the name of a residue pre-defined in the\n";
			print "ff99SB force field of AMBER. Such a ligand name cannot be used for calculation setup with FEW.\n";
			print "Please change the ligand name to a name not pre-defined in the ff99SB force field\n";
			print "of AMBER.\n";
		}
		exit;
	}

	if(($calling_module eq "TIW")&&($struct[$#struct] ne $c_in{'Lname_start'}.".mol2")&&($struct[$#struct] ne $c_in{'Lname_end'}.".mol2")){
		$consider_struct = 0;
	}

	if($consider_struct == 1){

		# Checking if antechamber is capable to handle provided mol2 file
		my $c_ante_check = "antechamber -i ".$s." -fi mol2 -o ".$s."_tmp -fo mol2 -at sybyl > log.txt 2>&1";
		call_ante_for_checking($c_ante_check, $s);
		unlink($s."_tmp");
		unlink glob "./ANTECHAMBER*";
		unlink("ATOMTYPE.INF");
		unlink("log.txt");

		# Checking of specific formatting features
		my $resname;
		my %atom_name;
		my $resid = 0;
		my $warn = 0;
		open(MOL, $s) || die "Cannot open file $s for reading.\n";
		while(my $mol = <MOL>){
			chomp($mol);
			if(($start_reading == 1)&&($mol eq "@<TRIPOS>BOND")){
				last;
			}

			if($start_reading == 1){

				$mol =~ s/^\s+//;
				my @mol = split(/\s+/, $mol);

				# Check whether the expected number of columns is present in the provided mol2-file
				if(@mol < 8){
					print "\nERROR:\n";
					print "The mol2-file ".$s."\n";
					print "does not contain the expected number of columns in the ATOM section.\n";
					print "Please ensure that all ligand mol2-files are correctly formatted.\n";
					print "Information about the recommended format for mol2-files can be found,\n";
					print "e.g. at http://www.tripos.com/data/support/mol2.pdf\n";
					print "In addition to the information that has to be provided in the ATOM\n";
					print "section of mol2-files, FEW requires the substructure ID and name,\n";
					print "i.e. the residue identifier and the residue name.\n";
					$check_mol = 0;
					last;
				}
				
				# Check whether special characters are present in ATOM section of mol2-file
				my $check_special_char = $mol;
				$check_special_char =~ s/[^a-zA-Z0-9\s\.\-\<\>]*//g;

				if($mol ne $check_special_char){
					print "\nERROR:\n";
					print "The mol2-file ".$s."\n";
					print "contains at least one special character in the ATOM section. It is expected\n";
					print "that this section contains only letters, digits, spcaces or full stops.\n";
					$check_mol = 0;
					last;
				}
				
				# Check if the same atom name was read before -> Ensure that atom names are unique
				if(exists $atom_name{$mol[1]}){
					print "\nERROR:\n";
					print "The atom names in the mol2-file ".$s."\n";
					print "are not unique. The atom name ".$mol[1]." is present at least twice.\n";
					print "Please ensure that the atom names in the ligand mol2-files are unique.\n";
				}
					
				if(length($mol[1]) > 3){
					print "\nERROR:\n";
					print "The mol2-file of the ligand ".$struct[$#struct]." contains at least one atom\n";
					print "which name comprises more than 3 characters. Only atom names\n";
					print "with up to three characters can be correctly handled by the\n";
					print "progrom. Please make sure that the atom names in the input\n";
					print "are not longer than three characters.\n";
					$check_mol = 0;
					last;
				}
				my $first_character = substr($mol[7], 0, 1);
				$first_character =~ m/(\d+)/;
				my $check_digit = $1;
				if($check_digit != 0){
					print "\nERROR:\n";
					print "The residue name in the mol2-file\n";
					print $s."\n"; 
					print "starts with a digit.\n";
					print "This format is not accepted by LEaP. Please specify a residue name\n";
					print "for the ligand that starts with a character.\n";
					$check_mol = 0;
					last; 
				}
				if($first_character eq "*"){
					print "\nERROR:\n";
					print "\nThe residue name in the mol2-file\n";
					print $s."\n";
					print "starts with a star.\n";
					print "This format is not allowed. Please specify a residue name for the ligand\n";
					print "that consists only of characters and digits and starts with a character.\n";
					$check_mol = 0;
					last;
				}
				
				# Check if number of residues in receptor structure is consistent with provided structure
				if(($resid != 0)&&($mol[6] =~ m/(\d+)/)&&($mol[6] != $resid)){
					print "\nERROR: The ligand file $s\n";
					print "contains more than one residue. FEW can only handle ligands that consist of\n";
					print "one residue.\n";
					$check_mol = 0;
					last;
				}			
				if(($resname ne "")&&($mol[7] ne $resname)){
					print "\nERROR: More than one residue name was found in the ligand file\n";
					print "$s\n";
					print "FEW can only handle ligands that consist of one residue.\n";
					$check_mol = 0;
					last;
				}
				
				if(($warn == 0)&&(length($mol[7]) > 3)){
					print "\nWARNING:\n";
					print "The residue name in the ligand file $s\n";
					print "you provided is longer than three characters and will be shortened\n"; 
					print "to the first three characters.\n";
					$warn = 1;
				}
				if(($warn == 0)&&($mol[8] != 0)){
					print "\nWARNING: The ligand mol2 file $s contains charges. If you conduct a charge\n";
					print "calculation this can result in additional warning messages. To avoid this,\n";
					print "please use only ligand mol2 files with zero atomic charges as input.\n";
					$warn = 1;
				}	
				if($mol[6] =~ m/(\d+)/){
					$resid = $mol[6];
				}
				$resname = $mol[7];
				$atom_name{$mol[1]} = 1;
			} 
		
			if($mol =~ /@<TRIPOS>ATOM/){
				$start_reading = 1;
			}
		}
	
		# Check if atom section was found at all
		if($start_reading != 1){
			print "ERROR:\n";
			print "The provided ligand file ".$s."\n";
			print "does not seem to be a mol2-formmatted file. Only mol2-formatted files\n";
			print "are allowed as ligand input files at this stage of the workflow.\n";
			$check_mol = 0;
		}
	}
	return $check_mol;	
}


# Subroutine for antechamber call for checking usability of provided mol2-files
sub call_ante_for_checking{
	my $command = shift;
	my $struct_file = shift;

	system($command) == 0
				or ((unlink glob "./ANTECHAMBER*") && (unlink "ATOMTYPE.INF") && (die "ERROR: The provided mol2-file '".$struct_file."'
cannot be processed by antechamber. Please ensure, that all provided mol2-files\nare in a format that can be handled by antechamber.\n"));
	if($? != 0){
		if(! -z "log.txt"){
			my $log_str = "";
			open(LOG, "log.txt") || die "Cannot open file log.txt for reading.\n";
			while(my $l_log = <LOG>){
				$log_str .= $l_log;
			}
			close LOG;
			print $log_str."\n";
		} 
		unlink("log.txt");
		exit;
	}
}
 

sub check_dG_params{
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my $check = 1;
		
	if((! $c_in{'root_path'})||($c_in{'root_path'} eq "")||(! -e $c_in{'root_path'})){
		print "\nERROR: No root_path information could not be found in the command file or the provided root_path does not exist.\n";
		$check = 0;
	}
	elsif((! $c_in{'chrg_meth'})||(($c_in{'chrg_meth'} ne "am1")&&($c_in{'chrg_meth'} ne "resp"))){
		print "\nERROR: No charge method information could be found in the command file or the provided\n";
		print "charge method information is not correct.\n";
		$check = 0;
	} 
	
	elsif((! $c_in{'Lname_start'})||($c_in{'Lname_start'} eq "")||(! $c_in{'Lname_end'})||($c_in{'Lname_end'} eq "")){
		print "\nERROR: The names of start- and/or end-structure was not provided in the command file.\n";
		print "Please ensure that names of start- and end-structure are provided in the command file.\n";
		$check = 0;
	}
	
	elsif((! $c_in{'Lalias_start'})||($c_in{'Lalias_start'} eq "")||(! $c_in{'Lalias_end'})||($c_in{'Lalias_end'} eq "")){
		print "\nERROR: The aliases of start- and/or end-structure was not provided in the command file.\n";
		print "Please ensure that aliases of start- and end-structure are provided in the command file.\n";
		$check = 0;
	}
	
	# Check existence of required folders / files	
	if($check == 1){
		my $transform_dir = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'};
		my $prod_dir = $transform_dir."/prod";
		
		if(! -e $transform_dir){
			print "\nERROR: The folder of the requested transformation to be analyzed cannot be found in the expected location at:\n";
			print $transform_dir."\n";
			$check = 0;
		}
		elsif(! -e $prod_dir){
			print "\nERROR: The folder containing the production simulations to be analyzed cannot be found in the expected location at:\n";
			print $prod_dir."\n";
			$check = 0;
		}
		foreach my $tag ("com", "lig"){
			my @sander_out_files = <$prod_dir/$tag/*.out>;
			if(@sander_out_files == 0){
				print "\nERROR: Output files requested to be analyzed cannot be found in the expected location at:\n";
				print $prod_dir."/".$tag."/\n";
				$check = 0;
			}
		}
	}
	
	if((! $c_in{'calc_meth'})||($c_in{'calc_meth'} == 0)){
		print "ERROR: You requested ddG computation, but did not specify a calculation method. Please set the\n";
		print "calculation method either to 1 or 2 in the command file.\n";
		$check = 0;
	}
	
	return $check;
}

1;
