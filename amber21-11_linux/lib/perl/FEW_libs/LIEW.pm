#/usr/bin/perl -w

# Perl module for the preparation of structures for MD input and for the determination
# of binding energies with the linerar interaction energy approach.

use File::Copy;
use File::Path;
use strict;

use global;
use check_input;
use separate_structures;
use common_prepare_MD;
use read_data;

sub LIEW{

	my $current_dir = shift;
	my $c_file = shift;
	my $traj_files_ref = shift;
	my $c_in_ref = shift;
	my @traj_files = @{$traj_files_ref};
	my %c_in = %{$c_in_ref};

	##################################
	# Specific check for LIEW module #
	##################################
	
	if($c_in{'sim'} == 1){
		my $check_liew = &check_input_liew($current_dir, \%c_in, \@traj_files);
		if($check_liew != 1){
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
		my $check_mol = check_mol2($s, "LIEW", \%c_in);
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

			create_crd_top("liew", $chrg_meth, $md_dir, $current_dir, \@struct_files, \%struct_b, \%c_in);

			####################################
			#  Create scripts for equilbration #
			####################################

			# Setup folders, input scripts and pbs script for equilibration
			if(($c_in{'equi'})&&(-e $c_in{'equi'})){
				$warnings = prepare_equilibration("liew", $chrg_meth, $md_dir, \@struct_files, \%struct_b, \%c_in);
			}


			#####################################
			#  Create scripts for MD-production #
			#####################################

			if(($c_in{'prod'})&&(-e $c_in{'prod'})&&($c_in{'prod_t'})){
				prepare_MD_production("liew", $chrg_meth, $md_dir, $warnings, \@struct_files, \%struct_b, \%c_in);
			}
		}
	}
	
	
	############################
	#  Setup LIE calculations  #
	############################
	if($c_in{'analyze'} == 1){

		# Set trajectory to 3 trajectory approach, since this is the only approach that
		# works for LIE analysis, receptor data are neglected.
		$c_in{'traj'} = 3;

		# Read trajectories into array 'traj_files' if 'all' was specified instead of trajectories.
		if(($traj_files[0] eq "all")||($c_in{'trj_file'} eq "all")){
			@traj_files = ();
			my $traj_files_ref = read_trajectory_files($struct_b{$struct_files[0]}, \%c_in);
			@traj_files = @{$traj_files_ref};
		}
		
		my $check_liew = &check_input_liew($current_dir, \%c_in, \@traj_files);
		if($check_liew != 1){
			exit;
		}
		
		# Check existence of files required for analysis
		my $check_files = check_file_existence("liew", \%c_in, \%struct_b);
		if($check_files != 1){
			exit;
		}

		# Check if specified trajectories needed for analysis exist
		my $check_traj = check_traj_files("liew", \@traj_files, \%struct_b, \%c_in);
		if($check_traj != 1){
			exit;
		}
		
		# When imaging is requested, check if imaged trajectories are already present
		if($c_in{'image'} == 1){
			check_imaged_traj_presence("liew", \@traj_files, \%struct_b, \%c_in);
		}
		
		# Check provided residue number of receptor structure using receptor structure prepared for first ligand
		my $rec_file = $c_in{'root_path'}."/leap/".$struct_b{$struct_files[0]}."/".$struct_b{$struct_files[0]}."_rec.pdb";
		if(! -e $rec_file){
			print "WARNING: Could not find receptor file $rec_file\n";
			print "for checking consistency of receptor residue information.\n";
		}
		else{
			my ($check_rec_resno, $ca_no, $resno) = check_resno($rec_file, $c_in{'rec_res'});
			if($check_rec_resno == 0){
				print "ERROR: $ca_no C-alpha atoms and $resno residues were detected in the\n";
				print "receptor structure ".$rec_file.".\n";
				print "The number of residues differs from the one specified in the command-file.\n"; 
				exit;
			}
		}
		
		# Generate directory for storing analysis results
		my $analyze_dir = $c_in{'root_path'}."/lie_".$c_in{'chrg_meth'};

		if(! -e $analyze_dir){
			mkdir($analyze_dir);
		}

		# Prepare submission script if batch setup is requested
		my $qsub_file = $analyze_dir."/qsub_LIE.sh";
		if(($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
			open(QSUB, ">$qsub_file") || die "Cannot open file $qsub_file for writing.\n";
		} 

		foreach my $s (@struct_files){

			print "\nGeneration of files for LIE calculation for structure ".$struct_b{$s}."\n";

			if(! -e "$analyze_dir/".$struct_b{$s}){
				mkdir("$analyze_dir/".$struct_b{$s});
			}

			my @comp = ();

			for my $tag ("com", "lig"){

				# Generate complex/ligand directory
				my $tag_dir = $analyze_dir."/".$struct_b{$s}."/".$tag;
				if(! -e $tag_dir){
					mkdir($tag_dir);
				}			
				
				# Define components to regard
				if($tag eq "lig"){
					@comp = ("tot", "lig", "wat");
				}
				if($tag eq "com"){
					if((defined $c_in{'calc_sasa'})&&($c_in{'calc_sasa'} == 1)){
						@comp = ("tot", "com", "res", "lig");
					}
					else{
						@comp = ("tot", "res", "lig");
					}
				}
				
				
				# Topology creation
				 
				# Generation of topology directory
				my $topo_dir = $analyze_dir."/".$struct_b{$s}."/".$tag."/topo";
				
				if(! -e $topo_dir){
					mkdir($topo_dir);
				}
  			
  				# Identify PDB-template file
				my $pdb_template = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_".$tag.".pdb";
	
				foreach my $comp (@comp){
    	
					# Generate pdb-files for topology setup
					my $pdb_out_str = "";
					my $remove_ter = 0;
					
					open(PDB_T, $pdb_template) || die "Cannot open file ".$pdb_template."\n"; 
					
					while(my $pdb_t = <PDB_T>){
						chomp($pdb_t);
						my @pdb_t = split(/\s+/, $pdb_t);
						
						# Remove superfluous parts of PDB structure 
						if(($pdb_t[0] eq "TER")&&($remove_ter == 1)){
							$remove_ter = 0;
							next;
						}
						if($comp eq "lig"){
							if(($pdb_t[3] ne "MOL")&&($pdb_t[3] ne "<1>")){
								$remove_ter = 1;
								next;
							}
						}
						if(($tag eq "com")&&($comp eq "res")){
							if(($pdb_t[3] eq "MOL")||($pdb_t[3] eq "<1>")){
								$remove_ter = 1;
								next;
							}
						}
						if(($tag eq "com")&&($comp eq "com")){
							if($pdb_t[4] > ($c_in{'rec_res'} + 1)){
								$remove_ter = 1;
								next;
							}
						}
						if((($tag eq "com")&&($comp eq "com"))||(($tag eq "com")&&($comp eq "res"))||(($tag eq "com")&&($comp eq "tot"))){
							if($pdb_t =~ m/ATOM/){
				                        	my($card,  $atnum, $atom,     $resname,    $resno,      $x, $y, $z) =
								unpack("a6     a5 x    a4     x   a3        x2 a4      x4   a8  a8  a8", $pdb_t);
								$resname =~ s/^\s+//g;
								$resname =~ s/\s+$//g;

								# Don't use ACE/NME atoms not recognized by leap
								if((($resname eq "ACE")&&($atom =~ m/HH3/))||(($resname eq "NME")&&($atom =~ m/^H/))){
									next;
								}
							}
						}
						if(($tag eq "lig")&&($comp eq "wat")){
							if(($pdb_t[3] eq "MOL")||($pdb_t[3] eq "<1>")){
								$remove_ter = 1;
								next;
							}
						}
						if($pdb_t[3] eq "MOL"){
							$pdb_t =~ s/MOL/<1>/;
						}
						$pdb_out_str .= $pdb_t."\n";
					}
					close PDB_T;
  
					open(PDB_FOR_TOPO, ">$topo_dir/$comp"."IN.pdb") || die "Cannot open file $topo_dir/$comp"."IN.pdb for writing.\n";
					print PDB_FOR_TOPO $pdb_out_str;
					close PDB_FOR_TOPO;
  
					# Create LEAP file for toplogy generation
					open(LEAP_S, ">$topo_dir/leap_".$comp.".in") || "Cannot open file $topo_dir/leap_".$comp.".in";
					if($c_in{'backwards'} == 0){
						print LEAP_S "source leaprc.protein.ff14SB\n";
                                                if($c_in{'water_model'} eq "TIP3P"){
                                                	print LEAP_S "source leaprc.water.tip3p\n";
						}
						if($c_in{'water_model'} eq "OPC"){
							print LEAP_S "source leaprc.water.opc\n";
						}
					}
					else{
						print LEAP_S "source ".$c_in{'amberhome'}."/dat/leap/cmd/oldff/leaprc.ff99SB\n";
                                                # need to source tip3p because we are removing tip3p parameters from the protein ff
	                                        print LEAP_S "loadAmberParams ".$c_in{'amberhome'}."/dat/leap/parm/frcmod.tip3p\n";

					}
					print LEAP_S "source leaprc.gaff\n";
					print LEAP_S "set default PBRadii bondi\n";
					if(($comp eq "lig")||($comp eq "tot")||($comp eq "com")){
						my $ligand_lib = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}."_".$c_in{'chrg_meth'}."_mod.lib";
						if(! -e $ligand_lib){
							$ligand_lib = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}."_".$c_in{'chrg_meth'}.".lib";
						}
						print LEAP_S "loadoff ".$ligand_lib."\n";
					}
					if(($c_in{'add_lib'} ne "")&&(($comp eq "tot")||($comp eq "res")||($comp eq "com"))){
						print LEAP_S "loadoff ".$c_in{'add_lib'}."\n";
					}
					if(($c_in{'add_frcmod'} ne "")&&(($comp eq "tot")||($comp eq "res")||($comp eq "com"))){
						print LEAP_S "loadAmberParams ".$c_in{'add_frcmod'}."\n";
					}
					print LEAP_S "SYS = loadpdb $topo_dir/$comp"."IN.pdb\n";
					if(($comp eq "lig")||($comp eq "tot")||($comp eq "com")){
						my $frcmod_file = $c_in{'root_path'}."/leap/".$struct_b{$s}."/".$struct_b{$s}.".frcmod";
						print LEAP_S "loadAmberParams ".$frcmod_file."\n";
					}
					if(($tag eq "com")&&(($comp eq "tot")||($comp eq "res")||($comp eq "com"))){
						# Check whether disulfide bonds were defined in MD setup
						my $old_leap_script = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/leap_script_1.in";
						my $cys_str = "";
						open(OLD_LEAP, $old_leap_script) || die "Cannot open file $old_leap_script for reading.\n";
						while(my $old_leap_line = <OLD_LEAP>){
							chomp($old_leap_line);
							if($old_leap_line =~ m/^bond REC/){
								$old_leap_line =~ s/REC/SYS/g;
								$cys_str .= $old_leap_line."\n";
							}
						}
						close OLD_LEAP;
						
						if($cys_str ne ""){
							print LEAP_S $cys_str;
						}
					}
					print LEAP_S "saveAmberParm SYS $topo_dir/$comp.top $topo_dir/$comp.crd\n";
					print LEAP_S "savepdb SYS $topo_dir/$comp.pdb\n";
					print LEAP_S "quit\n";
					close LEAP_S;
  
					# Create Topology
					my $topo_log = $topo_dir."/leap_".$comp.".log";
					my $c_topo = "tleap -f ".$topo_dir."/leap_".$comp.".in > ".$topo_log." 2>&1";
					call_prog($c_topo, $topo_log);
					unlink("leap.log");
				}

			
				# Perform imaging if requested
				if($c_in{'image'} == 1){
					foreach my $tag ("com", "lig"){
			
						# Set path to trajectory files
						my $traj_path = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod";

						# Determine path to topology files
						my $top = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_".$tag.".top";
			
						# Determine residues to image
						my $img_res;
						my $tag_name;
						if($tag eq "com"){
							$img_res = $c_in{'rec_res'}+1;
							$tag_name = "complex";
						}
						else{
							$img_res = 1;
							$tag_name = "ligand";
						}
			
						for(my $f=1; $f<=@traj_files; $f++){
							my @traj_file_s = split(/\./, $traj_files[$f-1]);
							my $imaged_file = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj_file_s[0]."_img.mdcrd.gz";
							if((! -e $imaged_file)&&(! -l $imaged_file)){
								print "Imaging trajectory $f of $tag_name.\n";
								my ($image_c, $image_log, $cpptraj_file, $imaged_trj) = image_traj($traj_files[$f-1], $f, $img_res, $traj_path, $top, $c_in{'image_backwards'});
								call_prog($image_c, $image_log);
								my $gzip_log = $c_in{'root_path'}."/gzip.log";
 								my $gzip_c = "gzip ".$imaged_trj." 2>".$gzip_log;
 								call_prog($gzip_c, $gzip_log);
								unlink($gzip_log);
							}
						}
					}
				}

			
				# Generate trajectories for energy calculation
				################################################
				
				# Identify ions present in system
				my $ions;
				my $pdb_template = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/cryst/".$struct_b{$s}."_solv_".$tag.".pdb";
				my $res_reached = 0; # Avoiding stripping of structural bound ions that are part of the receptor
					
				open(PDB_T, $pdb_template) || die "Cannot open file ".$pdb_template."\n"; 
					
				while(my $pdb_t = <PDB_T>){
					chomp($pdb_t);
					my @pdb_t = split(/\s+/, $pdb_t);
					
					# Determine which sort of ions is present
					if(($res_reached == 1)&&($pdb_t[3] eq "Cl-")){
						$ions = "Cl-";
						last;
					}
					if(($res_reached == 1)&&($pdb_t[3] eq "Na+")){
						$ions = "Na+";
						last;
					}
					if(($tag eq "com")&&($pdb_t[4] == $c_in{'rec_res'})){
						$res_reached = 1;
					}
					if(($tag eq "lig")&&($pdb_t[4] == 2)){
						$res_reached = 1;
					}
				}
				close PDB_T;
				
				# Consider each component separately		
				foreach my $comp (@comp){
					
					# Generate directory for specific trajectories without box-information
					my $snaps_dir = $tag_dir."/s_".$comp; 
					
					if(! -e $snaps_dir){
						mkdir($snaps_dir);
					}
						
					# Residue number of ligand
					my $lig_res_no = $c_in{'rec_res'} + 1;
					my $traj_count = 0;
				
					# Generate trajectory without box information
					for(my $f=1; $f<=@traj_files; $f++){
					
						my $first_snap = $traj_count * $c_in{'trj_snaps'} + 1;
						$traj_count++;
						my $last_snap = $traj_count * $c_in{'trj_snaps'};
				
						# When start snapshot in trajectory
						if(($c_in{'start'} <= $last_snap)&&($c_in{'stop'} >= $first_snap)){
					
							# Determine trajectory that shall be used for imaging
							my $traj_in;
							my @traj_file_s = split(/\./, $traj_files[$f-1]);
							my $traj_out = $snaps_dir."/".$traj_file_s[0]."_nobox.mdcrd";
							if($c_in{'image'} == 1){
								$traj_in = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj_file_s[0]."_img.mdcrd.gz";
							}
							else{
								$traj_in = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/".$tag."/prod/".$traj_files[$f-1];
							}
													
							# Generate cpptraj-script for creation of trajectory that
							# can be post-processed for energy determination
							my $actual_start;
							my $actual_end;
							
							# Determine actual start snapshot
							my $current_first_snap = (($traj_count - 1) * $c_in{'trj_snaps'}) + 1;
							if($c_in{'start'} > $current_first_snap){
								$actual_start = $c_in{'start'} - (($traj_count - 1) * $c_in{'trj_snaps'});
							}
							else{
								my $modulo_snaps_no = (($traj_count - 1) * $c_in{'trj_snaps'}) % $c_in{'offset'};
								if($modulo_snaps_no == 0){
									$actual_start = 1;
								}
								if($modulo_snaps_no > 0){
									$actual_start = ($c_in{'offset'} - $modulo_snaps_no);
								}
							}
							
							# Determine actual end snapshot				
							if($c_in{'stop'} < $last_snap){
								$actual_end = $c_in{'stop'} - (($traj_count - 1) * $c_in{'trj_snaps'});
							}
							else{
								$actual_end = $c_in{'trj_snaps'};
							}
							
							# Create cpptraj-script
							my $cpptraj_script_nobox_in = $snaps_dir."/remove_boxinfo_$f.cpptraj";
							open(CPPTRAJ_NOBOX_IN, ">$cpptraj_script_nobox_in") || die "Cannot open file $cpptraj_script_nobox_in for writing.\n";
							print CPPTRAJ_NOBOX_IN "trajin ".$traj_in." ".$actual_start." ".$actual_end." ".$c_in{'offset'}."\n";
							if($comp eq "lig"){
								if($tag eq "com"){
									print CPPTRAJ_NOBOX_IN "strip :1-".$c_in{'rec_res'}."\n";
								}
								print CPPTRAJ_NOBOX_IN "strip :WAT\n";
								if($ions ne ""){
									print CPPTRAJ_NOBOX_IN "strip :$ions\n";
								}
							}
							if($comp eq "com"){
								print CPPTRAJ_NOBOX_IN "strip :WAT\n";
								if($ions ne ""){
									print CPPTRAJ_NOBOX_IN "strip :$ions\n";
								}
							}	
							if($comp eq "res"){
								print CPPTRAJ_NOBOX_IN "strip :MOL\n";
								print CPPTRAJ_NOBOX_IN "strip :".$lig_res_no."\n";
							}
							if($comp eq "wat"){
								print CPPTRAJ_NOBOX_IN "strip :MOL\n";
								print CPPTRAJ_NOBOX_IN "strip :1\n";
							}
							print CPPTRAJ_NOBOX_IN "trajout ".$traj_out." nobox\n";
							close CPPTRAJ_NOBOX_IN;
					
							# Extract snapshots / trajectories
							my $cpptraj_log_nobox = $snaps_dir."/remove_boxinfo_$f.log";
							my $c_snaps_tot = "cpptraj -p ".$topo_dir."/tot.top -i ".$cpptraj_script_nobox_in." > ".$cpptraj_log_nobox;
							call_prog($c_snaps_tot, $cpptraj_log_nobox);
						}
					} # End trajectory
				} # End comp
			} # End tag
			
			
			# Determine root path for calculation 
			my $pbs_path;
			if(($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
				if(($c_in{'pbs_l_path'})&&($c_in{'pbs_l_path'} ne "")){
					$pbs_path = $c_in{'pbs_l_path'};
				}
				else{
					$pbs_path = $c_in{'root_path'};
				}
			}
			
			
			# Prepare input script for energy calculation
			my $ener_in;
			my $root_path = $c_in{'root_path'};
			
			$ener_in = $analyze_dir."/".$struct_b{$s}."/energy_calc.in";
			my $ener_str = "";
				
			open(C_FILE, $c_file) || die "Cannot open $c_file for reading.\n";
			while(my $c_line = <C_FILE>){
				chomp $c_line;
				$c_line =~ s/\r//g;
				my $first_char = substr($c_line, 0, 1);

				if(($c_line =~ m/trajectory_files/)&&($c_line =~ m/ all/)){
					$c_line = "";
					for(my $f=1; $f<=@traj_files; $f++){
						$c_line .= "trajectory_files       $traj_files[$f-1]\n";
					}
					chomp($c_line);
				}
				if(($c_line =~ m/$root_path/)&&($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
					$c_line =~ s/${root_path}/${pbs_path}/;
				}
				if(((! exists $c_in{'lie_exe'})||($c_in{'lie_exe'} eq ""))&&($c_line =~ m/lie_exe/)&&($first_char ne "#")){
					$c_line = "lie_executable          $current_dir/miscellaneous/LIE.pl";
				}
				$ener_str .= $c_line."\n";
			}
			close C_FILE;
				
			open(ENER_IN, ">$ener_in") || die "Cannot open file $ener_in for writing.\n";
			print ENER_IN $ener_str;
			close ENER_IN;
			
			# Create pbs-script
			if(($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
				if((!$c_in{'lie_exe'})||($c_in{'lie_exe'} eq "")){
					$c_in{'lie_exe'} = "$current_dir/miscellaneous/LIE.pl";
				}

				my $pbs_str = "";
				open(PBS_TEMPL, $c_in{'pbs_lie'}) || die "Cannot open file ".$c_in{'pbs_lie'}." for reading.\n";
				while(my $pbs_l = <PBS_TEMPL>){
					chomp($pbs_l);
					if($pbs_l =~ m/ -N/){
						$pbs_l .= " LIE_".$struct_b{$s};
					}  
					if($pbs_l =~ m/ -e/){
						$pbs_l .= " ".$struct_b{$s}.".ER";
					}
					if($pbs_l =~ m/ -o/){
						$pbs_l .= " ".$struct_b{$s}.".OU";
					}
					#START - For SLUM script support
					if($pbs_l =~ m/ --job-name=/){
						$pbs_l .= "LIE_".$struct_b{$s};
					}  
					if($pbs_l =~ m/ --error=/){
						$pbs_l .= $struct_b{$s}.".ER";
					}
					if($pbs_l =~ m/ --output=/){
						$pbs_l .= $struct_b{$s}.".OU";
					}
					#END - For SLUM script support
					if(($pbs_l =~ m/ PATH=/)||($pbs_l =~ m/ WORKDIR=/)){
						$pbs_l .= $pbs_path;
					}
					$pbs_str .= $pbs_l."\n";	
				}
				close PBS_TEMPL;
				
				$pbs_str .= "perl ".$c_in{'lie_exe'}." ".$pbs_path."/lie_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/energy_calc.in ".$struct_b{$s}."\n";
					
				if($c_in{'del_trajs'} == 1){
					$pbs_str .= "rm -rf ".$pbs_path."/".$c_in{'chrg_meth'}."/".$struct_b{$s}."/com/s_*\n";
					$pbs_str .= "rm -rf ".$pbs_path."/".$c_in{'chrg_meth'}."/".$struct_b{$s}."/lig/s_*\n";
				}
					
				my $pbs_file = $analyze_dir."/".$struct_b{$s}."/run_".$struct_b{$s}.".pbs";
				open(PBS_F, ">$pbs_file") || die "Cannot open file $pbs_file for writing.\n";
				print PBS_F $pbs_str;
				close PBS_F;
				chmod 0755, $pbs_file;	
			} # End: Create batch-script
			
			if(($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
				print QSUB "qsub ".$pbs_path."/lie_".$c_in{'chrg_meth'}."/".$struct_b{$s}."/run_".$struct_b{$s}.".pbs\n";
			}
		} # End: Process structures
		if(($c_in{'pbs_lie'})&&(-e $c_in{'pbs_lie'})){
			close QSUB;
			chmod 0755, $qsub_file;
		}
	} # End: Analysis
} # End: Subroutine

1;
