use strict;

# Generate input scripts for equilibration
sub generate_scripts_equi{
	my $tag = shift;
	my $equi_dir = shift;
	my $equi_files_ref = shift;
	my @equi_files = @{$equi_files_ref};
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};

	foreach my $equi_f (@equi_files){
		
		# Determine name of file to be copied
		my @file = split(/\//, $equi_f);
		my $target_file = $equi_dir."/".$file[$#file];

		copy($equi_f, $target_file) || die "Cannot copy $equi_f to $target_file!\n";
		if($c_in{'rec_res'}){
			assign_res_restrict($target_file, $tag, \%c_in);
		}
	}
}


# Define restrained residues in minimization and equilibration scripts
sub assign_res_restrict{
	my $equi_file = shift;
	my $tag = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $res;
	my $content;

	if($tag eq "com"){
		$res = $c_in{'rec_res'} + 1;
		
		# Add restraints for membrane, if membrane simulation setup was requested
		# Currently membrane setup is only available for 1-trajectory approach
		if(((exists $c_in{'withMem'})&&($c_in{'withMem'} == 1))&&(exists $c_in{'rstMem'} != 0)){
			$res = $res + $c_in{'rstMem'};
		}
	}
	if($tag eq "rec"){
		$res = $c_in{'rec_res'};
	}
	if($tag eq "lig"){
		$res = 1;
	}

	# Read old file
	open(E_FILE, $equi_file) || die "Cannot open $equi_file for writing.\n";

	while(my $e_line = <E_FILE>){
		chomp($e_line);
		$e_line =~ s/\r//;
		if($e_line =~ m/^RES/){
			$content .= "RES 1 ".$res."\n";
		}
		else{
			$content .= $e_line."\n";
		}
	}
	close EFILE;
	unlink($equi_file);

	# Write file with new residue restraint
	open(E_NEW, ">$equi_file") || die "Cannot open $equi_file for writing.\n";
	print E_NEW $content;
	close E_NEW;
}


# Generate PBS-Script for equilibration
sub generate_pbs_equi{

	my $equi_dir = shift;
	my $tag = shift;
	my $struct_b = shift;
	my $chrg_meth = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};

	if(($c_in{'pbs_e'})&&(-e $c_in{'pbs_e'})){
		open(PBS_E, $c_in{'pbs_e'}) || die "Cannot open PBS-Template $c_in{'pbs_e'}.\n";
		open(PBS_OUT, ">$equi_dir/$struct_b.pbs") || die "Cannot open $equi_dir/$struct_b.pbs for writing.\n";

		while(my $pbs_line = <PBS_E>){
			chomp($pbs_line);

			if($pbs_line =~ m/set SCRIPT/){
				if($tag eq "rec"){
					print PBS_OUT $pbs_line.$c_in{'MD_path'}."/MD_".$chrg_meth."/".$struct_b."/equi/".$struct_b.".pbs\n";;
				}
				else{
					print PBS_OUT $pbs_line.$c_in{'MD_path'}."/MD_".$chrg_meth."/".$struct_b."/".$tag."/equi/".$struct_b.".pbs\n";
				}
			}
			elsif(($pbs_line =~ m/set PATH/)||($pbs_line =~ m/set WORKDIR/)){
				if($tag eq "rec"){
					print PBS_OUT $pbs_line.$c_in{'MD_path'}."/MD_".$chrg_meth."/".$struct_b."\n";
				}
				else{
					print PBS_OUT $pbs_line.$c_in{'MD_path'}."/MD_".$chrg_meth."/".$struct_b."/".$tag."\n";
				}
			}
			elsif($pbs_line =~ m/set INPCRD/){
				if($tag eq "rec"){
					print PBS_OUT $pbs_line."../cryst/".$struct_b."_solv.crd\n";
				}
				else{
					print PBS_OUT $pbs_line."../../cryst/".$struct_b."_solv_".$tag.".crd\n";
				}
			}
			elsif($pbs_line =~m/set PRMTOP/){
				if($tag eq "rec"){
					print PBS_OUT $pbs_line."../cryst/".$struct_b."_solv.top\n";
				}
				else{
					print PBS_OUT $pbs_line."../../cryst/".$struct_b."_solv_".$tag.".top\n";
				}
			}
			elsif($pbs_line =~ m/ \-N/){
				print PBS_OUT $pbs_line." Equi_".$struct_b."\n";
			}
			elsif($pbs_line =~ m/#PBS \-o/){
				print PBS_OUT $pbs_line." Equi_".$struct_b.".OU\n";
			}
			elsif($pbs_line =~ m/#PBS \-e/){
				print PBS_OUT $pbs_line." Equi_".$struct_b.".ER\n";
			}
			#START - For SLUM script support
			elsif($pbs_line =~ m/ --job-name=/){
				print PBS_OUT $pbs_line."Equi_".$struct_b."\n";
			}
			elsif($pbs_line =~ m/ --error=/){
				print PBS_OUT $pbs_line."Equi_".$struct_b.".ER\n";
			}
			elsif($pbs_line =~ m/ --output=/){
				print PBS_OUT $pbs_line."Equi_".$struct_b.".OU\n";
			}
			#END - For SLUM script support
			else{
				print PBS_OUT $pbs_line . "\n";
			}
		}
		chmod 0755, "$equi_dir/$struct_b.pbs";
	}
}


# Generate input scripts for MD production based on template file
sub setup_MD_prod_in{
	my $MD_dir = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $nstlim;
	my $dt;

	# Determine total number of runs needed for setup of MD input files and for generation of pbs-script
	my $total_runs;

	open(PROD_TEMPL, $c_in{'prod'}) || die "Cannot open MD production template file $c_in{'prod'}.\n";

	while(my $line_prod = <PROD_TEMPL>){
		chomp($line_prod);

		if(($line_prod =~ m/nstlim/)||($line_prod =~ m/dt/)){		
			
			# Split by comma in case there is more than one flag per line
			my @line_prod = split(/\,/, $line_prod);
			
			foreach my $flag (@line_prod){
				if($flag =~ m/nstlim/){
					my @nstlim_flag = split(/\=/, $flag);
					$nstlim_flag[1] =~ m/(\d+)/;
					$nstlim = $1;
				}
				if($flag =~ m/dt/){
					my @dt_flag = split(/\=/, $flag);
					$dt_flag[1] =~ m/(\d\.\d+)/;
					$dt = $1;
				}
			}
		}
	}
	close PROD_TEMPL;

	my $runs_per_ns = 1000/($nstlim*$dt);
	$total_runs = $runs_per_ns * $c_in{'prod_t'};

	# Generate production inputs
	my $MD_in_file;

	for(my $n=1; $n<=$total_runs; $n++){
		if($n < 10){
			$MD_in_file = $MD_dir."/md_prod_00".$n.".in";
		}
		elsif($n < 100){
			$MD_in_file = $MD_dir."/md_prod_0".$n.".in";
		}
		else{
			$MD_in_file = $MD_dir."/md_prod_".$n.".in";
		}


		# Generate the scripts and modify the reference time
		my $ref_time = (($nstlim*$dt)*($n-1));
		$ref_time = $ref_time + $c_in{'equi_t'};

		open(MD_TEMPL, $c_in{'prod'}) || die "Cannot open MD production template file $c_in{'prod'}.\n";
		open(MD_IN, ">$MD_in_file") || die "Cannot open $MD_in_file for writing\n";

		while(my $md_line = <MD_TEMPL>){
			chomp($md_line);

			if(($md_line =~ /\s+t\s+/)||($md_line =~ /\s+t\=/)){
				my $new_md_line = "";
			
				# Split by comma in case there is more than one flag per line
				my @md_line = split(/\,/, $md_line);
				foreach my $flag (@md_line){
					if(($flag =~ /\s+t\s+/)||($flag =~ /\s+t\=/)){
						$flag =~ s/(\d+)/$ref_time/;
						my $check_digit = $1;
						if($check_digit eq ""){
							$flag .= $ref_time;
						}
					}
					$new_md_line .= $flag.",";
				}
				print MD_IN $new_md_line."\n";
			}
			else{
				print MD_IN $md_line."\n";
			}
		}
		close MD_TEMPL;
		close MD_IN;
	}
	return $total_runs;
}


# Create pbs-script for MD production based on PBS-template file
sub create_PBS_prod_script{
	my $MD_dir = shift;
	my $rel_MD_dir = shift;
	my $struct_b = shift;
	my $tag = shift;
	my $total_runs = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};


	open(PBS_T, $c_in{'pbs_p'}) || die "Cannot open PBS-Template for MD production $c_in{'pbs_p'}\n";
	open(PBS_NEW, ">$MD_dir/$struct_b.pbs") || die "Cannot open PBS-Script for writing\n";

	while(my $line_pbs_t = <PBS_T>){
		chomp($line_pbs_t);

		if($line_pbs_t eq "set SCRIPT="){
			if($tag eq "rec"){
				print PBS_NEW $line_pbs_t.$c_in{'MD_path'}.$rel_MD_dir."/".$tag."/prod/rec.pbs\n";
			}
			else{
				print PBS_NEW $line_pbs_t.$c_in{'MD_path'}.$rel_MD_dir."/".$struct_b."/".$tag."/prod/".$struct_b.".pbs\n";
			}
		}
		elsif($line_pbs_t eq "set BASE="){
			print PBS_NEW $line_pbs_t."md_prod\n";
		}
		elsif($line_pbs_t eq "set SY="){
			print PBS_NEW $line_pbs_t.$struct_b."\n";
		}
		elsif($line_pbs_t eq "set PRMTOP="){
			if(($c_in{'traj'} == 3)&&($tag eq "rec")){
				print PBS_NEW $line_pbs_t."../cryst/".$struct_b."_solv.top\n";
			}
			else{
				print PBS_NEW $line_pbs_t."../../cryst/".$struct_b."_solv_".$tag.".top\n";
			}
		}
		elsif($line_pbs_t eq "set WORK="){
			if($tag eq "rec"){
				print PBS_NEW $line_pbs_t.$c_in{'MD_path'}.$rel_MD_dir."/".$tag."/prod\n";
			}
			else{
				print PBS_NEW $line_pbs_t.$c_in{'MD_path'}.$rel_MD_dir."/".$struct_b."/".$tag."/prod\n";
			}
		}
		elsif($line_pbs_t eq "@ START="){
			print PBS_NEW $line_pbs_t."1\n";
		}
		elsif($line_pbs_t eq "@ END="){
			my $end = $total_runs + 1;
			print PBS_NEW $line_pbs_t.$end."\n";
		}
		elsif($line_pbs_t eq "set RESOLD="){
			print PBS_NEW $line_pbs_t."../equi/".$c_in{'MD_rst'}.".restrt\n";
		}
		elsif($line_pbs_t =~ / \-N/){
			print PBS_NEW $line_pbs_t." ".$struct_b."\n";
		}
		elsif($line_pbs_t =~ /#PBS \-o/){
			print PBS_NEW $line_pbs_t." ".$struct_b.".OU\n";
		}
		elsif($line_pbs_t =~ /#PBS \-e/){
			print PBS_NEW $line_pbs_t." ".$struct_b.".ER\n";
		}
		#START - For SLUM script support
		elsif($line_pbs_t =~ m/ --job-name=/){
			print PBS_NEW $line_pbs_t.$struct_b."\n";
		}
		elsif($line_pbs_t =~ m/ --error=/){
			print PBS_NEW $line_pbs_t.$struct_b.".ER\n";
		}
		elsif($line_pbs_t =~ m/ --output=/){
			print PBS_NEW $line_pbs_t.$struct_b.".OU\n";
		}
		#END - For SLUM script support
		else{
			print PBS_NEW $line_pbs_t . "\n";
		}
	}
	close PBS_NEW;
	close PBS_T;
	chmod 0755, "$MD_dir/$struct_b.pbs";
}

1;
