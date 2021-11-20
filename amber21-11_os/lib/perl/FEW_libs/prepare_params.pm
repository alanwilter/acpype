use Math::Trig;
use strict;

# Assign common residue name to ligand in mol2 files in case of AM1 charge calculation
sub check_set_lig_resname_mol2{
	my $mol2_file = shift;
	my $start_reading = 0;
	my $out_str = "";
	my $res_name = "";
	
	open(MOL, $mol2_file) || die "Cannot open file $mol2_file for reading.\n";
	while(my $line = <MOL>){
		chomp($line);
		if(($start_reading == 1)&&($line =~ m/^@<TRIPOS>BOND/)){
			$start_reading = 2;
		}
		if($start_reading == 1){
			my @line = split(/\s+/, $line);			
			$res_name = $line[8];
			if($res_name ne "<1>"){
				$line =~ s/$res_name/<1>/;
			}
		}
		if(($start_reading == 2)&&($line =~ m/ROOT/)){
			$line =~ s/$res_name/<1>/;
		}
		if($line =~ m/^@<TRIPOS>ATOM/){
			$start_reading = 1;
		}
		$out_str .= $line."\n";
	}
	close(MOL);
	
	$out_str .= 
	unlink($mol2_file);
	
	open(NEW_MOL, ">$mol2_file") || die "Cannot open file $mol2_file.tmp for writing.\n";
	print NEW_MOL $out_str;
	close(NEW_MOL);	
}

# Assign common residue name to ligand in complex and ligand input PDB files
sub check_set_residue_name_ligand{
	my $pdb_file = shift;
	my $out_str = "";
	open(PDB, $pdb_file) || die "Cannot open PDB file $pdb_file for reading.\n";
	
	while(my $line = <PDB>){
		chomp($line);
		if($line =~ /ATOM|HETATM/){
			my(  $card, $atnum, $atom, $alt, $resname,    $resno,     $x, $y, $z) =
			unpack("a6      a5 x   a4     a        a3  x2     a4 x4   a8  a8  a8", $line);
			if($resname ne "<1>"){
				$line =~ s/$resname/<1>/;
			}
			$out_str .= $line."\n";
		}
	}
	close(PDB);
	
	unlink($pdb_file);
	
	open(NEW_PDB, ">$pdb_file") || die "Cannot open file $pdb_file.tmp for writing.\n";
	print NEW_PDB $out_str;
	close(NEW_PDB);
}


# Assignment of total charge of ligand molecules specified in charge-file
sub define_charge_multi{
	my $struct_names = shift;
	my $c_ref = shift;
	my %s = %{$struct_names};
	my %c = %{$c_ref};
	my %charge;
	my %multi;

    open(CHARGE, $c{'c_file'}) || die "Cannot open " . $c{'c_file'} . " with charge information.\n";

	while(my $c_line = <CHARGE>){
		chomp($c_line);
		my $comment = substr($c_line, 0, 1);
		if($comment ne "#"){
			my @c_line = split(/\t/, $c_line);

			if($c_line[0] ne ""){
				foreach my $k (keys %s){

					if($c_line[0] eq $s{$k}){
						$charge{$c_line[0]} = $c_line[1];	# Store charge
						$multi{$c_line[0]} = $c_line[2];	# Store multiplicity
					}
				}
			}
		}
	}
	close CHARGE;
	return \%charge, \%multi;
}


# Charge and multiplicity command preparation - antechamber
sub charge_multi_ante{

	my $chg = shift;
	my $command = shift;
	my $struct_name = shift;
	my $ref_charge = shift;
	my $ref_multi = shift;
	my %charge = %{$ref_charge};
	my %multi = %{$ref_multi};

	if($chg == 1){
		if($charge{$struct_name}){
			$command .= " -nc " . $charge{$struct_name};
		}
		if($multi{$struct_name}){
			$command .= " -m " . $multi{$struct_name};
		}
	}
	return $command;
}


# Setup PBS-File for gaussian calculations
sub generate_pbs_gauss{

	my $struct_b = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	open(PBS_NEW, ">$c_in{'root_path'}/gauss/$struct_b/$struct_b.pbs")
		|| die "Cannot open $c_in{'root_path'}/gauss/$struct_b/$struct_b.pbs for writing.\n";
	open(PBS_T, $c_in{'pbs_t'}) || die "Cannot open PBS-Template $c_in{'pbs_t'}.\n";

	while(my $line_t = <PBS_T>){
		chomp($line_t);

		if($line_t =~ m/set JOB/){
			my @job = split(/\=/, $line_t);
			print PBS_NEW $job[0] . "=" . $struct_b . "\n";
		}
		elsif($line_t =~ m/set WORK/){
			my @work = split(/\=/, $line_t);
			print PBS_NEW $work[0] . "=" . $c_in{'pbs_w'} . "/gauss/" . $struct_b . "\n";
		}
		elsif($line_t =~ m/ -N$/){
			print PBS_NEW $line_t." ".$struct_b."\n";
		}
		else{
			print PBS_NEW $line_t . "\n";
		}
	}				
	close PBS_NEW;
	close PBS_T;
}


# It is checked whether a hydrogen bond is present in the structure optimized
# with gaussian
sub check_hbond_gauss{
	my $gout_file = shift;
	my $struct_b = shift;
	
	my ($atom_no, $opt_type, $ref_coords, $ref_connected) = read_gout($gout_file);
	my %coords = %{$ref_coords};
	my %connected = %{$ref_connected};
	
	my $ref_distances = measure_dist($atom_no, \%coords);
	my %distances = %{$ref_distances};
	
	my $ref_closest_neighbours = identify_closest_neighbours($atom_no, $opt_type, \%connected);
	my %closest_neighbours = %{$ref_closest_neighbours};
	
	# Print out potential hydrogen bonds
	my $neighbour_N_or_O = 0;
	my $closest_equal_to_heavy = 0;
	my $angle_d_h_a = 0;
	my $warning = 0;
	
	for(my $n=1; $n<=$atom_no; $n++){
		for(my $m=1; $m<=$atom_no; $m++){
			# Distance criterion - Distance between H and N or O < 2.5 Angstrom
			if((exists $distances{$n}->{$m})&&($distances{$n}->{$m} < 2.5)){
			
				# One of the neighbouring atoms of the hydrogen participating 
				# in a hydrogen bond must be either a nitrogen or an oxygen
				for my $i (0...$#{ $closest_neighbours{$n}}){
					if(($coords{$closest_neighbours{$n}[$i]}->{'type'} == 8)||($coords{$closest_neighbours{$n}[$i]}->{'type'} == 7)){
						$neighbour_N_or_O = $closest_neighbours{$n}[$i];
					}
				}
				
				# Acceptor is not directly connected to hydrogen
				for my $i (0...$#{ $closest_neighbours{$n}}){
					if($closest_neighbours{$n}[$i] == $m){
						$closest_equal_to_heavy = 1;
					}
				}
								
				# Angle criterion
				$angle_d_h_a = determine_angle($neighbour_N_or_O, $n, $m, \%coords); 
				
				# Check hydrogen bond presence				
				if((($coords{$m}->{'type'} == 8)||($coords{$m}->{'type'} == 7))&&(($neighbour_N_or_O > 0)
				  &&($closest_equal_to_heavy != 1)&&($angle_d_h_a > 150))){
					
					if($warning == 0){
						print "\n";
						print "########################################################################\n";
						print "WARNING:\n"; 
						print "The optimized structure of Ligand $struct_b contains an intramolecular\n";
						print "hydrogen bond. Please consider that this can significantly affect the\n";
						print "calculated RESP charges and should be avoided, whenever possible.\n";
						print "########################################################################\n";
						$warning = 1;
					}
					print "Hydrogen bond (hydrogen - acceptor):  Atom $n\t Atom $m\n";
				}
			}
			$closest_equal_to_heavy = 0;
		}
		$neighbour_N_or_O = 0;
	}
	
	write_xyz($gout_file, $atom_no, \%coords);
}


# Read XYZ- and internal-coordinates from gaussian output file
sub read_gout{
	my $gout_file = shift;
	
	my $start = 0;
	my $start_read_connect = 0;
	my $check_connect;
	my $start_read_coords = 0;
	my $check_coords;
	my $end_read_coords = 0;
	my %coords;
	my $atom_no = 0;
	my %connected;
	my $gout_type = "";
	
	open(GOUT, $gout_file) || die "Cannot open file $gout_file.\n";
	while(my $gout_l = <GOUT>){
		chomp $gout_l;
		if($gout_l =~ m/Optimization completed/){
			$start = 1;
		}
		if(($start == 1)&&($gout_l =~ m/Distance matrix/)){
			$start = 0;
		}
		
		if($start == 0){
			next;
		}
		else{
			
			# Read data of optimized structure
			###################################
			
			# Identify output type - gcrt or gzmatrix output
			if($gout_l =~ m/Name/){
				if($gout_l =~ m/Definition/){
					$gout_type = "opt";
				}
				else{
					$gout_type = "popt";
				}
			}
			
			# Read connectivity from free optimization output 
			# file until "-" is present as second character of 
			# line (input: gcrt or gzmat)
			####################################################
			
			# Start criterion
			if(($gout_type eq "opt")&&($gout_l =~ m/R1 /)){
				$start_read_connect = 1;
			}

			# Stop criterion
			$check_connect = substr($gout_l, 1, 1);
			if(($start_read_connect == 1)&&($check_connect eq "-")){
				$start_read_connect = 0;
			}
			
			# Read		
			if($start_read_connect == 1){
				my @gout_l = split(/\s+/, $gout_l);
				my $connected = substr($gout_l[3], 2, -1);
				my @connected = split(/\,/, $connected);
				$connected{$gout_l[2]} = [ @connected ];
			}
			
			# Read coordinates from gcrt and gzmat output file
			###################################################
			
			# Start reading when the word "Coordinates" is present in line
			if($gout_l =~ m/Coordinates/){
				$start_read_coords = 1;
			}
			
			# Count number of "-" as second character
			my $check_coords = substr($gout_l, 1, 1);
			if(($start_read_coords == 1)&&($check_coords eq "-")){
				$end_read_coords++;
			}
			
			# End of coordinate matrix if "-" appears the second time after "Coordinates"
			if(($start_read_coords == 1)&&($end_read_coords == 2)){
				$start_read_coords = 0;
			}	
			
			# Read atom type and coordinates
			if(($start_read_coords == 1)&&(($check_coords ne "C")&&($check_coords ne "N")&&($check_coords ne "-"))){
			
				my @gout_l = split(/\s+/, $gout_l);
				
				$atom_no++;
				$coords{$gout_l[1]}->{'type'} = $gout_l[2];
				$coords{$gout_l[1]}->{'x'} = $gout_l[4];
				$coords{$gout_l[1]}->{'y'} = $gout_l[5];
				$coords{$gout_l[1]}->{'z'} = $gout_l[6];
			}
		}
	}
	close GOUT;
	
	# Read connectivity information from head of gout file if
	# partial optimization was performed.
	##########################################################
	if($gout_type eq "popt"){
		my $z_atom_no = 0;
	
		open(GOUT, $gout_file) || die "Cannot open file $gout_file.\n";
		while(my $gout_l = <GOUT>){
			chomp $gout_l;
			if($gout_l =~ m/Symbolic Z-matrix/){
				$start_read_connect = 1;
			}
		
			my @check_connect = split(/\s+/, $gout_l);
		
			if($check_connect[1] eq "Variables:"){
				$start_read_connect = 0;
			}
		
			if(($start_read_connect == 1)&&(($check_connect[1] ne "Symbolic")&&($check_connect[1] ne "Charge"))){
				my @gout_l = split(/\s+/, $gout_l);
				my @connected;
				$z_atom_no++;
				
				for(my $g=1; $g<@gout_l; $g++){
					if($g == 2){
						push(@connected, $gout_l[$g]);
					}
					if($g == 4){
						push(@connected, $gout_l[$g]);
					}
					if($g == 6){
						push(@connected, $gout_l[$g]);
					}
				}

				if($connected[0] ne ""){
					$connected{$z_atom_no} = [ @connected ];
				}
			}
		}
	}
	return($atom_no, $gout_type, \%coords, \%connected);
}


# Measure distance between hydogen atom and oxygen and nitrogen atoms
sub measure_dist{
	my $atom_no = shift;
	my $ref_coords = shift;
	my %coords = %{$ref_coords};
	my %distances;
	
	for(my $n=1; $n<=$atom_no; $n++){	
		if($coords{$n}->{'type'} == 1){
			for(my $m=1; $m<=$atom_no; $m++){
				if(($coords{$m}->{'type'} == 7)||($coords{$m}->{'type'} == 8)){
					my $diff_x = ($coords{$m}->{'x'} - $coords{$n}->{'x'})**2;
					my $diff_y = ($coords{$m}->{'y'} - $coords{$n}->{'y'})**2;
					my $diff_z = ($coords{$m}->{'z'} - $coords{$n}->{'z'})**2;
					my $distance = sqrt($diff_x + $diff_y + $diff_z);
					$distances{$n}->{$m} = $distance;
				}
			}
		}
	}
	return \%distances;
}	


# Identify atoms directly connected to atom X
sub identify_closest_neighbours{
	my $atom_no = shift;
	my $opt_type = shift;
	my $ref_connected = shift;
	my %connected = %{$ref_connected};
	my $hit = 0;
	my $intern_connect_type;
	my @neighbours;
	my %closest_neighbours;
	my $first = 0;
	my $torsion_connect = 0;
	
	if($opt_type eq "opt"){
		for(my $n=1; $n<=$atom_no; $n++){	
			foreach my $k (keys %connected){
				for my $i (0...$#{ $connected{$k}}){
					if($connected{$k}[$i] == $n){
						$hit = 1;
					}
				}
				if($hit == 1){
					$intern_connect_type = substr($k, 0, 1);
				
					# Case bond connection
					if($intern_connect_type eq "R"){
						for my $i (0...$#{ $connected{$k}}){				
							if($connected{$k}[$i] != $n){
								if(grep {$_ eq $connected{$k}[$i]} @neighbours){}
								else{
									push(@neighbours, $connected{$k}[$i]);
								}
							}
						}
					}
				}
				$hit = 0;
			}		
			$closest_neighbours{$n} = [ @neighbours ];			
			@neighbours = ();
		}
	}
	
	
	if($opt_type eq "popt"){
		for(my $n=1; $n<=$atom_no; $n++){
			for(my $k=1; $k<=$atom_no; $k++){
				if($n == $k){
					if(grep {$_ == $connected{$k}[0]} @neighbours){}
					else{
						if($connected{$k}[0] ne ""){
							push(@neighbours, $connected{$k}[0]);
						}
					}
				}
				else{
					if(grep{$_ == $n} @{$connected{$k}}){
						# Atom number in bond position
						if($connected{$k}[0] == $n){
							if(grep {$_ == $k} @neighbours){}
							else{
								push(@neighbours, $k);
							}
							if(grep {$_ == $connected{$k}[1]} @neighbours){}
							else{
								if($connected{$k}[1] ne ""){
									push(@neighbours, $connected{$k}[1]);
								}
							}
						}
					}
				}
			}
			$closest_neighbours{$n} = [ @neighbours ];
			@neighbours = ();
		}
	}		
	return \%closest_neighbours;
}


sub determine_angle{
	my $d = shift;
	my $h = shift;
	my $a = shift;
	my $ref_coords = shift;
	my %coords = %{$ref_coords};
	my %h_a_vector;
	my %h_d_vector;
	my $length_h_a_vector = 0;
	my $length_h_d_vector = 0;
	my $scalar = 0;
	my $angle = 0;
	my $pi = 3.14159265358979;

	# Define vectors
	for my $xyz ("x", "y", "z"){
		$h_a_vector{$xyz} = ($coords{$a}->{$xyz} - $coords{$h}->{$xyz});
		$h_d_vector{$xyz} = ($coords{$d}->{$xyz} - $coords{$h}->{$xyz});
	}

	# Determine length of vectors
	foreach my $k (keys %h_a_vector){
		$length_h_a_vector = ($length_h_a_vector + ($h_a_vector{$k}**2));
		$length_h_d_vector = ($length_h_d_vector + ($h_d_vector{$k}**2));
	}

	$length_h_a_vector = sqrt($length_h_a_vector);
	$length_h_d_vector = sqrt($length_h_d_vector);
	
	# Calculate scalar product
	foreach my $k (keys %h_a_vector){
		$scalar = $scalar + ($h_a_vector{$k} * $h_d_vector{$k});
	}

	# Determine angle
	$angle = ($scalar) / ($length_h_a_vector * $length_h_d_vector);

	$angle = acos($angle);
	$angle = ($angle/$pi) * 180;
	return $angle;
}


# Write XYZ
sub write_xyz{
	my $gout_file = shift;
	my $atom_no = shift;
	my $ref_coords = shift;
	my %coords = %{$ref_coords};
	
	open(XYZ, ">$gout_file.xyz") || die "Cannot open file ".$gout_file.".xyz for writing.\n";
	print XYZ "  $atom_no\n";
	print XYZ "\n"; # Comment line
	
	for(my $n=1; $n<=$atom_no; $n++){
		
		if($coords{$n}->{'type'} == 1){
			print XYZ "H";
		}
		elsif($coords{$n}->{'type'} == 6){
			print XYZ "C";
		}
		elsif($coords{$n}->{'type'} == 7){
			print XYZ "N";
		}
		elsif($coords{$n}->{'type'} == 8){
			print XYZ "O";
		}
		elsif($coords{$n}->{'type'} == 9){
			print XYZ "F";
		}
		elsif($coords{$n}->{'type'} == 15){
			print XYZ "P";
		}
		elsif($coords{$n}->{'type'} == 16){
			print XYZ "S";
		}
		elsif($coords{$n}->{'type'} == 17){
			print XYZ "Cl";
		}
		elsif($coords{$n}->{'type'} == 35){
			print XYZ "Br";
		}
		elsif($coords{$n}->{'type'} == 53){
			print XYZ "I";
		}	
		else{ 
			die "Error: XYZ-file of gaussian output cannot be generated due to unknown atom type.\n"; 
		}
		

		$coords{$atom_no}->{'x'} = sprintf("%.6f", $coords{$atom_no}->{'x'});
		$coords{$atom_no}->{'y'} = sprintf("%.6f", $coords{$atom_no}->{'y'});
		$coords{$atom_no}->{'z'} = sprintf("%.6f", $coords{$atom_no}->{'z'});
		
		printf XYZ "%13s%12s%12s\n", $coords{$n}->{'x'}, $coords{$n}->{'y'}, $coords{$n}->{'z'};	
		
	}
}


# Copying atom names from original mol2-file to ensure atom name matching,
# if setup with RESP charges is performed
sub set_atom_names{
	my $struct_path = shift;
	my $leap_dir = shift;
	my $struct_b = shift;
	
	my $mol2_gauss = $leap_dir."/".$struct_b."/".$struct_b."_resp.mol2";
	my $mol2_org = $struct_path."/".$struct_b.".mol2";
	my $temp_file = $leap_dir."/".$struct_b."/".$struct_b."_resp.mol2_tmp";
	my %atom_org;
	my $start=0;

	# Read atom names from original mol2-file
	open(ORG, $mol2_org) || die "Cannot open file $mol2_org for preparation of mol2-file for leap.\n";
	open(NEW, ">$temp_file") || die "Cannot open file $temp_file for writing of new mol2-file for leap\n";

	while(my $line_org = <ORG>){
		chomp($line_org);

		if($line_org eq "@<TRIPOS>BOND"){
			$start = 0;
		}

		if($start == 1){
			my @line_org = split(/\s+/, $line_org);

			if(($line_org[1] =~ m/(\d+)/)&&($line_org[2] =~ m/(\w)/)){
				$atom_org{$line_org[1]} = $line_org[2];
			}
		}

		if($line_org eq "@<TRIPOS>ATOM"){
			$start = 1;
		}
	}


	# Replace atom names in mol2-file generated from gaussian output
	if(! -e $mol2_gauss){
		print "Error: Mol2-file needed for further steps was not created from Gaussian output\n";
		exit;
	}
	else{
		open(G_MOL, $mol2_gauss) || die "Cannot open file $mol2_gauss for modification of atom types.\n";

		while(my $line_gmol = <G_MOL>){
			chomp($line_gmol);

			if($line_gmol eq "MOL"){
				print NEW "<1>\n";
				next;
			}

			if($line_gmol eq "@<TRIPOS>BOND"){
				$start = 0;
			}

			if($start == 1){
				my @line_gauss = split(/\s+/, $line_gmol);

				printf NEW "%7s%1s%-4s%14s%10s%10s%1s%-2s%9s%4s%14s\n", $line_gauss[1], " ", $atom_org{$line_gauss[1]},
							$line_gauss[3], $line_gauss[4], $line_gauss[5], " ", $line_gauss[6], $line_gauss[7], 
							"<1>", $line_gauss[9];
			}
			else{
				print NEW $line_gmol."\n";
			}

			if($line_gmol eq "@<TRIPOS>ATOM"){
				$start = 1;
			}
		}
		close G_MOL;
		close NEW;
					
		unlink($mol2_gauss);
		copy($temp_file, $mol2_gauss)|| die "Cannot copy $temp_file to $mol2_gauss.\n";
		unlink($temp_file);

	}
	
	unlink("esout");
	unlink("punch");
	unlink("qout");
	unlink("QOUT");
}


# Read isomer pairs for which charges shall be averaged
sub read_isomer_pairs{
	my $isomer_file = shift;
	my %iso_pairs;
	
	open(ISO, $isomer_file) || die "Cannot open file $isomer_file for reading\n";
	while(my $iso_l = <ISO>){
		chomp($iso_l);
		my $check_comment = substr($iso_l, 0, 1);
		
		if($check_comment ne "#"){
			my @iso_l = split(/\t/, $iso_l);
			
			my $check_mol2_ending = substr($iso_l[0], -5);
			if($check_mol2_ending ne ".mol2"){
				print "Steroisomer ligand structures for which charges shall be averaged\n";
				print "not correctly specified. Please make sure you give the complete\n";
				print "names of the corresponding mol2 files\n";
				exit;
			}
			
			my $iso_one = substr($iso_l[0], 0, -5);
			my $iso_two = substr($iso_l[1], 0, -5);
			
			$iso_pairs{$iso_one} = $iso_two;
		}
	}
	close ISO;
	return \%iso_pairs;
}


# Average the charges in the mol2-files generated with antechamber
# for those steroisomers for which averaging is requested
sub average_charges{
	my $leap_dir = shift;
	my $ref_iso_pairs = shift;
	my $ref_c_in = shift;
	my %iso_pairs = %{$ref_iso_pairs};
	my %c_in = %{$ref_c_in};

	# Determine charge method
	my $charge_method;
	
	if($c_in{'am1'} == 1){
		$charge_method = "am1";
	}
	if($c_in{'r_2'} == 1){
		$charge_method = "resp";
	}


	foreach my $k (keys %iso_pairs){
		my $first_file = $leap_dir."/".$k."/".$k."_".$charge_method.".mol2";
		my $second_file = $leap_dir."/".$iso_pairs{$k}."/".$iso_pairs{$k}."_".$charge_method.".mol2";
		my $start_read;
		my $end_read;
		my %charges;
		my %save_terminal_str;
		my %save_atom_str;
		my %atom_no;
		my $atom_no;
		
		if((! -e $first_file)||(! -e $second_file)){
			print "At least for one of the isomer pairs you specified for averaging no charges\n";
			print "were calculated. Please ensure that there are only stereoisomer pairs in the\n";
			print "list for charge averaging, for which charge calculation was requested.\n";
			exit;
		}
		
		# Read mol2-files of steroisomers
		for my $file ($first_file, $second_file){
			
			$start_read = 0;
			$end_read = 0;
			$atom_no = 0;
	
			open(FILE, $file) || die "Cannot open file $file for reading\n";
			
			while(my $l = <FILE>){
				chomp($l);
				if($l eq "@<TRIPOS>BOND"){
					$end_read = 1;
				}
				if(($start_read == 1)&&($end_read == 0)){
					$atom_no++;
					$charges{$atom_no}->{$file} = substr($l, -9);
					$save_atom_str{$atom_no}->{$file} = substr($l, 0, -9);
				}
				if(($start_read == 0)&&($end_read == 0)){
					$save_terminal_str{$file}->{'start'} .= $l."\n";
				}
			
				if(($start_read == 1)&&($end_read == 1)){
					$save_terminal_str{$file}->{'end'} .= $l."\n";
				}
			
				if($l eq "@<TRIPOS>ATOM"){
					$start_read = 1;
				}
			}
			close FILE;
		
			# Save atom number
			$atom_no{$file} = $atom_no;
		}
		
		# Check atom number
		if($atom_no{$first_file} != $atom_no{$second_file}){
			print "Atom numbers of $k and $iso_pairs{$k} do not match.\n";
			print "Please ensure you specified the correct steroisomer pairs\n";
			exit;
		}
		
		# Average charges and store them
		print "Charges for $k and $iso_pairs{$k} stereoisomer pair to check:\n";
		printf("%12s%12s%12s\n", $k, $iso_pairs{$k}, "Averaged");
		for(my $n=1; $n<=$atom_no; $n++){
			$charges{$n}->{'avg'} = ($charges{$n}->{$first_file} + $charges{$n}->{$second_file}) / 2;
			$charges{$n}->{'avg'} = sprintf("%.6f", $charges{$n}->{'avg'});
			printf("%12s%12s%12s\n", $charges{$n}->{$first_file}, $charges{$n}->{$second_file}, $charges{$n}->{'avg'});
		} 
		
		# Copy old mol2 files
		for my $f ($first_file, $second_file){
			my $rename_file = substr($f, 0, -5);
			$rename_file = $rename_file."_org.mol2";
			copy($f, $rename_file);
			unlink $f;	
		}
		
		# Replace charges by new charges
		for my $file ($first_file, $second_file){
			open(NEW, ">$file") || die "Cannot open file $file for writing.\n";
			print NEW $save_terminal_str{$file}->{'start'};
			for(my $n=1; $n<=$atom_no; $n++){
				printf NEW "%67s%9s\n", $save_atom_str{$n}->{$file}, $charges{$n}->{'avg'};
			}
			print NEW $save_terminal_str{$file}->{'end'};
		}
	}	
}

# Create complex from Ligand and Receptor files
sub create_complex{
	my $compl_file = shift;
	my $pdb_lig_file = shift;
	my $rec_file = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $atom_no;	
					
	open(LIG, $pdb_lig_file) || die "Cannot open PDB-file $pdb_lig_file for complex preparation.\n";
	open(REC, $rec_file) || die "Cannot open receptor file $rec_file.\n";
	open(COM, ">$compl_file") || die "Cannot open Complex-file $compl_file for writing\n";

	my $check_wat = 0;
	my @cryst_wat;	# Array for storing information about crystal waters
	my $res_no;
	my $str_com = "";
	my $last_rec_line = "";

	# Process receptor file
	while(my $rec_line = <REC>){
		chomp($rec_line);
		$rec_line =~ s/\r//;
		$last_rec_line = $rec_line;
		
		if($rec_line =~ m/END/){
			last;
		}
		elsif(($rec_line =~ m/TER/)&&($check_wat == 0)){
			$str_com .= "$rec_line\n";
			next;
		}
		elsif(($rec_line =~ m/TER/)&&($check_wat == 1)){
			push(@cryst_wat, $rec_line);
			next;
		}
		elsif(($rec_line =~ m/WAT/)||($rec_line =~ m/HOH/)){
			if(($c_in{'wat'})&&($c_in{'wat'} == 1)){
				push(@cryst_wat, $rec_line);
				$check_wat = 1;
			}
			next;
		}
		else{
			$str_com .= "$rec_line\n";
			my @rec_line = split(/\s+/, $rec_line);
			$res_no = $rec_line[4];
		}

		# Determine actual atom number
		my @rec_line = split(/\s+/, $rec_line);
		$atom_no = $rec_line[1];
	}
	close REC;

	# Check if receptor entry ends with TER otherwise 
	# insert TER before ligand entry
	if($last_rec_line !~ /TER/){
		my @split_com_str = split(/\n/, $str_com);
		# Consider case that PDB file can end with TER\nEND
		if($split_com_str[$#split_com_str - 1] != /TER/){
			$str_com .= "TER\n";
		}
	}

	# Print receptor part of complex file
	print COM $str_com;
	
	# Process ligand file
	$res_no = $res_no+1; #Residue number of ligand

	while(my $lig_line = <LIG>){
		chomp($lig_line);

		if($lig_line =~ m/ATOM/){
			$atom_no++;
			my @l_line = unpack("a6 a5 x a5 a3 x2 a4 x4 a8 a8 a8", $lig_line);

			if($l_line[2] < 5){
				printf COM "%6s%5s%1s%4s%3s%2s%4s%4s%8s%8s%8s%6s%6s\n",
							"HETATM", $atom_no, " ", $l_line[2], $l_line[3], " ", $res_no, " ", $l_line[5],
							$l_line[6], $l_line[7], "1.00", "0.00";
			}
			else{
				printf COM "%6s%5s%1s%5s%3s%2s%4s%4s%8s%8s%8s%6s%6s\n",
							"HETATM", $atom_no, " ", $l_line[2], $l_line[3], " ", $res_no, " ", $l_line[5],
							$l_line[6], $l_line[7], "1.00", "0.00";
			}
		}
	}


	# Add information about crystal waters
	if(($c_in{'wat'})&&($c_in{'wat'} == 1)){
		print COM "TER\n";

		foreach my $k (@cryst_wat){
			if($k =~ m/TER/){
				print COM $k."\n";
			}
			else{
				$atom_no++;
				my @wat_line = unpack("a6 a5 x a5 a3 a2 a4 x4 a8 a8 a8", $k);
				my $new_res_no = $wat_line[5]+1;
				if($new_res_no == 10000){
					$new_res_no = 0;
				}
				printf COM "%6s%5s%1s%5s%3s%2s%4s%12s%8s%8s%6s%6s\n",
							"HETATM", $atom_no, " ", $wat_line[2], $wat_line[3], $wat_line[4], $new_res_no, $wat_line[6],
							$wat_line[7], $wat_line[8], "1.00", "0.00";
			}
		}
	}
	
	print COM "END\n";

	close LIG;
	close COM;
}	


# Add membrane, water, and ions to receptor and complex PDB file, if 'membrane'
# was specified. Water molecules are only added, if they are further away from
# ligand than the cutoff distance 
sub add_membrane{
	my $rec_pdb_file = shift;
	my $compl_pdb_file = shift;
	my $membrane_file = shift;
	my $lig_cutoff = shift;
	my $lig_pdb_file = shift;
	
	my $rec_str_out = "";
	my $compl_str_out = "";
	my $wat_start = 0;
	my $wat_stop = 0;
	my @wat_res;
	my $keep_wat = 1;
	my %lig_coords;
	my $lig_atom_no = 0;
	
	# Read ligand coordinates
	open(LIG, $lig_pdb_file) || die "Cannot open ligand PDB file $lig_pdb_file for reading.\n";
	while(my $lig_l = <LIG>){
		chomp($lig_l);
		if(($lig_l =~ m/ATOM/)||($lig_l =~ m/HETATM/)){
			$lig_atom_no++;
		
				my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
			unpack("a6    a5  x   a4     a     a3    x1  a5   x4   a8  a8  a8", $lig_l);

			$x =~ s/^\s+//g;
			$x =~ s/\s+$//g;
			$y =~ s/^\s+//g;
			$y =~ s/\s+$//g;
			$z =~ s/^\s+//g;
			$z =~ s/\s+$//g;
		
			$lig_coords{$lig_atom_no}->{'x'} = $x;
			$lig_coords{$lig_atom_no}->{'y'} = $y;
			$lig_coords{$lig_atom_no}->{'z'} = $z;
		}
	}
	close LIG;
	
	
	# Read PDB files
	my $last_rec_entry = "";
	open(REC_PDB, $rec_pdb_file) || die "Cannot open file $rec_pdb_file for reading\n";
	while(my $pdb_l = <REC_PDB>){
		chomp($pdb_l);
		if($pdb_l =~ m/^END/){
			last;
		}
		$rec_str_out .= $pdb_l."\n";
		$last_rec_entry = $pdb_l;
	}
	close REC_PDB;
	
	my $last_compl_entry = "";
	open(COMPL_PDB, $compl_pdb_file) || die "Cannot open file $compl_pdb_file for reading\n";
	while(my $pdb_l = <COMPL_PDB>){
		chomp($pdb_l);
		if($pdb_l =~ m/^END/){
			last;
		}
		$compl_str_out .= $pdb_l."\n";
		$last_compl_entry = $pdb_l;
	}
	close COMPL_PDB;
	
	# Add TER card between old PDB part and new PDB part
	if($last_rec_entry !~ m/TER/){
		my @split_rec_str = split(/\n/, $rec_str_out);
		# Consider case that PDB file can end with TER\nEND
		if($split_rec_str[$#split_rec_str - 1] !~ /TER/){
			$rec_str_out .= "TER\n";
		}
	}
	if($last_compl_entry !~ m/TER/){
		my @split_compl_str = split(/\n/, $compl_str_out);
		# Consider case that PDB file can end with TER\nEND
		if($split_compl_str[$#split_compl_str - 1] !~ /TER/){
			$compl_str_out .= "TER\n";
		}
	}
	
	# Read membrane file
	my $mem_str_out = "";
	my $last_mem_entry = "";
	open(MEM, $membrane_file) || die "Cannot open file $membrane_file for reading.\n";
	while(my $mem_l = <MEM>){
		chomp($mem_l);
		
		if(($mem_l =~ m/WAT/)||($mem_l =~ m/HOH/)){
			$wat_start = 1;
		}
		
		if(($wat_start == 1)&&($wat_stop == 0)&&((! ($mem_l =~ m/TER/))&&(! ($mem_l =~ m/HOH/))&&(! ($mem_l =~ m/WAT/)))){
			$wat_stop = 1;
		}
		
		if(($mem_l =~ /TER/)&&($keep_wat == 0)){
			next;
		}
		
		if(($wat_start == 1)&&($wat_stop == 0)&&((($mem_l =~ m/WAT/)||($mem_l =~ m/HOH/))&&(! ($mem_l =~ m/TER/)))){
			if(@wat_res < 3){
				push(@wat_res, $mem_l);
			}
			if(@wat_res < 3){
				next;
			}
			else{
				$keep_wat = determine_water_distance($lig_cutoff, $lig_atom_no, \@wat_res, \%lig_coords);
				if($keep_wat == 1){
					foreach my $wat (@wat_res){
						$mem_str_out .= $wat."\n";
					}
					@wat_res = ();
					next;
				}
				if($keep_wat == 0){
					@wat_res = ();
					next;
				}
			}
		}
		
		if(($mem_l =~ m/ATOM/)||($mem_l =~ m/HETATM/)||($mem_l =~ m/TER/)||($mem_l =~ m/END/)){
			$mem_str_out .= $mem_l."\n";
		}
		$last_mem_entry = $mem_l;
	}
	close(MEM);
	
	unlink($rec_pdb_file);
	unlink($compl_pdb_file);
	
	# Write new files
	open(REC_WITH_MEM, ">$rec_pdb_file") || die "Cannot open file $rec_pdb_file for writing PDB with membrane.\n";
	print REC_WITH_MEM $rec_str_out;
	print REC_WITH_MEM $mem_str_out;
	if(! ($last_mem_entry =~ m/END/)){
		print REC_WITH_MEM "END\n";
	}
	close REC_WITH_MEM;
	
	open(COMPL_WITH_MEM, ">$compl_pdb_file") || die "Cannot open file $compl_pdb_file for writing PDB with membrane.\n";
	print COMPL_WITH_MEM $compl_str_out;
	print COMPL_WITH_MEM $mem_str_out;
	if(! ($last_mem_entry =~ m/END/)){
		print COMPL_WITH_MEM "END\n";
	}
	close COMPL_WITH_MEM;
}


# Determine distance between water molecule and ligand
sub determine_water_distance{
	my $lig_cutoff = shift;
	my $lig_atom_no = shift;
	my $wat_res_ref = shift;
	my $lig_coords_ref = shift;
	
	my @wat_res = @{$wat_res_ref};
	my %lig_coords = %{$lig_coords_ref};
	my $check = 1;
	my %wat_coords;
	my $temp_resno = 0;
	
	# Determine coordinates of atoms of water molecule
	my $wat_atom_count = 0;
	foreach my $wat_atom (@wat_res){
		$wat_atom_count++;
		my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
			unpack("a6    a5  x   a4     a     a3    x1  a5   x4   a8  a8  a8", $wat_atom);
			
			$x =~ s/^\s+//g;
			$x =~ s/\s+$//g;
			$y =~ s/^\s+//g;
			$y =~ s/\s+$//g;
			$z =~ s/^\s+//g;
			$z =~ s/\s+$//g;
		
			$temp_resno = $resno;
			$wat_coords{$wat_atom_count}->{'x'} = $x;
			$wat_coords{$wat_atom_count}->{'y'} = $y;
			$wat_coords{$wat_atom_count}->{'z'} = $z;
	}

	for(my $lig_atom_count=1; $lig_atom_count<=$lig_atom_no; $lig_atom_count++){
		for(my $wat_atom=1; $wat_atom <=3; $wat_atom++){
			my $diff_x = ($lig_coords{$lig_atom_count}->{'x'} - $wat_coords{$wat_atom}->{'x'})**2;
			my $diff_y = ($lig_coords{$lig_atom_count}->{'y'} - $wat_coords{$wat_atom}->{'y'})**2;
			my $diff_z = ($lig_coords{$lig_atom_count}->{'z'} - $wat_coords{$wat_atom}->{'z'})**2;
			my $distance = sqrt($diff_x + $diff_y + $diff_z);
			
			if($distance < $lig_cutoff){
				print "Deleted WAT ".$temp_resno."\n";
				$check = 0;
				last;
			}
		}
		if($check == 0){
			last;
		}
	}

	return($check);
}


# Check if frcmod file contains "ATIN, need revision" statement
# indicative of parameters that could not be determined by parmchk
sub check_frcmod{
	my $frc_file = shift;
	my $frc_check = 1;
	
	open(FRCMOD, $frc_file) || die "Cannot open file $frc_file for reading.\n";
	while(my $frc_line = <FRCMOD>){
		chomp($frc_line);
		if($frc_line =~ m/need revision/){
			print "ERROR: Missing parameters!\n";
			print "Some parameters for the ligand could not be determined by the program 'parmchk'\n";
			print "Please have a look at the entries maked as 'ATIN, need revision' in the file\n";
			print $frc_file."and set the parameters manually before proceeding with the simulation setup.\n";
			$frc_check = 0;
		}
	}
	return $frc_check;
}


# Generate leap-input file for Library creation
sub setup_leap_for_lib{
	my $lib_leap_in = shift;
	my $lib_file = shift;
	my $leap_dir = shift;
	my $chrg_meth = shift;
	my $s = shift;
	my $ref_c_in = shift;
	my %c_in = %{$ref_c_in};
	my $ref_struct_b = shift;
	my %struct_b = %{$ref_struct_b};

	# Identify residue name
	my $start_reading = 0;
	my $res_name;
	
	open(MOL, "$leap_dir/$struct_b{$s}/$struct_b{$s}_$chrg_meth.mol2") || die "Cannot open mol2-file with charges for reading\n";
	while(my $mol = <MOL>){
		chomp($mol);
		
		if($start_reading == 1){
			my @mol = split(/\s+/, $mol);
			$res_name = $mol[8];
			last;
		}
		
		if($mol eq "@<TRIPOS>ATOM"){
			$start_reading = 1;
		}
	}
	close MOL;	
	
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
	print LEAP_LIB "source leaprc.water.tip3p\n";
	print LEAP_LIB "loadAmberParams $leap_dir/$struct_b{$s}/$struct_b{$s}.frcmod\n";
	print LEAP_LIB "$res_name = loadmol2 $leap_dir/$struct_b{$s}/$struct_b{$s}_$chrg_meth.mol2\n";
	print LEAP_LIB "check $res_name\n";
	print LEAP_LIB "saveoff $res_name $lib_file\n";
	print LEAP_LIB "quit\n";
}

1;
