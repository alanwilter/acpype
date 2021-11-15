use strict;

# Separation of multi-ligand file
# Single mol2 structures are generated
# from multi-ligand mol2 file
sub separate{
    my $c_ref = shift;
    my %c = %{$c_ref};
	
	# Check if file extension of ligand file was provided
	my $ext = substr($c{'fi_lig'}, -4);
	my $fi_lig;
	
	if($ext eq "mol2"){
		$fi_lig = substr($c{'fi_lig'}, 0, -5);
	}
	elsif($ext eq ".sdf"){
		$fi_lig = substr($c{'fi_lig'}, 0, -4);
	}
	else{
		$fi_lig = $c{'fi_lig'};
	}
	
    my $file_base = $c{'i_path'} . "/" . $fi_lig;

    # Create directory for separated structures
	my $sep_dir = $c{'root_path'} . "/structs";
	mkdir($sep_dir);

	# Process the ligand file
	my $mol = 0;
	my $first = 0;
	my $n_lig = 0;
	my @names;         # Array for storing names of ligands specified in mol2

	open(SEP_FILE, "$file_base.mol2") || die "Cannot open file $file_base.mol2 for structure separation.\n";

	while(my $s_line = <SEP_FILE>){
		chomp($s_line);
		$s_line =~ s/\r//;

		if($mol == 1){
			if($s_line ne ""){  # Use names for file specification
				# Ensure that names do not contain spaces; replace them by "_"
				$s_line =~ s/(\s+)/_/g;
				# Remove brackets
				$s_line =~ s/\(//g;
				$s_line =~ s/\)//g;
			
				if($first == 0){
					push(@names, $s_line);
					open(STRUCT, ">$sep_dir/Lig_$s_line.mol2") || die "Cannot open ligand structure file " . $sep_dir .
																		  "/Lig_" . $s_line .".mol2 for writing.";
				}
				if($first == 1){
					# Check if name was already used
					if(grep {$s_line eq $_} @names){
						print "Error: Ligand names are not unique.\n";
						exit;
					}
					else{
						open(STRUCT, ">$sep_dir/Lig_$s_line.mol2") || die "Cannot open ligand structure file " . $sep_dir .
																		  "/Lig_" . $s_line .".mol2 for writing.";
					}
				}
			}
			else{	# Automatic numbering of files if no name is specified
				$n_lig++;
				open(STRUCT, ">$sep_dir/Lig_$n_lig.mol2")|| die "Cannot open ligand structure file " . $sep_dir .
																 "/Lig_" . $n_lig . ".mol2 for writing.\n";
				$s_line = "Lig_$n_lig";
			}
			print STRUCT "@<TRIPOS>MOLECULE\n";
			$mol = 0;
			$first = 1;
		}

		if($s_line eq "@<TRIPOS>MOLECULE"){
		    close(STRUCT);
            $mol = 1;
		}

		if(($mol == 0)&&($first == 1)){
			print STRUCT "$s_line\n";
		}
	}
}
1;
