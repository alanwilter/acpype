#/usr/bin/perl

use Math::Trig;
use strict;

use Chemistry::Mol;
use Chemistry::Bond::Find ':all';
use Chemistry::Ring::Find ':all';
use lib "../additional_libs";
use FreezeThaw qw(cmpStr);
use AtomWithType;
use MOL2;
use RMSD_Kabsch;
use global;

sub gen_match_list{
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	
	# Set input-info
	my $mol1_file = $c_in{'i_path'}."/".$c_in{'Lname_start'}.".mol2";
	my $mol2_file = $c_in{'i_path'}."/".$c_in{'Lname_end'}.".mol2";
	
	my $sc1 = $c_in{'mask0'};
	my $sc2 = $c_in{'mask1'};

	# Read molecules
	my $mol1 = Chemistry::Mol->read($mol1_file);
	my $mol2 = Chemistry::Mol->read($mol2_file);

	# Determine bonding pattern / connectivity of atoms
	my ($bonded_atoms1_ref) = determine_bonded_atoms($mol1);
	my ($bonded_atoms2_ref) = determine_bonded_atoms($mol2);
	my @bonded_atoms1 = @{$bonded_atoms1_ref};
	my @bonded_atoms2 = @{$bonded_atoms2_ref};


	# Generate lists of soft core and none soft core atoms
	my (%sc1, %sc2);
	my @tmp_sc1 = split(/\@/, $sc1);
	my @tmp_sc2 = split(/\@/, $sc2);
	my @sc1_atm_names = split(/\,/, $tmp_sc1[1]);
	my @sc2_atm_names = split(/\,/, $tmp_sc2[1]);
	my (@none_sc1, @none_sc2);

	foreach my $a1 (@bonded_atoms1){
		if(grep {$a1->sprintf("%n") eq $_} @sc1_atm_names){
			$sc1{$a1} = $a1;
		}
		else{
			push(@none_sc1, $a1);
		}
	}

	foreach my $a2 (@bonded_atoms2){
		if(grep {$a2->sprintf("%n") eq $_} @sc2_atm_names){
			$sc2{$a2} = $a2;
		}
		else{
			push(@none_sc2, $a2);
		}
	}


	# Determine potential matching atoms first to reduce the search space.
	#
	# Potentially matching atoms are determined by checking the identity of the
	# types of atoms in the direct neighborhood. 
	my %match_list;
	my $match_list_ref;
	my $error_message = 0;

	foreach my $none_sc1 (@none_sc1){
		# Store atom information in matching list
		$match_list{$none_sc1}->{'mol1_atom'} = $none_sc1;

		# Get neighbors for each atom of molecule 1 that are not part of soft core
		my @neighbors1 = grep{!$sc1{$_}} $none_sc1->neighbors($none_sc1);
	
		# Determine atom types of neighboring atoms
		my @n1_type = ();
		foreach my $neighbor1 (@neighbors1){
			my $n1_type = $neighbor1->get_type();
			push(@n1_type, $n1_type);
		}
	
		foreach my $none_sc2 (@none_sc2){
			if($none_sc1->get_type() eq $none_sc2->get_type()){
				# Get neighbors for each atom of molecule 2 that are not part of soft core
				my @neighbors2 = grep{!$sc2{$_}} $none_sc2->neighbors($none_sc2);
		
				my @n2_type;
				foreach my $neighbor2 (@neighbors2){
					my $n2_type = $neighbor2->get_type();
					push(@n2_type, $n2_type);
				}
		
				# Sort types of neighboring atoms
				my (@n1_type_sorted, @n2_type_sorted);
				@n1_type_sorted = sort {$a cmp $b} @n1_type;
				@n2_type_sorted = sort {$a cmp $b} @n2_type;
		
				my $equal_type_neighbors = cmpStr(\@n1_type_sorted, \@n2_type_sorted) == 0 ? 1 : 0;
		
				if($equal_type_neighbors == 1){
					$match_list{$none_sc1}->{$none_sc2} = $none_sc2;
				}
			}
		}
	}

	# For at least one atom of molecule 1 no corresponding atom could be found in molecule 2
	my $use_elements = 0;
	foreach my $atm1 (keys %match_list){
		my $numberOfAssociated = keys %{$match_list{$atm1}};
		if(($numberOfAssociated == 1)&&($error_message == 0)){
			print "\nWARNING:\n";
			print "For atom ".$match_list{$atm1}->{'mol1_atom'}->sprintf("%n")." of molecule\n";
			print "$mol1_file\n";
			print "no atom with equivalent closest neighbors was found in molecule 2\n";
			print "$mol2_file\n"; 
			print "based on a sybyl atom type comparison.\n";
			print "You might wish to enlarge the soft core part, such that it comprises all atoms of\n";
			print "different atom type. It will now be tried to match the atoms based elements\n";
			print "instead of atom types.\n\n";
			$use_elements = 1;
		}
	}

	# In case for at least one atom of molecule 1 no corresponding atom was found 
	# in molecule 2 use elements instead of atom types for matching
	if($use_elements == 1){
		
		# Remove all elements from match_list
		foreach my $atm1 (keys %match_list){
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				delete($match_list{$atm1}->{$atm2});
			}
			delete($match_list{$atm1});
		}
		
		# Generate new list based on element matching
		foreach my $none_sc1 (@none_sc1){
			# Store atom information in matching list
			$match_list{$none_sc1}->{'mol1_atom'} = $none_sc1;

			# Get neighbors for each atom of molecule 1 that are not part of soft core
			my @neighbors1 = grep{!$sc1{$_}} $none_sc1->neighbors($none_sc1);
	
			# Determine atom types of neighboring atoms
			my @n1_element = ();
			foreach my $neighbor1 (@neighbors1){
				my $n1_element = $neighbor1->symbol();
				push(@n1_element, $n1_element);
			}
	
			foreach my $none_sc2 (@none_sc2){
				if($none_sc1->symbol() eq $none_sc2->symbol()){
					# Get neighbors for each atom of molecule 2 that are not part of soft core
					my @neighbors2 = grep{!$sc2{$_}} $none_sc2->neighbors($none_sc2);
		
					my @n2_element;
					foreach my $neighbor2 (@neighbors2){
						my $n2_element = $neighbor2->symbol();
						push(@n2_element, $n2_element);
					}
		
					# Sort types of neighboring atoms
					my (@n1_element_sorted, @n2_element_sorted);
					@n1_element_sorted = sort {$a cmp $b} @n1_element;
					@n2_element_sorted = sort {$a cmp $b} @n2_element;
		
					my $equal_element_neighbors = cmpStr(\@n1_element_sorted, \@n2_element_sorted) == 0 ? 1 : 0;
		
					if($equal_element_neighbors == 1){
						$match_list{$none_sc1}->{$none_sc2} = $none_sc2;
					}
				}
			}
		}
		
		# For at least one atom of molecule 1 no corresponding atom could be found in molecule 2
		# even if elements are used -> Wrong soft core definition.
		foreach my $atm1 (keys %match_list){
			my $numberOfAssociated = keys %{$match_list{$atm1}};
			if(($numberOfAssociated == 1)&&($error_message == 0)){
				print "\nERROR:\n";
				print "For atom ".$match_list{$atm1}->{'mol1_atom'}->sprintf("%n")." of molecule\n";
				print "$mol1_file\n";
				print "no atom with equivalent closest neighbors was found in molecule 2\n";
				print "$mol2_file\n"; 
				print "based on an element comparison.\n";
				print "You might need to enlarge the soft core part, such that it comprises all atoms that are different\n";
				print "with respect to their environment in molecules 1 and 2.\n\n";
				exit;
			}
		}
	}

	# Determine which of the atoms of molecule 2 that are not unique assigned yet
	# has an equivalent connectivity as the corresponding atom of molecule 1
	foreach my $atm1 (keys %match_list){
		my (%n1, %n2);
		my $level = 1;
		my %checked_connect;
	
		# Conduct mapping only if match list entry is not already unique
		my @no_match_list_entries = keys %{$match_list{$atm1}};
	
		if(@no_match_list_entries > 2){
	
			# Initialize neighbor list for molecule 1
			$atm1 = $match_list{$atm1}->{'mol1_atom'};
	
			foreach my $atm2 (keys %{$match_list{$atm1}}){
		
				# Do not consider entry of atm1 in match_list
				if($atm2 eq 'mol1_atom'){
					next;
				}
			
				# Which of the potentially mapping atoms has an equivalent connectivity ?
				my ($atm1_to_check, $atm2_to_check);				
				my @atm2_to_check = ();
				my (@n1,@n2);
				my %n1;
				my @checked_connect; 	# Array for storing atoms that were successfully checked
				my %map;                # Hash for storing information about neighbors of atm2
										# potentially matching the neighbors of atm1
				my $atm2_to_check_once_matched = 0;
				my $current_ref = ""; 	# Used to store atom which neighbors are currently checked. 
										# This is neccesary to ensure ordered sequential checking
										# of the atoms in the structures.

				# Initialization for the atm1 / atm2 pair
				$atm1_to_check = $atm1;
				$atm2_to_check = $match_list{$atm1}->{$atm2};
			
				while(@checked_connect < @none_sc1){
				
					# Determine neighbors for atom1 to check
					@n1 = grep{!$sc1{$_}} $atm1_to_check->neighbors($atm1_to_check);
					$n1{$atm1_to_check} = [ @n1 ];
				
					# Determine neighbors for atom2 to check
					@n2 = grep{!$sc2{$_}} $atm2_to_check->neighbors($atm2_to_check);
				
					# Compare each neighboring atom n1 with each neighboring atom n2
					my (@env1, @env2);
					my $equal_neighbors_found = 0;

					if(@n1 == @n2){				
						foreach my $n1 (@n1){
							my @pot_match_found = ();
							foreach my $n2 (@n2){
								if($use_elements == 0){		# Atom type comparison					
									if($n1->get_type eq $n2->get_type){  # Atom found potentially matching to n1
										push(@pot_match_found, $n2);
									}
								}
								if($use_elements == 1){ 	# Element comparison
									if($n1->symbol() eq $n2->symbol()){ # Atom found potentially matching to n1
										push(@pot_match_found, $n2);
									}
								}
							}
					
							if(@pot_match_found > 0){
								$equal_neighbors_found++;
								foreach my $pot (@pot_match_found){
									$map{$atm1_to_check}->{$n1}->{$pot} = $pot;
								}
							}
							else{
								last;
							}
						}
					}
				
				
					# All neighbors of checked atm2 are equivalent to the neighbors of atm1
					if($equal_neighbors_found == @n1){
				
						# Set reference, if no atom was accepted before
						if(@checked_connect == 0){
							$current_ref = $atm1_to_check;
						}
				
						# Add checked atom to list of checked if it is not already present 
						my %checked_connect;
						map{$checked_connect{$_} = 1} @checked_connect;
						if(! $checked_connect{$atm1_to_check}){
							push(@checked_connect, $atm1_to_check);
						}
					
						# In case two or more alternative atm2_to_check exist, this variable ensures
						# that match_list entry is only deleted if none atm2_to_check was matched 
						$atm2_to_check_once_matched = 1;
					}
					
					# Search was conducted for all potential atms2, but for no atm2 all atm2
					# neighbors are equivalent to the neighbors of atm1
					elsif(($equal_neighbors_found != @n1)&&(@atm2_to_check == 0)){
						if((exists $match_list{$atm1}->{$atm2})&&($atm2_to_check_once_matched == 0)){
							delete($match_list{$atm1}->{$atm2});
							last;
						}
					}
				
					# Initialize new atm1 / atm2 pair
					my $new_atm_to_check_found = 0;
					
					# Case1 : Still atom2 present that was not checked
					if(@atm2_to_check > 0){
						$atm2_to_check = shift(@atm2_to_check);
						$new_atm_to_check_found = 1;
					}
					# Case 2: Next complete initialization
					else{
						$atm2_to_check_once_matched = 0;
					
						# First check all direct neighbors of previous accepted atom, before
						# going through the list of accepted atoms again to identify unconsidered
						# neighbors.
						if(@{$n1{$current_ref}} > 0){
					
							# Initialize atm1_to_check
							$atm1_to_check = shift @{$n1{$current_ref}};
						
							# Ensure that only atoms that have not been checked
							# before are regarded as new choice for atm1_to_check
							if(grep {$atm1_to_check eq $_} @checked_connect){
								while(@{$n1{$current_ref}} > 0){
									$atm1_to_check = shift @{$n1{$current_ref}};
									if(grep {$atm1_to_check eq $_} @checked_connect){
										next;
									}
									else{
										$new_atm_to_check_found = 1;
										last;
									}
								}
							}
							else{
								$new_atm_to_check_found = 1;
							}
						
							# Initialize atm2_to_check
							if($new_atm_to_check_found == 1){
								foreach my $k (keys %{$map{$current_ref}->{$atm1_to_check}}){
									push(@atm2_to_check, $map{$current_ref}->{$atm1_to_check}->{$k});
								}
								$atm2_to_check = shift(@atm2_to_check);
							}
						}
					
						# If no direct neighbor was found for initialization, reprocess all already accepted atoms	
						if((@{$n1{$current_ref}} == 0)&&($new_atm_to_check_found == 0)){
							my @map_entries = ();
							foreach my $check (@checked_connect){
								if(@{$n1{$check}} > 0){
									$atm1_to_check = shift @{$n1{$check}};
							
									# Ensure that only atoms that have not been checked
									# before are regarded as new choice for atm1_to_check
									if(grep {$atm1_to_check eq $_} @checked_connect){
										while(@{$n1{$check}} > 0){
											$atm1_to_check = shift @{$n1{$check}};
											if(grep {$atm1_to_check eq $_} @checked_connect){
												next;
											}
											else{
												$new_atm_to_check_found = 1;
												last;
											}
										}
									}
									else{
										$new_atm_to_check_found = 1;
									}
								}
							
								# Initialize atm2_to_check
								if($new_atm_to_check_found == 1){
									$current_ref = $atm1_to_check;
									foreach my $k (keys %{$map{$check}->{$atm1_to_check}}){
										push(@atm2_to_check, $map{$check}->{$atm1_to_check}->{$k});
									}
									last;
								}
							}
							$atm2_to_check = shift(@atm2_to_check);
						}					
					}
				

					if(($new_atm_to_check_found == 1)&&(! defined $atm2_to_check)){
						$new_atm_to_check_found = 0;
					}
					if(($new_atm_to_check_found == 0)&&(@checked_connect < @none_sc1)){
						delete($match_list{$atm1}->{$atm2});
						last;
					}			
				}	
			}
		}
	}


	####################################
	# Handling none-unique assignments #
	####################################

	my $unique_assignment_found = 0;

	# Remove ambiguous entries that have already been uniquely matched to another atom of molecule 1
	($unique_assignment_found, $match_list_ref) = remove_ambiguous(\%match_list);
	%match_list = %{$match_list_ref};

	if($unique_assignment_found != 1){

		# Case 1A: Two possible solutions for heavy atom in ring
		# Case 1B: Two possible solutions for heavy atoms connected to the same atom, which
		#          has three or four direct neighbors (e.g. 2 C bonded to N)
		($error_message, $match_list_ref) = rmsd_matching_heavy($error_message, \%match_list);
		%match_list = %{$match_list_ref};
	}

	# Remove ambiguous entries that have already been uniquely matched to another atom of molecule 1
	($unique_assignment_found, $match_list_ref) = remove_ambiguous(\%match_list);
	%match_list = %{$match_list_ref};

	if($unique_assignment_found != 1){

		# Case 2A: Two possible solutions for hydrogen and oxygen atoms connected to the same atom, i.e. having 
		#          the same neighbor and neighbor has 4 direct neighbors, e.g.: -CH2- groups.
		# Case 2B: Three possible solutions for hydrogens that are part of methyl groups or oxygens that are part
		#          of PO3-groups.
		($error_message, $match_list_ref) = spatproduct_based_matching($error_message, \@sc1_atm_names, \%match_list);
		%match_list = %{$match_list_ref};
	}

	# Remove ambiguous entries that have already been uniquely matched to another atom of molecule 1
	($unique_assignment_found, $match_list_ref) = remove_ambiguous(\%match_list);
	%match_list = %{$match_list_ref};

	if($unique_assignment_found != 1){

		# Case 3: Two possible solutions for hydrogens bonded to the same neighbor, and the neighbor has 3 direct
		#         neighbors, e.g. -NH2 groups or -CH2- groups directly flanking soft core region.
		($error_message, $match_list_ref) = rmsd_matching_hydrogens($error_message, \%match_list);
		%match_list = %{$match_list_ref};
	}

	# Remove ambiguous entries that have already been uniquely matched to another atom of molecule 1
	($unique_assignment_found, $match_list_ref) = remove_ambiguous(\%match_list);
	%match_list = %{$match_list_ref};

	if($unique_assignment_found != 1){

		# Repeat spatproduct matching, if unique assignment was not found yet
		($error_message, $match_list_ref) = spatproduct_based_matching($error_message, \@sc1_atm_names, \%match_list);
		%match_list = %{$match_list_ref};
	}

	# Check assignment: Ensure that not one atom of molecule 2 was assigned to two
	# or more atoms of molecule 1
	$error_message = check_assignment($error_message, \%match_list);

	# Print final matching list
	my $match_list_file;
	if(($c_in{'match_list'})&&($c_in{'match_list'} ne "")){
		$match_list_file = $c_in{'match_list'};
	}
	else{
		$match_list_file = $c_in{'root_path'}."/TI_".$c_in{'chrg_meth'}."/".$c_in{'Lname_start'}."_".$c_in{'Lname_end'}."/match_list.txt";
	}
	write_matching_list($mol1, \@none_sc1, $match_list_file, \%match_list);
	
	my $setup_successful = ($error_message == 1) ? 0 : 1;
	return $setup_successful;
}


#################################################################################
#  Subroutines

# Subroutine for the determination of the bonding connectivity
# of the atoms based on the bonding functions provided in the
# perl subroutine Chemistry::Mol in order to identify the atoms
# involved in bonds.
sub determine_bonded_atoms{
	my $molecule = shift;

	find_bonds($molecule);
	assign_bond_orders($molecule);

	# Determine bond and atom number
	my $bond_no = $molecule->sprintf("%b");

	# Store bond ID's
	my @bonds = ();
	for(my $i=1; $i<=$bond_no; $i++){
    	push(@bonds, $molecule->bonds($i));
	}

	# Generate connection list
	my %bonded_list;
	my @bonded = ();

	foreach my $bond (@bonds){
    	my @atoms_in_bond = $bond->atoms();

    	my %b;
		map { $b{$_} = 1 } @bonded;
		foreach my $atom (@atoms_in_bond){
			if(!exists $b{$atom}){
				push(@bonded, $atom);
			}
		}
	}
	return \@bonded;
}


# Check if none uniquely identified heavy atoms exist that are part of rings and could 
# therefore have the same connectivity. These atoms need to be matched by RMSD fitting.
# Furthermore RMSD fitting is carried out for two heavy atoms with equivalent connectivity
# that are bonded to the same neighbor.
sub rmsd_matching_heavy{
	my $error_message = shift;
	my $match_list_ref = shift;
	my %match_list = %{$match_list_ref};
	my %none_unique_case1;
	foreach my $atm1 (keys %match_list){
		my @entries = keys %{$match_list{$atm1}};
		if(@entries == 3){
			my @none_unique_ring = ();
			my @none_unique_heavy = ();
			my $caseB = 0;
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				my $check_ring = find_ring($match_list{$atm1}->{$atm2});
				if(($atm2 ne 'mol1_atom')&&($check_ring != 0)){
					push(@none_unique_ring, $match_list{$atm1}->{$atm2});
				}
				if(($atm2 ne 'mol1_atom')&&($match_list{$atm1}->{$atm2}->symbol() ne "H")&&($check_ring == 0)){
					push(@none_unique_heavy,  $match_list{$atm1}->{$atm2})
				}
			}
			# Do not unambiguously matched atoms have the same neighbor (case B)?
			if(@none_unique_heavy == 2){
				my @neighbors_heavy1 = $none_unique_heavy[0]->neighbors();
				my @neighbors_heavy2 = $none_unique_heavy[1]->neighbors();
				my $found_equal_neighbor;
				foreach my $n_heavy1 (@neighbors_heavy1){
					foreach my $n_heavy2 (@neighbors_heavy2){
						if($n_heavy2 eq $n_heavy1){
							$found_equal_neighbor = $n_heavy2;
						}
					}
				}
				if($found_equal_neighbor){
					my @n_equal_neighbor = $found_equal_neighbor->neighbors();
					my $n = @n_equal_neighbor;
					if((@n_equal_neighbor == 3)||(@n_equal_neighbor == 4)){
						$caseB = 1;
					}
				}
			}
		
			# Are not unambiguously matched atoms part of ring (case A)?
			if(@none_unique_ring == 2){
				# Store all atoms of molecule 2 that fulfil the above criteria
				$none_unique_case1{$atm1}->{$none_unique_ring[0]} = $none_unique_ring[0];
				$none_unique_case1{$atm1}->{$none_unique_ring[1]} = $none_unique_ring[1];
			}
			elsif($caseB ==1){
				$none_unique_case1{$atm1}->{$none_unique_heavy[0]} = $none_unique_heavy[0];
				$none_unique_case1{$atm1}->{$none_unique_heavy[1]} = $none_unique_heavy[1];
			}
		}
	}
	
	# Search for corresponding entries in none_unique_case1 list, i.e. map entries with the 
	# same assignment solution.
	my $map_list_ref = search_corresponding(\%match_list, \%none_unique_case1);
	my %map_list = %{$map_list_ref};

	# Calculate rmsd for all possible assignments
	# Currently only assignments per atom of molecule 1 are regarded in fitting process,
	# i.e. the best assignment for each atom is chosen based on the smallest RMSD found 
	# for all possible solutions of that atom.
	# An unambigous assignment of all atoms could also be performed by regarding all 
	# 2^^n possible solutions of assignments in the fitting procedure -> Recursive analysis 
	# (Currently not implemented).
	my %pot_match;
	my $pot_match_count = 0;
	my %rmsd;
	my %wrong_assign;

	foreach my $atm1 (keys %none_unique_case1){
		foreach my $atm2 (keys %{$none_unique_case1{$atm1}}){
			$pot_match_count++;
			$pot_match{$pot_match_count}->{$atm1} = $atm2;
		
			# Generate input lists for rmsd calculation
			my (@pot_list_mol1, @pot_list_mol2); 
			foreach my $a1 (keys %match_list){
				my @entries_mol2 = keys %{$match_list{$a1}};
				if($atm1 eq $a1){
					push(@pot_list_mol1, $match_list{$atm1}->{'mol1_atom'});
					push(@pot_list_mol2, $match_list{$atm1}->{$atm2});
					
					# Assign second possible solution based on mapping list
					push(@pot_list_mol1, $match_list{$map_list{$atm1}}->{'mol1_atom'});
					foreach my $a2 (keys %{$match_list{$map_list{$atm1}}}){
						if(($a2 ne 'mol1_atom')&&($a2 ne $atm2)){
							push(@pot_list_mol2, $match_list{$map_list{$atm1}}->{$a2});
						}
					}
				}
				else{
					if(@entries_mol2 == 2){
						foreach my $a2 (keys %{$match_list{$a1}}){
							if($a2 eq 'mol1_atom'){
								push(@pot_list_mol1, $match_list{$a1}->{'mol1_atom'});
							}
							else{
								push(@pot_list_mol2, $match_list{$a1}->{$a2});
							}
						}
					}
				}
			}
		
			# Ensure that lists contain equivalent number of atoms
			my $n_atoms_list1 = @pot_list_mol1;
			my $n_atoms_list2 = @pot_list_mol2;
			if($n_atoms_list1 != $n_atoms_list2){
				print "ERROR in creation of potential matching lists. \n";
				print "Please generate the matching-list manually.\n";
				exit;
			}
		
			$rmsd{$pot_match_count} = calc_rmsd(\@pot_list_mol1, \@pot_list_mol2);
		}
	
		# Determine smallest rmsd
		my $count_equal = 0;
		my $min_rmsd = $rmsd{'1'};
		my $best_solution = 1;
	
		for(my $pot_c=2; $pot_c<=$pot_match_count; $pot_c++){
			if(($min_rmsd != 0)&&($min_rmsd == $rmsd{$pot_c})){
				$count_equal++;
			}
			if($min_rmsd >= $rmsd{$pot_c}){
				$min_rmsd = $rmsd{$pot_c};
				$best_solution = $pot_c;
			}
		}
	
		if($count_equal >= 1){
			print "WARNING: RMSD matching was ambiguous, i.e. the same RMSD was found for\n";
			print "at least two potential solutions. Please check the matching-list.\n";
		}
		
		# Record wrong assignments
		for(my $pot_c=1; $pot_c<=$pot_match_count; $pot_c++){
			if($pot_c != $best_solution){
				my $atm2 = $pot_match{$pot_c}->{$atm1};
				$wrong_assign{$atm1}->{$atm2} = $match_list{$atm1}->{$atm2};
			}
		}
		$pot_match_count = 0;
	}

	# Delete wrong assignments from matching list
	foreach my $atm1 (keys %wrong_assign){
		foreach my $atm2 (keys %{$wrong_assign{$atm1}}){
			delete($match_list{$atm1}->{$atm2});
		}
	}

	# Check for all fitted atoms if a hydrogen connectivity exists that can
	# now be unambiguously matched.
	foreach my $atm1 (keys %none_unique_case1){
		my @entries_mol2 = keys %{$match_list{$atm1}};
		if(@entries_mol2 == 2){ 	# Entry now uniquely matched ?
			# Identify direct neighbors of uniquely matched heavy atoms
			my @neighbors_atm1 = $match_list{$atm1}->{'mol1_atom'}->neighbors();
			my @neighbors_atm2 = ();
			foreach my $atm2 (@entries_mol2){
				if($atm2 ne 'mol1_atom'){
					@neighbors_atm2 = $match_list{$atm1}->{$atm2}->neighbors();
				}
			}
			# Remove all hydrogen entries of molecule 2 from matching list that are
			# not direct neighbors of matched heavy atom
			foreach my $n1 (@neighbors_atm1){
				if($n1->symbol() eq "H"){
					foreach my $associated_mol2 (keys %{$match_list{$n1}}){
						if($associated_mol2 ne 'mol1_atom'){
							if(grep {$associated_mol2 eq $_->sprintf("%i")} @neighbors_atm2){
								next;
							}
							else{
								delete($match_list{$n1}->{$associated_mol2});
							}
						}
					}
				}
			}
		}
	}
	return($error_message, \%match_list);
}


# Search for corresponding none-unique entries in match-list 
sub search_corresponding{
	my $match_list_ref = shift;
	my $none_unique_ref = shift;
	my %match_list = %{$match_list_ref};
	my %none_unique = %{$none_unique_ref};
	my @mapped;
	my %map_list;
	foreach my $atm1 (keys %none_unique){
		if((@mapped)&&(grep {$atm1 eq $_} @mapped)){
			next;
		}
		else{
			my @ambiguous_at1 = ();
			my @ambiguous_sorted_at1 = ();
			foreach my $atm2 (keys %{$none_unique{$atm1}}){
				push(@ambiguous_at1, $none_unique{$atm1}->{$atm2}->sprintf("%i"));
			}
			@ambiguous_sorted_at1 = sort {$a cmp $b} @ambiguous_at1;
					
			foreach my $at1 (keys %none_unique){
				my @ambiguous_at2 = ();
				my @ambiguous_sorted_at2 = ();
				
				if($at1 ne $atm1){
					foreach my $at2 (keys %{$none_unique{$at1}}){
						push(@ambiguous_at2, $none_unique{$at1}->{$at2}->sprintf("%i"));
					}
					@ambiguous_sorted_at2 = sort {$a cmp $b} @ambiguous_at2;

					my $equal_type_entries = cmpStr(\@ambiguous_sorted_at1, \@ambiguous_sorted_at2) == 0 ? 1 : 0;
					if($equal_type_entries == 1){
						push(@mapped, $atm1);
						push(@mapped, $at1);
						$map_list{$atm1} = $at1;
						$map_list{$at1} = $atm1;
						last;
					}
				}
			}
		}
	}
	return \%map_list;
}


# Check for all none uniquely identified hydrogens (e.g. -CH2-) and oxygens (case e.g. -SO2-) if they are 
# neighbors of the same atom. If this is the case, use the spatproduct to uniquely assign the corresponding 
# atoms.
sub spatproduct_based_matching{
	my $error_message = shift;
	my $sc1_atm_names_ref = shift;
	my @sc1_atm_names = @{$sc1_atm_names_ref};
	my $match_list_ref = shift;
	my %match_list = %{$match_list_ref};
	
	foreach my $atm1 (keys %match_list){
		my @entries = keys %{$match_list{$atm1}};
		if((@entries == 3)||(@entries == 4)){
			my @none_unique_mol2 = ();
			my $element;
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				if(($atm2 ne 'mol1_atom')&&($match_list{$atm1}->{$atm2}->symbol() eq "H")){
					$element = "H";
					push(@none_unique_mol2, $match_list{$atm1}->{$atm2});
				}
				elsif(($atm2 ne 'mol1_atom')&&($match_list{$atm1}->{$atm2}->symbol() eq "O")){
					$element = "O";
					push(@none_unique_mol2, $match_list{$atm1}->{$atm2});
				}
				else{
					if($atm2 ne 'mol1_atom'){
						$element = $match_list{$atm1}->{$atm2}->symbol();
						if($error_message == 0){
							print "ERROR: Atom of element $element was not uniquely matched.\n"; 
							print "Please check the matching list and correct it manually.\n";
							$error_message = 1;
						}
						next;
					}
				}
			}
			my ($check_neighbor, $central_atom_mol2) = check_same_neighbor_none_unique(\@none_unique_mol2);
		
			# Calculate spatproduct for not unambigously indentified atoms of molecule 2		
			if($check_neighbor == 1){
		
				# Determine central atom of molecule 1
				my $central_atom_mol1;
				my @neighbor_atom_mol1 = $match_list{$atm1}->{'mol1_atom'}->neighbors();
		
				if(@neighbor_atom_mol1 > 1){
					if($error_message == 0){
						print "For the non-uniquely matched atom ".$match_list{$atm1}->{'mol1_atom'}->sprintf("%n")." of molecule 1\n";
						print "more than one neighbor was found. The matching process was not successful. Please check and modify\n";
						print "the matching list manually.\n";
						$error_message =1;
					}
					next;
				}
				else{
					$central_atom_mol1 = shift(@neighbor_atom_mol1);
				}
		
				# Determine neighbors of central atom of molecule 1 that are not of checked element type
				my @plane_atoms_mol1 = ();
				my @none_unique_mol1 = ();
				my ($plane_atom1_mol1, $plane_atom2_mol1, $plane_atom1_mol2, $plane_atom2_mol2);
				my @neighbors_central_mol1 = $central_atom_mol1->neighbors();
				if(@neighbors_central_mol1 == 4){
					foreach my $n1 (@neighbors_central_mol1){
						if($n1->symbol() ne $element){
							push(@plane_atoms_mol1, $n1);
						}
						else{
							push(@none_unique_mol1, $n1);
						}
					}
			
					# Two neighboring plane atoms need to be found in case 2A
					if((@plane_atoms_mol1 != 2)&&(@none_unique_mol2 == 2)){
						foreach my $none_unique_mol1 (@none_unique_mol1){
							if(grep {$none_unique_mol1->sprintf("%n") eq $_} @sc1_atm_names){
								if($error_message == 0){
									print "At least one of the atoms with equivalent connectivity could not be unambiguously\n";
									print "matched, because the central atom to which it is bound directly contacts the soft\n";
									print "core part. Please closely inspect the matching list and assign the correct atoms\n";
									print "to the unambigiously matched entries.\n";
									$error_message = 1;
								}
								last;
							}
						}
						if($error_message == 0){
							print "Something went wrong in the detection of atoms in plane.\n";
							print "Please setup the matching list manually.\n";
							$error_message = 1;
						}
						next;
					}
			
					# Assignment of plane atoms for case 2B
					elsif((@plane_atoms_mol1 == 1)&&(@none_unique_mol2 == 3)){
			
						# Determination of plane atom 1
						$plane_atom1_mol1 = shift(@plane_atoms_mol1);
						
						if(grep {$plane_atom1_mol1->sprintf("%n") eq $_} @sc1_atm_names){
					 		if($error_message == 0){
								print "At least one of the atoms with equivalent connectivity could not be unambiguously\n";
								print "matched, because the central atom to which it is bound directly contacts the soft\n";
								print "core part. Please closely inspect the matching list and assign the correct atoms\n";
								print "to the unambigiously matched entries.\n";
								$error_message = 1;
							}
							next;
						}
					
						my @corresponding_plain_atom = keys %{$match_list{$plane_atom1_mol1}};
						if(@corresponding_plain_atom != 2){
							if($error_message == 0){
								print "ERROR: Problem in atom assignment procedure. Please check and manually\n";
								print "modify the created matching list\n";
								$error_message = 1;
							}
							next;
						}
						foreach my $p (keys %{$match_list{$plane_atom1_mol1}}){
							if($p ne 'mol1_atom'){
								$plane_atom1_mol2 = $match_list{$plane_atom1_mol1}->{$p};
							}
						}
				
						# Select plane atom 2 based on smallest difference in torsion angle between
						# one none_unique_atom of molecule 1 and one none_unique_atom of molecule 2.
				
						# Determine terminal torsion atom
						my @plain_atom1_neighbor = $plane_atom1_mol1->neighbors();
						my ($torsion_atom_mol1, $torsion_atom_mol2);
						foreach my $p (@plain_atom1_neighbor){
							if($p ne $central_atom_mol1){
								if(grep {$p->sprintf("%n") eq $_} @sc1_atm_names){
									next;
								}
								else{
									$torsion_atom_mol1 = $p;
								}
							}
						}
					
						if(! $torsion_atom_mol1){
					 		if($error_message == 0){
								print "At least one of the atoms with equivalent connectivity could not be unambiguously\n";
								print "matched. Please closely inspect the matching list and assign the correct atoms\n";
								print "to the unambigiously matched entries.\n";
								$error_message = 1;
							}
							next;
						}
					
						my @corresponding_mol2 = keys %{$match_list{$torsion_atom_mol1}};
						if(@corresponding_mol2 != 2){
							if($error_message == 0){
								print "ERROR: Atom selected for torsion angle calculation was not unambiguously\n";
								print "matched before. Please check and manually modify the matching list.\n";
								$error_message = 1;
							}
							next;
						}
						else{
							foreach my $atm2 (keys %{$match_list{$torsion_atom_mol1}}){
								if($atm2 ne "mol1_atom"){
									$torsion_atom_mol2 = $match_list{$torsion_atom_mol1}->{$atm2};
								}
							}
						}
				
						# Calculate torsion angles
						my (%tors_mol1, %tors_mol2);
						foreach my $none_unique_mol1 (@none_unique_mol1){
							$tors_mol1{$none_unique_mol1} = calc_torsion($none_unique_mol1, $central_atom_mol1, $plane_atom1_mol1, $torsion_atom_mol1);
						}
						foreach my $none_unique_mol2 (@none_unique_mol2){
							$tors_mol2{$none_unique_mol2} = calc_torsion($none_unique_mol2, $central_atom_mol2, $plane_atom1_mol2, $torsion_atom_mol2);
						}
						my ($error_message, $plane_atom2_mol1_ref, $plane_atom2_mol2_ref) = min_torsion($error_message, \%tors_mol1, \%tors_mol2);

						foreach my $none_unique_mol1 (@none_unique_mol1){
							if($none_unique_mol1->sprintf("%i") eq $plane_atom2_mol1_ref){
								$plane_atom2_mol1 = $none_unique_mol1;
							}
						}

						foreach my $none_unique_mol2 (@none_unique_mol2){
							if($none_unique_mol2->sprintf("%i") eq $plane_atom2_mol2_ref){
								$plane_atom2_mol2 = $none_unique_mol2;
							}
						}
					}
			
					# Assignment of plane atoms for case 2A	
					else{
						$plane_atom1_mol1 = shift(@plane_atoms_mol1);
						$plane_atom2_mol1 = shift(@plane_atoms_mol1);
	
						if((grep {$plane_atom1_mol1->sprintf("%n") eq $_} @sc1_atm_names)
					 	||(grep {$plane_atom2_mol1->sprintf("%n") eq $_} @sc1_atm_names)){
					 		if($error_message == 0){
								print "At least one of the atoms with equivalent connectivity could not be unambiguously\n";
								print "matched, because the central atom to which it is bound directly contacts the soft\n";
								print "core part. Please closely inspect the matching list and assign the correct atoms\n";
								print "to the unambigiously matched entries.\n";
								$error_message = 1;
							}
							next;
						}
				
						my @entries_plane_atom1_mol2 = keys %{$match_list{$plane_atom1_mol1}};
						my @entries_plane_atom2_mol2 = keys %{$match_list{$plane_atom2_mol1}};
					
						if((@entries_plane_atom1_mol2 == 2)&&(@entries_plane_atom2_mol2 == 2)){
							foreach my $p (keys %{$match_list{$plane_atom1_mol1}}){
								if($p ne 'mol1_atom'){
									$plane_atom1_mol2 = $match_list{$plane_atom1_mol1}->{$p};
								}
							}
				
							foreach my $p (keys %{$match_list{$plane_atom2_mol1}}){
								if($p ne 'mol1_atom'){
									$plane_atom2_mol2 = $match_list{$plane_atom2_mol1}->{$p};
								}
							}
						}	
					}
			
					# Calculate spatproduct for molecule 1 - case 2A
					my $spat_mol1;
					$spat_mol1 = spatproduct($central_atom_mol1, $plane_atom1_mol1, $plane_atom2_mol1, $match_list{$atm1}->{'mol1_atom'});

					# Calculate spatproduct for molecule 2 and delete wrong assigned atoms
					my %spat_mol2;
					# Case 2B : Plane atom 2 used as reference is equivalent to atm1 
					if($match_list{$atm1}->{'mol1_atom'} eq $plane_atom2_mol1){
						foreach my $atm2 (keys %{$match_list{$atm1}}){
							if(($match_list{$atm1}->{$atm2} ne $plane_atom2_mol2)&&($atm2 ne "mol1_atom")){
								delete($match_list{$atm1}->{$atm2});
							}	
						}
					}
					# Case 2A and 2B : Identify wrongly assigned atoms based on spat product
					else{
						foreach my $none_unique_mol2 (@none_unique_mol2){
							if($none_unique_mol2 ne $plane_atom2_mol2){
								$spat_mol2{$none_unique_mol2} = spatproduct($central_atom_mol2, $plane_atom1_mol2, $plane_atom2_mol2, $none_unique_mol2);
							}
						}
						foreach my $atm2 (keys %spat_mol2){
							if((($spat_mol2{$atm2} > 0)&&($spat_mol1 < 0))
			 	 			||(($spat_mol2{$atm2} < 0)&&($spat_mol1 > 0))){
								delete($match_list{$atm1}->{$atm2});
							}
						}
						# If case 2B, delete also plane reference from list
						if(@none_unique_mol2 == 3){
							delete($match_list{$atm1}->{$plane_atom2_mol2});
						}
					}
				}
			}
		}
	}
	return $error_message, \%match_list;
}


# Check if not uniquely identified hydrogens have the same neighbor
sub check_same_neighbor_none_unique{
	my $none_unique_ref = shift;
	my @none_unique = @{$none_unique_ref};
	my %neighbors;
	my $identical_neighbor = 1;
	my $check_neighbor;
	my $check;
	
	foreach my $none_unique (@none_unique){
		my @neighbors = $none_unique->neighbors($none_unique);
		$neighbors{$none_unique} = [ @neighbors ];
	}

	while(@{$neighbors{$none_unique[0]}} > 0){

		# Define neighboring atom to check
		$check_neighbor = shift(@{$neighbors{$none_unique[0]}});
		for my $n (1...$#none_unique){
			if(grep {$check_neighbor eq $_} @{$neighbors{$none_unique[$n]}}){
				next;
			}
			else{
				$identical_neighbor = 0;
			}
		}
		if($identical_neighbor == 1){
			$check = 1;
			last;
		}
	}
	return ($check, $check_neighbor);
}


# Measure angles between neighboring atom and hydrogens in both molecules
sub spatproduct{
	my $central_atom = shift;
	my $plane_atom1 = shift;
	my $plane_atom2 = shift;
	my $checked_atom = shift;
	
	# Assign coordinates and calculate spat-product
	my $spat;
	my %coords;
	for my $c ($checked_atom, $central_atom, $plane_atom1, $plane_atom2){
		my $coords = $c->coords;
		($coords{$c}->{'x'}, $coords{$c}->{'y'}, $coords{$c}->{'z'}) = $coords->array;
	}
	$spat = determine_spat($checked_atom, $central_atom, $plane_atom1, $plane_atom2, \%coords);
	undef(%coords);
	return $spat;
}


# Subroutine determining the spat product based on coordinates of two plane atoms, 
# one central atom, and one atom, which position shall be checked.
sub determine_spat{
	my $checked_atom = shift;
	my $central_atom = shift;
	my $plane_atom1 = shift;
	my $plane_atom2 = shift;
	my $ref_coords = shift;
	my %coords = %{$ref_coords};
	my %central_p1_vector;
	my %central_p2_vector;
	my %central_checked_vector;
	my $scalar = 0;

	# Define vectors
	for my $xyz ("x", "y", "z"){
		$central_p1_vector{$xyz} = ($coords{$plane_atom1}->{$xyz} - $coords{$central_atom}->{$xyz});
		$central_p2_vector{$xyz} = ($coords{$plane_atom2}->{$xyz} - $coords{$central_atom}->{$xyz});
		$central_checked_vector{$xyz} = ($coords{$checked_atom}->{$xyz} - $coords{$central_atom}->{$xyz});
	}

	# Determine vector product
	my %vector_prod;
	$vector_prod{'x'} = (($central_p1_vector{'y'} * $central_p2_vector{'z'}) - ($central_p1_vector{'z'} * $central_p2_vector{'y'}));
	$vector_prod{'y'} = (($central_p1_vector{'z'} * $central_p2_vector{'x'}) - ($central_p1_vector{'x'} * $central_p2_vector{'z'}));
	$vector_prod{'z'} = (($central_p1_vector{'x'} * $central_p2_vector{'y'}) - ($central_p1_vector{'y'} * $central_p2_vector{'x'}));
	
	# Calculate scalar product
	foreach my $xyz ("x", "y", "z"){
		$scalar = $scalar + ($vector_prod{$xyz} * $central_checked_vector{$xyz});
	}

	my $spat = $scalar;
	
	# Determine length of vectors
	my $length_vector_prod = 0;
	my $length_checked_vector = 0;
	
	foreach my $xyz ("x", "y", "z"){
		$length_vector_prod = ($length_vector_prod + ($vector_prod{$xyz}**2));
		$length_checked_vector = ($length_checked_vector + ($central_checked_vector{$xyz}**2));
	}

	$length_vector_prod = sqrt($length_vector_prod);
	$length_checked_vector = sqrt($length_checked_vector);
	
	return $spat;
}


# Conduct rmsd matching for two hydrogens with equivalent connectivity bonded to the same neighbor 
sub rmsd_matching_hydrogens{
	my $error_message = shift;
	my $match_list_ref = shift;
	my %match_list = %{$match_list_ref};
	my %none_unique_case3;
	
	foreach my $atm1 (keys %match_list){
		my @entries = keys %{$match_list{$atm1}};
		my $found_equal_neighbor;
		if((@entries == 3)&&($match_list{$atm1}->{'mol1_atom'}->symbol() eq "H")){
			my @none_unique_hydrogens = ();
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				if($atm2 ne 'mol1_atom'){
					push(@none_unique_hydrogens, $match_list{$atm1}->{$atm2});
				}
			}
			# Do not unambiguously matched atoms have the same neighbor?
			if(@none_unique_hydrogens == 2){
				my @neighbors_h1 = $none_unique_hydrogens[0]->neighbors();
				my @neighbors_h2 = $none_unique_hydrogens[1]->neighbors();
				foreach my $n_h1 (@neighbors_h1){
					foreach my $n_h2 (@neighbors_h2){
						if($n_h2 eq $n_h1){
							$found_equal_neighbor = $n_h2;
						}
					}
				}
			}
		
			# Are not unambiguously matched atoms part of ring (case A)?
			if($found_equal_neighbor){
				# Store all atoms of molecule 2 that fulfil the above criteria
				$none_unique_case3{$atm1}->{$none_unique_hydrogens[0]} = $none_unique_hydrogens[0];
				$none_unique_case3{$atm1}->{$none_unique_hydrogens[1]} = $none_unique_hydrogens[1];
			}
		}
	} 

	# Search for corresponding entries in none_unique_case1 list, i.e. map entries with the 
	# same assignment solution.
	my $map_list_ref = search_corresponding(\%match_list, \%none_unique_case3);
	my %map_list = %{$map_list_ref};

	# Calculate rmsd for all possible assignments
	# Currently the best assignment for an atom is chosen based on the smallest RMSD found 
	# for all possible solutions of that atom. This is possible as long as only two solutions 
	# per atom are accepted.
	my %pot_match;
	my $pot_match_count = 0;
	my %rmsd;
	my %wrong_assign;

	foreach my $atm1 (keys %none_unique_case3){
		foreach my $atm2 (keys %{$none_unique_case3{$atm1}}){
			$pot_match_count++;
			$pot_match{$pot_match_count}->{$atm1} = $atm2;
		
			# Generate input lists for rmsd calculation
			my (@pot_list_mol1, @pot_list_mol2); 
			foreach my $a1 (keys %match_list){
				my @entries_mol2 = keys %{$match_list{$a1}};
				if($atm1 eq $a1){
					push(@pot_list_mol1, $match_list{$atm1}->{'mol1_atom'});
					push(@pot_list_mol2, $match_list{$atm1}->{$atm2});
					
					# Assign second possible solution based on mapping list
					push(@pot_list_mol1, $match_list{$map_list{$atm1}}->{'mol1_atom'});
					foreach my $a2 (keys %{$match_list{$map_list{$atm1}}}){
						if(($a2 ne 'mol1_atom')&&($a2 ne $atm2)){
							push(@pot_list_mol2, $match_list{$map_list{$atm1}}->{$a2});
						}
					}
				}
				else{
					if(@entries_mol2 == 2){
						foreach my $a2 (keys %{$match_list{$a1}}){
							if($a2 eq 'mol1_atom'){
								push(@pot_list_mol1, $match_list{$a1}->{'mol1_atom'});
							}
							else{
								push(@pot_list_mol2, $match_list{$a1}->{$a2});
							}
						}
					}
				}
			}
		
			# Ensure that lists contain equivalent number of atoms
			my $n_atoms_list1 = @pot_list_mol1;
			my $n_atoms_list2 = @pot_list_mol2;
			if($n_atoms_list1 != $n_atoms_list2){
				print "ERROR in creation of potential matching lists.\n";
				print "Please generate the matching-list manually.\n";
				exit;
			}
		
			$rmsd{$pot_match_count} = calc_rmsd(\@pot_list_mol1, \@pot_list_mol2);
		}
	
		# Determine smallest rmsd
		my $count_equal = 0;
		my $min_rmsd = $rmsd{'1'};
		my $best_solution = 1;
	
		for(my $pot_c=2; $pot_c<=$pot_match_count; $pot_c++){
			if(($min_rmsd != 0)&&($min_rmsd == $rmsd{$pot_c})){
				$count_equal++;
			}
			if($min_rmsd >= $rmsd{$pot_c}){
				$min_rmsd = $rmsd{$pot_c};
				$best_solution = $pot_c;
			}
		}
	
		if(($count_equal >= 1)&&($error_message == 0)){
			print "WARNING: RMSD matching was ambiguous, i.e. the same RMSD was found for\n";
			print "at least two potential solutions. Please check the matching-list.\n";
			$error_message = 1;
		}

		# Record wrong assignments
		for(my $pot_c=1; $pot_c<=$pot_match_count; $pot_c++){
			if($pot_c != $best_solution){
				my $atm2 = $pot_match{$pot_c}->{$atm1};
				$wrong_assign{$atm1}->{$atm2} = $match_list{$atm1}->{$atm2};
			}
		}
		$pot_match_count = 0;
	}

	# Delete wrong assignments from matching list
	foreach my $atm1 (keys %wrong_assign){
		foreach my $atm2 (keys %{$wrong_assign{$atm1}}){
			delete($match_list{$atm1}->{$atm2});
		}
	}
	
	return($error_message, \%match_list);
}


# Checking the generated match list for ambiguous assignments
sub check_assignment{
	my $error_message = shift;
	my $ref_match_list = shift;
	my %match_list = %{$ref_match_list};
	my @assigned_mol2 = ();
	
	foreach my $atm1 (keys %match_list){
		my @entries = keys %{$match_list{$atm1}};
	
		# Check if there exists an assignment containing more than two atom entries for molecule 2,
		# i.e. three or more atoms with the same connectivity. In such cases the user should check
		# and manually adjust the generated matching list.
		if(@entries > 2){
			print "\n\nATTENTION: More than two atoms of molecule 2 were assigned to atom ";
			print $match_list{$atm1}->{'mol1_atom'}->sprintf("%n")." of molecule 1.\n";
			print "Please manually check the generated matching list for the correct assignment for\n";
			print "atom ".$match_list{$atm1}->{'mol1_atom'}->sprintf("%n")." of molecule 1\n";
			$error_message = 1;
		}
		
		if(@entries == 2){
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				if(($atm2 ne 'mol1_atom')&&(grep {$_ eq $atm2} @assigned_mol2)){
					print "ATTENTION: Atom ".$match_list{$atm1}->{$atm2}->sprintf("%n")." of molecule 2\n";
					print "was assigned to at least two atoms of molecule 1\n";
					print "Please manually check the generated matching list for the correct matching\n";
					print "of atom ".$match_list{$atm1}->{$atm2}->sprintf("%n")." of molecule 2\n";
				}
			}
		}
	}
	return $error_message;
}


# Calculation of the torsion angle between four sequetially connected atoms
sub calc_torsion{
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my %coords;
	my %b_a_vector;
	my %c_d_vector;
	my $length_b_a_vector = 0;
	my $length_c_d_vector = 0;
	my $scalar = 0;
	my $torsion_angle = 0;
	my $pi = 3.14159265358979;

	# Determine coordinates
	for my $atom ($a, $b, $c, $d){
		my $coords = $atom->coords;
		($coords{$atom}->{'x'}, $coords{$atom}->{'y'}, $coords{$atom}->{'z'}) = $coords->array;
	}

	# Define vectors
	for my $xyz ("x", "y", "z"){
		$b_a_vector{$xyz} = ($coords{$a}->{$xyz} - $coords{$b}->{$xyz});
		$c_d_vector{$xyz} = ($coords{$d}->{$xyz} - $coords{$c}->{$xyz});
	}

	# Determine length of vectors
	foreach my $k (keys %b_a_vector){
		$length_b_a_vector = ($length_b_a_vector + ($b_a_vector{$k}**2));
		$length_c_d_vector = ($length_c_d_vector + ($c_d_vector{$k}**2));
	}

	$length_b_a_vector = sqrt($length_b_a_vector);
	$length_c_d_vector = sqrt($length_c_d_vector);
	
	# Calculate scalar product
	foreach my $k (keys %b_a_vector){
		$scalar = $scalar + ($b_a_vector{$k} * $c_d_vector{$k});
	}

	# Determine angle
	$torsion_angle = ($scalar) / ($length_b_a_vector * $length_c_d_vector);

	$torsion_angle = acos($torsion_angle);
	$torsion_angle = ($torsion_angle/$pi) * 180;
	return $torsion_angle;
}


# Determination of the minimial difference between two sets of torsion angles
sub min_torsion{
	my $error_message = shift;
	my $tors_ref1 = shift;
	my $tors_ref2 = shift;
	my %tors1 = %{$tors_ref1};
	my %tors2 = %{$tors_ref2};
	my $init = 1;
	my $min_diff = 0;
	my $best_atom_mol1;
	my $best_atom_mol2;
	
	foreach my $t1 (keys %tors1){
		foreach my $t2 (keys %tors2){
			my $diff = abs($tors1{$t1} - $tors2{$t2});
			if($init == 1){
				$min_diff = $diff;
				$best_atom_mol1 = $t1;
				$best_atom_mol2 = $t2;
				$init = 0;
			}
			else{
				if($diff < $min_diff){
					$min_diff = $diff;
					$best_atom_mol1 = $t1;
					$best_atom_mol2 = $t2;
				}
			}
		}
	}
	
	# Ensure that the smallest difference found can be used as criterion for
	# a unique assignment
	my $count_smallest_diff = 0;
	foreach my $t1 (keys %tors1){
		foreach my $t2 (keys %tors2){
			my $diff = $tors1{$t1} - $tors2{$t2};
			if($min_diff == $diff){
				$count_smallest_diff++;
			}
		}
	}
	if($count_smallest_diff > 1){
		print "Identification of unique reference for assignment of three oxygens\n";
		print "or hydrogens bonded to the same atom, was not possible.\n";
		print "Please check and manually adjust the created matching list.\n";
		exit;
	}
	return($error_message, $best_atom_mol1, $best_atom_mol2);
}


# Remove ambiguous entries from matching list
# Deletes all entries from none-unique possible assignments for which a unique solution 
# was already found.
sub remove_ambiguous{
	my $match_list_ref = shift;
	my %match_list = %{$match_list_ref};
	my @unique_matched = ();
	my @entries_mol1 = keys %match_list;
	my $number_common_atoms = @entries_mol1;
	my $count_unique = 0;
	my $unique_assignment_found = 0;
	
	foreach my $atm1 (keys %match_list){
		my @entries_mol2 = keys %{$match_list{$atm1}};
		if(@entries_mol2 == 2){
			$count_unique++;
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				if($atm2 ne 'mol1_atom'){
					push(@unique_matched, $match_list{$atm1}->{$atm2});
				}
			}
		}
	}
	
	if($number_common_atoms == $count_unique){
		$unique_assignment_found = 1;
	}
	
	if($unique_assignment_found != 1){
		foreach my $atm1 (keys %match_list){
			my @entries_mol2 = keys %{$match_list{$atm1}};
			if(@entries_mol2 > 2){
				foreach my $atm2 (keys %{$match_list{$atm1}}){
					if(($atm2 ne 'mol1_atom')&&(grep {$match_list{$atm1}->{$atm2} eq $_} @unique_matched)){
						delete($match_list{$atm1}->{$atm2});
					}
				}
			}
		}
	}
	
	return($unique_assignment_found, \%match_list);
}

# Print matching list
sub write_matching_list{
	my $mol1 = shift;
	my $none_sc1_ref = shift;
	my @none_sc1 = @{$none_sc1_ref};
	my $match_list_file = shift;
	my $match_list_ref = shift;
	my %match_list = %{$match_list_ref};
	my $count_atoms = 0;
	
	open(MATCH_LIST, ">$match_list_file") || die "Cannot open file $match_list_file for writing the matching list\n";
	print MATCH_LIST "# Atom matching list\n";
	print MATCH_LIST "#\n";
	print MATCH_LIST "# Position in mol1"."\t"."Atom name mol1"."\t"."Atom name mol2"."\n";
	
	# Print initial matching list
	foreach my $atm1 ($mol1->atoms()){
		if(grep {$atm1 eq $_} @none_sc1){
			$count_atoms++;
			print MATCH_LIST $count_atoms."\t".$match_list{$atm1}->{'mol1_atom'}->sprintf("%n")."\t";
			foreach my $atm2 (keys %{$match_list{$atm1}}){
				if($atm2 ne 'mol1_atom'){
					print MATCH_LIST $match_list{$atm1}->{$atm2}->sprintf("%n")."\t";
				}
			}
			print MATCH_LIST "\n";
		}
	}
}

1;
