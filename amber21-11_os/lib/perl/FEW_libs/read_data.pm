# Perl module for reading data from command file

use strict;

sub read_data{
	my $c_file = shift;
	
	my %c_in;           # Hash for storing command-file information
	my $n_traj_inter;   # Number of trajectory intervals to analyze
	my %traj_inter;     # Hash for storing information about trajectory intervals
	my @traj_files;     # Array for storing names of trajectory files
	my $count = 0;

	open(C_FILE, $c_file) || die "\nCannot open command file!\n";

	# Check if correct input file was provided
	my $header = <C_FILE>;

	while(my $c_line = <C_FILE>){
		chomp($c_line);
		my $chk_comment = substr($c_line, 0, 1);

		if(($chk_comment ne "#")&&($chk_comment ne "")){
			my @c = split(/\s+/, $c_line);

			if($c[0] ne ""){
				if(($c[0] eq "tot_inter")||($c[0] eq "total_no_of_intervals")){
					$n_traj_inter = $c[1];
				}
				elsif(($n_traj_inter >= 0)&&(($c[0] eq "Start")||($c[0] eq "first_PB_snapshot"))){
					
					# If first 'Start', but not 'tot_inter' then set total 
					# number of intervals per default to 1
					if($n_traj_inter == 0){
						$n_traj_inter = 1;
					}
								
					$count++;
					$traj_inter{$count}->{'Start'} = $c[1];
				}
				elsif(($n_traj_inter > 0)&&(($c[0] eq "Stop")||($c[0] eq "last_PB_snapshot"))){
					$traj_inter{$count}->{'Stop'} = $c[1];
				}
				elsif(($n_traj_inter > 0)&&(($c[0] eq "Offset")||($c[0] eq "offset_PB_snapshots"))){
					$traj_inter{$count}->{'Offset'} = $c[1];
				}
				elsif(($c[0] eq "trj_file")||($c[0] eq "trajectory_files")){
					push(@traj_files, $c[1]);
					#START LW_03_2018
					$c_in{$c[0]} = $c[1] # $c_in{'trj_file'} would default to "all" otherwise
					#LW_03_2018
				}
				else{
					$c_in{$c[0]} = $c[1];   # Common case
				}
			}
		}
	}
	return($header, \%c_in, \@traj_files, $n_traj_inter, \%traj_inter);
}


# Read trajectories into array 'traj_files' if 'all' was specified instead of trajectories.
# File names are determined based on trajectories present in complex MD directory of ligand 1,
# since trajectories from this directory are required for 1-trajectory MM-PB(GBSA), 3-trajectory
# MM-PB(GB)SA, and LIE analyses.
sub read_trajectory_files{
	my $first_struct = shift;
	my $c_in_ref = shift;
	my %c_in = %{$c_in_ref};
	my @traj_files = ();
	
	my $gzipped = 0;
	my $traj_path = $c_in{'root_path'}."/MD_".$c_in{'chrg_meth'}."/".$first_struct."/com/prod/";
	my @mdcrd_files = <$traj_path/*.mdcrd>;
	if(@mdcrd_files == 0){
		@mdcrd_files = <$traj_path/*.mdcrd.gz>;
		$gzipped = 1;
	}
	if(@mdcrd_files == 0){
		print "\nERROR: No trajectory files were found in the directory ".$traj_path.".\n";
		print "Please ensure that all trajectories are present in the default locations\n";
		print "and have the file extension '.mdcrd' or '.mdcrd.gz', when you specify 'all'\n";
		print "under 'trj_file'\n";
		exit;
	}
			
	my @mdcrd_no;
	foreach my $mdcrd_file (@mdcrd_files){
		my @filename_components = split(/\//, $mdcrd_file);
		if($filename_components[$#filename_components] =~ m/md\_prod\_(\d+)\./){
			push(@mdcrd_no, $1);
		}
	}
			
	my @sorted_mdcrd_no = sort {$a <=> $b} @mdcrd_no;
			
	@traj_files = ();
	foreach my $no (@sorted_mdcrd_no){
		if($gzipped == 1){
			push(@traj_files, "md_prod_".$no.".mdcrd.gz");
		}
		else{
			push(@traj_files, "md_prod_".$no.".mdcrd");
		}
	}
	return \@traj_files;
}

1;
