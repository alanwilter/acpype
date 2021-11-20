#/usr/bin/perl

use strict;
use Chemistry::Mol;

sub calc_rmsd{
	my $ref_atomlist_mol1 = shift;
	my $ref_atomlist_mol2 = shift;
	my @atomlist_mol1 = @{$ref_atomlist_mol1};
	my @atomlist_mol2 = @{$ref_atomlist_mol2};
	my $n_atoms = @atomlist_mol1;

	my (%coords1, %coords2);
	for(my $n=0; $n<$n_atoms; $n++){
		my $coords1 = $atomlist_mol1[$n]->coords;
		my $coords2 = $atomlist_mol2[$n]->coords;
		($coords1{$n}->{'x'}, $coords1{$n}->{'y'}, $coords1{$n}->{'z'}) = $coords1->array;
		($coords2{$n}->{'x'}, $coords2{$n}->{'y'}, $coords2{$n}->{'z'}) = $coords2->array;

	}

	# Determine geometric centers
	my @gcenter_mol1 = calc_geom_center($n_atoms, \%coords1);
	my @gcenter_mol2 = calc_geom_center($n_atoms, \%coords2);

	# Calculate difference between geometric centers
	my @diff;
	for(my $i=0; $i<3; $i++){
		$diff[$i] = $gcenter_mol2[$i] - $gcenter_mol1[$i];
	}

	# Determine tensor
	my @dtensor = ();
	my @dcoord_mol1 = ();
	my @dcoord_mol2 = ();

	# Initialize
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			$dtensor[$i][$j] = 0; 
		}
	}

	for(my $n=0; $n<$n_atoms; $n++){
		$dcoord_mol1[0] = $coords1{$n}->{'x'} - $gcenter_mol1[0];
		$dcoord_mol1[1] = $coords1{$n}->{'y'} - $gcenter_mol1[1];
		$dcoord_mol1[2] = $coords1{$n}->{'z'} - $gcenter_mol1[2];
		$dcoord_mol2[0] = $coords2{$n}->{'x'} - $gcenter_mol2[0];
		$dcoord_mol2[1] = $coords2{$n}->{'y'} - $gcenter_mol2[1];
		$dcoord_mol2[2] = $coords2{$n}->{'z'} - $gcenter_mol2[2];

		for(my $i=0; $i<3; $i++){
			for(my $j=0; $j<3; $j++){
				$dtensor[$i][$j] += ($dcoord_mol1[$i] * $dcoord_mol2[$j]);
			}
		}
	}

	# Generate turn matrix
	my $ref_tmatrix = turn_matrix(\@dtensor);
	my @tmatrix = @{$ref_tmatrix};

	for(my $n=0; $n<$n_atoms; $n++){
	
		# Rotate
		my ($coord_x, $coord_y, $coord_z) = rotate_point($n, $gcenter_mol2[0], $gcenter_mol2[1], $gcenter_mol2[2], \@tmatrix, \%coords2);
	
		$coords2{$n}->{'x'} = $coord_x;
		$coords2{$n}->{'y'} = $coord_y;
		$coords2{$n}->{'z'} = $coord_z;
	
    	# Translate
		$coords2{$n}->{'x'} = $coords2{$n}->{'x'} - $diff[0];
		$coords2{$n}->{'y'} = $coords2{$n}->{'y'} - $diff[1];
		$coords2{$n}->{'z'} = $coords2{$n}->{'z'} - $diff[2];	
	}

	my $rms = calc_rms($n_atoms, \%coords1, \%coords2);

	# Clean up
	@dtensor = ();
	@tmatrix = ();
	
	return $rms;
}



# Subroutine calculating the geometric center based on a set of coordinates
sub calc_geom_center{
	my $natom = shift;
	my $coords_ref = shift;
	my %coords = %{$coords_ref};
	
	my $gcenter_x = 0;
	my $gcenter_y = 0;
	my $gcenter_z = 0;
	
	foreach my $atom_id (keys %coords){
		$gcenter_x += $coords{$atom_id}->{'x'};
		$gcenter_y += $coords{$atom_id}->{'y'};
		$gcenter_z += $coords{$atom_id}->{'z'};
	}
	
	$gcenter_x = $gcenter_x / $natom;
	$gcenter_y = $gcenter_y / $natom;
	$gcenter_z = $gcenter_z / $natom;
	
	return($gcenter_x, $gcenter_y, $gcenter_z);
}

# Create turn matrix subroutine
#################################################################################################
#
# In accordance to: Th Mietzner
#
# Aim of the subroutine: Optimal adaption of positions of two molecules in 3D coordinate system
#
# Literature: W. Kabsch: A solution for the best rotation to relate two sets of vectors,
#             acta cryst. (1976) a32, 922.
#
#             W. Kabsch: A discussion of the solution for the best rotation to relate two
#             sets of vectors, acta cryst. (1976) a34,827-828. 
#
# Input parameter: Normed matrix generated from the coordinates of the molecules and a weighting
#                  vector providing information about the relative position of the molecules 
#                  with respect to each other.
# Output: Turn matrix

sub turn_matrix{
	my $dtensor_ref = shift;
	my @r = @{$dtensor_ref};
	my (@rtr, @b, @tmatrix);
	my $gen = 0.00001;
	my $fak;
	my $sum;
	
	# Initialize
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			$rtr[$i][$j] = 0;
			$b[$i][$j] = 0;
		}
	}
	
	# Creation of the help-matrix, standardization
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			for(my $izwi=0; $izwi<3; $izwi++){
				$rtr[$i][$j] = $rtr[$i][$j] + ($r[$izwi][$i] * $r[$izwi][$j]);
			}
		}
	}
	
	my ($evector_ref, $evalue_ref) = eigsys(\@rtr);
	my @evector = @{$evector_ref};
	my @evalue = @{$evalue_ref};
	
	# Create vectors B to the eigenvectors
	
	# Determine B1
	for(my $i=0; $i<3; $i++){
		$b[$i][0] = ((($r[$i][0]*$evector[0][0]) + ($r[$i][1]*$evector[1][0]) + ($r[$i][2]*$evector[2][0])) / (sqrt($evalue[0])));
	}

	# Determine B2
	if(abs($evalue[0]) > $gen) {
		my $fakt = sqrt($evalue[1]);
		for(my $i=0; $i<3; $i++){
			$b[$i][1] = ((($r[$i][0]*$evector[0][1]) + ($r[$i][1]*$evector[1][1]) + ($r[$i][2]*$evector[2][1])) / $fakt);
		}
 	}
	else{
		# Determine B2 from B1 by search for the smallest element of B1
		my @zv;
		$zv[0] = 0;
		$zv[1] = 1;
		$zv[2] = 2;
		for(my $i=1; $i<3; $i++){
			if(abs($b[$i][0]) < abs($b[0][0])) {
				my $zwi = $zv[0];
				$zv[0] = $i;
				$zv[$i] = $zwi;
			}
		}
		$b[$zv[0]][1] = 0;
		$b[$zv[1]][1] = -($b[$zv[2]][0]);
		$b[$zv[2]][1] =   $b[$zv[1]][0];
		my $fak = 1/(sqrt(($b[0][1]*$b[0][1]) + ($b[1][1]*$b[1][1]) + ($b[2][1]*$b[2][1])));
		for(my $i=0; $i<3; $i++){
			$b[$i][1] = $b[$i][1] * $fak;
		}
 	}
	
	# Determine B3
	$b[0][2] = ($b[1][0]*$b[2][1]) - ($b[2][0]*$b[1][1]);
	$b[1][2] = ($b[2][0]*$b[0][1]) - ($b[0][0]*$b[2][1]);
	$b[2][2] = ($b[0][0]*$b[1][1]) - ($b[1][0]*$b[0][1]);

	# Creation of the turn matrix
	for(my $p=0; $p<3; $p++){
		for(my $q=0; $q<3; $q++){
			$sum = 0;
			for(my $i=0; $i<3; $i++){
				$sum = $sum + ($evector[$q][$i] * $b[$p][$i]);
			}
			$tmatrix[$p][$q] = $sum;
		}
	}

	# Delete unnecessary data
	@rtr = ();
	@evector = ();
	@b = ();
	
	return \@tmatrix;
}


# Subroutine for the calculation of the eigenvectors and eigenvalues of a symetrical matrix
sub eigsys{
	my $matrix_ref = shift;
	my @r = @{$matrix_ref};
	my $gen = 0.00001;
	my $pi = 3.1415926535898;
	my $rmax = 0;
	my @evalue;
	my @evector;
	
	# Initialization of eigenvector matrix
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			$evector[$i][$j] = 0;
		}
	}
	
	# Standardization of the matrix r
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			if($rmax < abs($r[$i][$j])){
				$rmax = abs($r[$i][$j]);
			}
		}
	}
	
	if($rmax == 0){
		print "Error in subroutine eigsys: r is nullmatrix";
		exit;
	}
	
	$rmax = 1/$rmax;
	
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			$r[$i][$j] = $r[$i][$j] * $rmax;
		}
	}

	# Creation of the normal form of the cubic polynom
	my($nk1, $nk2, $nk3, $nk4);

	$nk1 = 1;

	$nk2 =  - ( $r[1][1] + $r[0][0] + $r[2][2] );

	$nk3 =  - ( $r[2][0] * $r[0][2] +
				$r[1][2] * $r[2][1] +
				$r[1][0] * $r[0][1] -
				$r[0][0] * $r[1][1] -
				$r[1][1] * $r[2][2] -
				$r[0][0] * $r[2][2] );

	$nk4 =  - ( $r[0][0] * $r[1][1] * $r[2][2] +
				$r[0][1] * $r[1][2] * $r[2][0] +
				$r[0][2] * $r[1][0] * $r[2][1] -
				$r[2][0] * $r[0][2] * $r[1][1] -
				$r[2][1] * $r[1][2] * $r[0][0] -
				$r[1][0] * $r[0][1] * $r[2][2] );

	# Generation of the reduced form of the cubic polynom
	my ($rk1, $rk2, $rk3);

	$rk1 = 1.;

	$rk2 = ((3*$nk3) - ($nk2*$nk2))/3;

	$rk3 = ((2*($nk2**3))/27) - (($nk2*$nk3)/3) + $nk4;

	# Cardano's formula for 3 real solutions
	my ($y1, $y2, $y3, $rho, $arg, $phi);

	if ($rk2 == 0){
		$y1 = 0;
		$y2 = 0;
		$y3 = 0;
	}
	else{
		$rho = sqrt((-$rk2**3)/27);
		$arg = ((-$rk3/(2*$rho)) < 1.0) ? (-$rk3/(2*$rho)) : 1.0;
		$arg = ($arg > -1.0) ? $arg : -1.0;
		$phi = acos($arg);
		$y1 = (2*($rho**(1/3))) * cos($phi/3);
		$y2 = (2*($rho**(1/3))) * cos(($phi/3)+((2*$pi)/3));
		$y3 = (2*($rho**(1/3))) * cos(($phi/3)+((4*$pi)/3));
	}

	# Null terms of the normal form
	$evalue[0] = $y1 - ($nk2/3);
	$evalue[1] = $y2 - ($nk2/3);
	$evalue[2] = $y3 - ($nk2/3);
	
	# Calculation of the eigenvectors of matrix r
	my @zv = ();
	
	# Sorting eigenvalues according to size
	for(my $i=0; $i<2; $i++){
		for(my $j=1; $j<3; $j++){
			if (abs($evalue[$i]) < abs($evalue[$j])){
				my $zwi = $evalue[$i];
				$evalue[$i] = $evalue[$j];
				$evalue[$j] = $zwi;
			}
		}
	}

	# Determine 1. eigenvector
	my $eig = $evalue[0];
	my $m = 0;
	my $kenn = 0;
	my $evector_gauss_ref;
	
	($m, $kenn, $evector_gauss_ref) = gauss($eig, $gen, $m, $kenn, \@evector, \@r);
	@evector = @{$evector_gauss_ref};
	
	# Standardization of first eigenvector
	my $fak = 1/sqrt(($evector[0][0]*$evector[0][0]) + ($evector[1][0]*$evector[1][0]) + ($evector[2][0]*$evector[2][0]));
	for(my $i=0; $i<3; $i++){
		$evector[$i][0] = $evector[$i][0] * $fak;
	}
	
	# Determine 2. eigenvector ...
	if($kenn != 0){
		my $sp = ($evector[0][0] * $evector[0][1]) + ($evector[1][0] * $evector[1][1]) + ($evector[2][0] * $evector[2][1]);
		for(my $i=0; $i<3; $i++){
			$evector[$i][1] = $evector[$i][1] - ($sp * $evector[$i][0]);
		}
	}
	else{
		if (abs($evalue[1]) > $gen){
			# ... from the second eigenvalue
			$eig = $evalue[1];
			$m = 1;
			($m, $kenn, $evector_gauss_ref) = gauss($eig, $gen, $m, $kenn, \@evector, \@r);
			@evector = @{$evector_gauss_ref};
		}
		else{
			# ... or from the first eigenvalue -> Search for smallest element of first eigenvalue
			$zv[0] = 0;
			$zv[1] = 1;
			$zv[2] = 2;
			
			for(my $i=1; $i<3; $i++){
				if(abs($evector[$i][0]) < abs($evector[0][0])){
					my $zwi = $zv[0];
					$zv[0] = $i;
					$zv[$i] = int($zwi);
				}
			}
			$evector[$zv[0]][1] = 0;
			$evector[$zv[1]][1] = -($evector[$zv[2]][0]);
			$evector[$zv[2]][1] = $evector[$zv[1]][0];
		}
	}

	# Standardize the 2. eigenvector
	my $fak = 1/(sqrt(($evector[0][1]*$evector[0][1]) + ($evector[1][1]*$evector[1][1]) + ($evector[2][1]*$evector[2][1])));
	for(my $i=0; $i<3; $i++){
		$evector[$i][1] = $evector[$i][1] * $fak;
	}
	
	# 3. eigenvector : vectorproduct of the two other eigenvectors
	$evector[0][2] = ($evector[1][0]*$evector[2][1]) - ($evector[2][0]*$evector[1][1]);
	$evector[1][2] = ($evector[2][0]*$evector[0][1]) - ($evector[0][0]*$evector[2][1]);
	$evector[2][2] = ($evector[0][0]*$evector[1][1]) - ($evector[1][0]*$evector[0][1]);

	for(my $p=0; $p<3; $p++){
		$evalue[$p] = $evalue[$p]/$rmax;
	}

	return(\@evector, \@evalue);
}


# This subroutine calculates for a symmetric, positive semidefinite matrix
# with eigenvalue EIG the corresponding eigenspace by factorization according
# to Gauss and complete Pivot search
sub gauss{
	my $eig = shift;
	my $gen = shift;
	my $m = shift;
	my $kenn = shift;
	my $evector_ref = shift;
	my $r_ref = shift;
	my @evector = @{$evector_ref};
	my @r = @{$r_ref};
	my (@tau, @a);
	
	for(my $j=0; $j<3; $j++){
		$tau[$j] = $j;
	}
		
	for(my $i=0; $i<3; $i++){
		for(my $j=0; $j<3; $j++){
			$a[$i][$j] = $r[$i][$j];
		}
		$a[$i][$i] = $a[$i][$i] - $eig;
	}

	my $zeil = 3;
		
	# Total Pivot search and triangle shape of matrix a
	for(my $k=0; $k<3; $k++){

		# Search for largest element
		my $izwi = $k;
		my $jzwi = $k;
		my $gelem = abs($a[$k][$k]);

		for(my $i=$k; $i<3; $i++){
			for(my $j=$k; $j<3; $j++){
				if(abs($a[$i][$j]) > $gelem){
					$izwi = $i;
					$jzwi = $j;
					$gelem = abs($a[$i][$j]);
				}
			}
		}
			
		if($gelem <= $gen){
			last;	
		}
		else{
			$zeil = $zeil - 1;

			# Zeilentausch
			if($izwi != $k){
				for(my $j=0; $j<3; $j++){
					my $zwi = $a[$k][$j];
					$a[$k][$j] = $a[$izwi][$j];
					$a[$izwi][$j] = $zwi;
				}
			}

			# Spaltentausch
			if($jzwi != $k){
				for(my $i=0; $i<3; $i++){
					my $zwi = $a[$i][$k];
					$a[$i][$k] = $a[$i][$jzwi];
					$a[$i][$jzwi] = $zwi;
				}
				my $tauzwi = $tau[$k];
				$tau[$k] = $tau[$jzwi];
				$tau[$jzwi] = $tauzwi;
			}
				
			# Triangular shape
			for(my $l=$k+1; $l<3; $l++){
				my $zwi = $a[$l][$k]/$a[$k][$k];
				for(my $j=$k; $j<3; $j++){
					$a[$l][$j] = $a[$l][$j] - ($zwi * $a[$k][$j]);
				}
			}
		}
	}

	# Solving equation system
		
	# Rang A = 0; null matrix
	if ($zeil == 3){
		for(my $i=0; $i<3; $i++){
			for(my $j=0; $j<3; $j++){
				$evector[$i][$j] = 0;
			}
			$evector[$i][$i] = 1.;
		}
		return($m, $kenn, \@evector);
	}
		
	# Rang A = 1; two null lines
	if($zeil == 2){
		$evector[0][$m] = -($a[0][2]/$a[0][0]);
		$evector[1][$m] = 0;
		$evector[2][$m] = 1;
		$evector[0][$m+1] = -($a[0][1]/$a[0][0]);
		$evector[1][$m+1] = 1;
		$evector[2][$m+1] = 0;
		$kenn = 1;
	}

	# Rang A = 2; one null line
	if($zeil == 1){
		$evector[1][$m] = -($a[1][2]/$a[1][1]);
		$evector[0][$m] = -((($a[0][1]*$evector[1][$m]) + $a[0][2])/$a[0][0]);
		$evector[2][$m] = 1;
	}
		
	# Rang A = 3; not possible, since (M - lambda*E)*evector = 0, evector != 0
	if ($zeil == 0){
		print "Error in subroutine gauss: No null lines.\n";
		exit;
	}

	# Sort eigenvector
	my $terminate = 0;
	while($terminate == 0){
		$a[$tau[0]][0] = $evector[0][$m];
		$a[$tau[1]][0] = $evector[1][$m];
		$a[$tau[2]][0] = $evector[2][$m];
		$evector[0][$m] = $a[0][0];
		$evector[1][$m] = $a[1][0];
		$evector[2][$m] = $a[2][0];

		if(!(($m == 0) && $kenn == 1)){
			last;
		}
		$m = 1;
	}
	return($m, $kenn, \@evector);
}


# Subroutine for changing coordinates according to the rotate matrix
sub rotate_point{
	my $atom_id = shift;
	my $ref_x = shift;
	my $ref_y = shift;
	my $ref_z = shift;
	my $tmatrix_ref_help = shift;
	my @tmatrix_help = @{$tmatrix_ref_help};
	my $coords2_ref = shift;
	my %coords2 = %{$coords2_ref};
	my $coord_x = $coords2{$atom_id}->{'x'};
	my $coord_y = $coords2{$atom_id}->{'y'};
	my $coord_z = $coords2{$atom_id}->{'z'};

	my $x = $coord_x - $ref_x;
	my $y = $coord_y - $ref_y;
	my $z = $coord_z - $ref_z;

	$coord_x = $tmatrix_help[0][0] * $x + $tmatrix_help[0][1] * $y + $tmatrix_help[0][2] * $z;
	$coord_y = $tmatrix_help[1][0] * $x + $tmatrix_help[1][1] * $y + $tmatrix_help[1][2] * $z;
	$coord_z = $tmatrix_help[2][0] * $x + $tmatrix_help[2][1] * $y + $tmatrix_help[2][2] * $z;
	
	$coord_x += $ref_x;
	$coord_y += $ref_y;
	$coord_z += $ref_z;

	return($coord_x, $coord_y, $coord_z);
}


# Subroutine for the calculation of the RMSD
sub calc_rms{
	my $n_atom = shift;
	my $coords1_ref = shift;
	my $coords2_ref = shift;
	my %coords1 = %{$coords1_ref};
	my %coords2 = %{$coords2_ref};
	my $rms = 0;
	
	for(my $n=0; $n<$n_atom; $n++){
		my $sqr_distx = ($coords1{$n}->{'x'} - $coords2{$n}->{'x'} == 0) ? 0 : (($coords1{$n}->{'x'} - $coords2{$n}->{'x'})**2);
		my $sqr_disty = ($coords1{$n}->{'y'} - $coords2{$n}->{'y'} == 0) ? 0 : (($coords1{$n}->{'y'} - $coords2{$n}->{'y'})**2);
		my $sqr_distz = ($coords1{$n}->{'z'} - $coords2{$n}->{'z'} == 0) ? 0 : (($coords1{$n}->{'z'} - $coords2{$n}->{'z'})**2);
		$rms = $rms + $sqr_distx + $sqr_disty + $sqr_distz;
	}
	
	$rms = sqrt($rms / $n_atom);
	$rms = sprintf("%.4f", $rms);
	return $rms;
}

1;
