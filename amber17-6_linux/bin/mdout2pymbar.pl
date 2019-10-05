#!/usr/bin/perl -w
# Amber output to PyMBAR input conversion script. Written by Thomas Steinbrecher 05032010
# Use at your own risk

$calorie = 4.184;

open FILE, $ARGV[0];

print STDERR "Processing file $ARGV[0]\n";
print "# Amber FEP energy data extracted from file $ARGV[0]\n";
print "#\n# Multiple Found XYZ messages in the following are ok\n";

# Figure out dt and bar_intervall

$dt=0;
$intervall=0;

while ($line=<FILE>) {
    if ($line =~ /dt\s+=\s+(\d+.\d+)/) {
	$dt = $1;
	print "# Found timestep = $dt fs\n";
    }
    if ($line =~ /bar_intervall\s+=\s+(\d+.\d+)/) {
	$intervall = $1;
	print "# Found bar_intervall = $intervall steps\n";
    }
    if ($line =~ /clambda\s+=\s+(\d+.\d+)/) {
	$clambda = $1;
	print "# Found lambda = $clambda\n";
    }
    if ($line =~ /TOTAL SIZE OF NONBOND/) {
	last;
    }
}

#cycle through the output file and collect FEP energy values

# on first go, figure out number of lambdas and compare to 'real' lambda
while ($line=<FILE>) {
    if ($line =~ /MBAR Energy analysis/) {
	$i=0;
	while ( ($line=<FILE>) =~ /Energy at (\d+.\d+) =\s*(-?\d+.\d+)/) {
	    $energy[0][$i] = $2;
	    $lambda_value[$i++] = $1;
	}
	$lambda_max=$i; # one higher than actual number of points
	print "# Found total of $lambda_max l-values\n";
	last;
    }
}

for ($i=0 ; $i<$lambda_max ; $i++) {
    if ($lambda_value[$i] == $clambda) {
	printf "# This output file is assumed to be number %4i of $lambda_max\n", $i+1;
	$position = $i;
	last;
    }
}
    
# now go through the rest and check if the right number of lambdas is found on each block
$step=1;
while ($line=<FILE>) {
    if ($line =~ /MBAR Energy analysis/) {
	$i=0;
	while ( ($line=<FILE>) =~ /Energy at \d+.\d+ =\s*(-?\d+.\d+)/) {
	    $energy[$step][$i++] = $1;
	}
	if ($i == $lambda_max) {
	    #ok
	}else {
	    die ("Lambda value count differs on block $step\n");
	}
	$step++;
    }
}
print STDERR "Found $step output blocks\n";

#output all the values found as differences to the current energy, in kJ/mol

print "# Energy values:";

$time=0;
for ($i=0 ; $i< $step ; $i++) {
    $time += $dt*$intervall/1000;
    printf "\n%12.4f%14.6f", $time, 0.0;
    for ($j=0 ; $j < $lambda_max ; $j++) {
	printf "%16.7e", ($energy[$i][$j]-$energy[$i][$position])*$calorie ;
	#printf "%16.7e", $energy[$i][$j] ;
    }
    print "   1.0000000e+00"; # no volume data, leave no final newline
}
