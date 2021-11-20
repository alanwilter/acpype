#!/usr/bin/env perl

# Main perl-script for automated calculation of free energies of binding for ligands
# with tools available in AMBER16

use FindBin;
use File::Find;

BEGIN{
	my $amber_home =  $ENV{"AMBERHOME"};
	if($amber_home eq ""){
        my $current_dir = $FindBin::Bin;
        print "AMBERHOME is not set\n";
        if($current_dir =~ m/bin/){
            $amber_home = $current_dir;
            $amber_home =~ s/bin//g;
            print "Trying to set AMBERHOME to $amber_home\n";
            $ENV{'AMBERHOME'} = $amber_home;
        }else{
            exit
        }
	}
	else{
		# Search for system specific installation path of PerlMol modules
		my $perlmol_lib_path = $amber_home."/lib/perl/";
		find(\&sub_dir_names, $perlmol_lib_path);
		sub sub_dir_names{
			my $file_or_dir = $_;
			my $curr_path_name = $File::Find::name;
			if((-d $file_or_dir)&&($file_or_dir eq "Chemistry")){
				my @lib_path = split(/\//, $curr_path_name);
				my $perl_five_lib_path;
				for my $sub_path (0...(@lib_path-2)){
					$perl_five_lib_path .= $lib_path[$sub_path] ."/";
				}
				if(-e $perl_five_lib_path){
					unshift @INC, $perl_five_lib_path;
				}
			}
		}
	}
}

# after checking/setting AMBERHOME
use lib "$ENV{'AMBERHOME'}/lib/perl/lib/perl5/site_perl";
use lib "$ENV{'AMBERHOME'}/lib/perl/FEW_libs";

my $current_dir = $FindBin::Bin;

######################################
# Check presence of required modules # 
######################################

# 1. Modules installed with perl and not provided with FEW
my @modules = ("File::Copy", "File::Path");

foreach my $module (@modules){
	if(eval "require $module" ne 1){
		print "\nThe Perl module '$module' required for execution of FEW was not found.\n";
		print "Please install this module before you run FEW.\n\n";
		exit;
	}
}

if(eval "require Chemistry::Mol" ne 1){
	my $amber_home =  $ENV{"AMBERHOME"};
	my $chem_module_path = $amber_home."/AmberTools/src/FEW/additional_libs/PerlMol";
	my $chem_bundle_path = $amber_home."/AmberTools/src/FEW/additional_libs/PerlMol-0.3500";
	print "The Perl module 'Chemistry::Mol' required for execution of FEW was not found.\n";
	print "Please check whether the PerlMol modules were installed under:\n";
	print $chem_module_path."\n";
	print "If this is not the case please install them under this location or somewhere on\n";
	print "your system following the instructions at http://www.perlmol.org/download.html\n";
	print "A version of PerlMol compatible with FEW ready for installation can be found under:\n";
	print "$chem_bundle_path.\n\n";
	print "==================================================================================\n";
	exit;
}


# 2. Modules distributed with FEW
my @addlib_files = ("AtomWithType", "Descriptive", "Distributions", "FreezeThaw", "PointEstimation", "ReadBackwards", "Smoother", "Normality");
my @pm_files = ("atom_matching", "check_input", "common_prepare_MD", "global", "LIEW", "MOL2", "prepare_params", "prepare_top_crd_MD",
				"read_data", "RMSD_Kabsch", "separate_structures", "setup_MD", "setup_PB_GB", "setup_TI_MD", "TIW", "WAMM");

foreach my $addlib_file (@addlib_files){
	if((! -e $current_dir."/additional_libs/".$addlib_file.".pm")
	 &&(! -e $amber_home."/AmberTools/src/FEW/additional_libs/".$addlib_file.".pm")
	 && !(eval 'require '.$addlib_file.' ; 1;')) {
		print "\nThe Perl-module '$addlib_file' required for execution of FEW was not found.\n";
		print "Please ensure that the folder 'additional_libs' is located in the same folder\n";
		print "as the FEW.pl script.\n\n";
		exit;
	}
}

foreach my $pm_file (@pm_files){
	if(!(-e $current_dir."/libs/".$pm_file.".pm" ) && !(eval 'require '.$pm_file.' ; 1;')){
		print "\nThe Perl-module '$pm_file' required for execution of FEW was not found.\n";
		print "Please ensure that the folder 'libs' is located in the same folder\n";
		print "as the FEW.pl script.\n\n";
		exit;
	}
}

#################################################
# Check existence of required external programs #
#################################################

my @program_calls = ("antechamber", "tleap", "ambpdb", "cpptraj", "parmchk2");
foreach my $program_call (@program_calls){
	$call = "which $program_call >log.txt 2>error.txt";
	system($call);
	
	if(-z "log.txt"){
		print "\nERROR: The program '$program_call' cannot be invoked. Please ensure that\n";
		print "all programs of AmberTools can be invoked by their name.\n";
		print "In case of a bash shell this can e.g. be done by specifying\n\n";
		print "PATH=\$AMBERHOME/bin:\"\${PATH}\"\n\n";
		print "in the .bashrc file.\n\n";
		unlink("log.txt");
		unlink("error.txt");
		exit;
	}
	else{
		unlink("log.txt");
		unlink("error.txt");
	}
}


my $libs_path;
my $additional_libs_path;
BEGIN {
	my $amber_home = $ENV{"AMBERHOME"};
	my $called_few_dir = $FindBin::Bin;
	# In case link in $AMBERHOME/bin is called
	my $new_few_dir = "";
	if($called_few_dir =~ m/bin/){
		$new_few_dir .= $amber_home."/AmberTools/src/FEW";
	}
	else{
		$new_few_dir = $called_few_dir;
	}
	$libs_path = $new_few_dir."/libs";
	$additional_libs_path = $new_few_dir."/additional_libs";
	
	my $perl_five_lib_path = $amber_home."/AmberTools/src/FEW/additional_libs/PerlMol";
}
use lib $libs_path;
use lib $additional_libs_path;


use strict;
require read_data;
require global;
require check_input;
require WAMM;
require LIEW;
require TIW;

print "\n\t\t#####################################################\n";
print   "\t\t#                                                   #\n";
print   "\t\t#                        FEW                        #\n";
print   "\t\t#               Free Energy Workflow                #\n";
print   "\t\t#                                                   #\n";
print   "\t\t#####################################################\n";
print   "\n";

my $calc_procedure = $ARGV[0];
my $c_file = $ARGV[1];

if((!$calc_procedure)||(!$c_file)){
	print "Usage: perl FEW.pl <calculation procedure> <command file>\n\n";
	exit;
}
if(! -e $c_file){
	print "The specified command file does not exist.\n\n";
	exit;
}

#####################
# Read command file #
#####################

my ($c_header, $c_in_ref, $traj_files_ref, $n_traj_inter, $traj_inter_ref) = read_data($c_file);
my %c_in = %{$c_in_ref};
my @traj_files = @{$traj_files_ref};
my %traj_inter = %{$traj_inter_ref};

# Associate new to old keywords
my $c_in_ref = keyword_association(\%c_in);
my %c_in = %{$c_in_ref};

#################
#  Check input  #
#################

# Check consistency between command file and requested procedure
if((($calc_procedure eq "MMPBSA")||($calc_procedure eq "MMGBSA"))&&(!($c_header =~ m/\@WAMM/))){
	print "The provided command-file is not the file required for the requested procedure.\n";
	print "Please provide a file containing the commands for the Workflow for automated MM-PBSA\n";
	print "and MM-GBSA calculations (WAMM).\n";
	exit;
}
if(($calc_procedure eq "TI")&&(!($c_header =~ /\@TIW/))){
	print "The provided command-file is not the file required for the requested procedure.\n";
	print "Please provide a file containing the commands for the thermodynamic integration workflow (TIW)\n";
	exit;
}
if(($calc_procedure eq "LIE")&&(!($c_header =~ /\@LIEW/))){
	print "The provided command-file is not the file required for the requested procedure.\n";
	print "Please provide a file containing the commands for the linear interaction energy workflow (LIEW)\n";
	exit;
}

# Set default values
my $c_in_ref = &set_defaults($current_dir, $calc_procedure, \%c_in);
%c_in = %{$c_in_ref};

# Common check
my $check_common = &check_input_common(\%c_in);

# End program if check failed
if($check_common != 1){
	exit;
}


if(($calc_procedure eq "MMPBSA")||($calc_procedure eq "MMGBSA")){
	WAMM($current_dir, $n_traj_inter, \@traj_files, \%traj_inter, \%c_in);
}
elsif($calc_procedure eq "LIE"){
	LIEW($current_dir, $c_file, \@traj_files, \%c_in);
}
elsif($calc_procedure eq "TI"){
	TIW($current_dir, \%c_in);
}
else{
	print "\n\nThe procedure you specified in the program call does not exist.\n";
}
