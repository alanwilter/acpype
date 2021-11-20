#
# Module to create coordinates for mm_pbsa
#
# Holger Gohlke: 18.10.2001
#
# Last update: 14.02.2014
#
########################################################################

package mm_pbsa_createcoords;
require Exporter;

@ISA         = qw(Exporter);
@EXPORT      = qw(create_coords);
@EXPORT_OK   = qw();
$VERSION     = 1.00;

########################################################################

use strict;
use mm_pbsa_global qw($HTMLPATH);

########################################################################

sub create_coords(){
####################

  # Parameters: \%PRO,\$TRA
  my $r_pro = shift;
  my $r_tra = shift;

  print "\n";
  print "=>> Creating coordinates\n";
    
  my $cat;
  if($$r_tra =~ /\.gz/){
    $cat = "gzip -cd";
  }
  elsif($$r_tra =~ /\.Z/){
    die("    Specify path to zcat in subroutine create_coords in mm_pbsa_createcoords.pm\n    For details see: $HTMLPATH#spec_path_to_zcat");
    #$cat = ... 
  }
  else{
    $cat = "/bin/cat";
  }

  # For more than one trajectory file, remove title line for the 2nd, 3rd, ...
  # For this to work, all title lines must be the same!
  my @tmp;
  my $no = scalar(@tmp = split/ +/,$$r_tra);
  print "    Executing makecrd\n";
  open(FROM,"$cat $$r_tra |")           || die("    Open FROM failed: $!\n    For details see: $HTMLPATH#makecrd_open_from_failed");
  open(TO  ,"| $r_pro->{\"MAKE_CRD\"}") || die("    Open TO failed: $!\n    For details see: $HTMLPATH#makecrd_open_to_failed");
  my $line;
  my $title = "";
  while(defined($line = <FROM>)){
    if($title eq ""){
      $title = $line;
    }
    elsif("$line" eq "$title"){
      print "        Skipped title_line\n";
      next;
    }
    print TO $line;
  }
  close(TO)   || die("    Close TO failed: $?/$!\n    For details see: $HTMLPATH#makecrd_close_from_failed");
  close(FROM) || die("    Close FROM failed: $?/$!\n    For details see: $HTMLPATH#makecrd_close_to_failed");
}

1; # Necessary for package function
