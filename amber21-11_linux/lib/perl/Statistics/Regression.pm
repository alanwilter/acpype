package Statistics::Regression;

$VERSION = '0.15';

use strict;

################################################################
use constant TINY => 1e-8;
################################################################

=head1 NAME

  Regression.pm - weighted linear regression package (line+plane fitting)

=head1 SYNOPSIS

  use Statistics::Regression;

  # Create regression object
  my $reg = Statistics::Regression->new( 
    3, "sample regression", 
    [ "const", "someX", "someY" ] 
  );

  # Add data points
  $reg->include( 2.0, [ 1.0, 3.0, -1.0 ] );
  $reg->include( 1.0, [ 1.0, 5.0, 2.0 ] );
  $reg->include( 20.0, [ 1.0, 31.0, 0.0 ] );
  $reg->include( 15.0, [ 1.0, 11.0, 2.0 ] );

  # Print the result
  $reg->print(); 

  # Prints the following:
  # ****************************************************************
  # Regression 'sample regression'
  # ****************************************************************
  # Theta[0='const']=       0.2950
  # Theta[1='someX']=       0.6723
  # Theta[2='someY']=       1.0688
  # R^2= 0.808, N= 4
  # ****************************************************************

  # Or, to get the values of the coefficients and R^2
  my @theta = $reg->theta;
  my $rsq   = $reg->rsq;

=head1 DESCRIPTION

Regression.pm is a multivariate linear regression package.  That is, it
estimates the c coefficients for a line-fit of the type

y= c(0)*x(0) + c(1)*x1 + c(2)*x2 + ... + c(k)*xk

given a data set of N observations, each with k independent x variables and one
y variable.  Naturally, N must be greater than k---and preferably considerably
greater.  Any reasonable undergraduate statistics book will explain what a
regression is.  Most of the time, the user will provide a constant ('1') as
x(0) for each observation in order to allow the regression package to fit an
intercept.

=head1 ALGORITHM

=head2 Original Algorithm (ALGOL-60):

	W.  M.  Gentleman, University of Waterloo, "Basic Description
	For Large, Sparse Or Weighted Linear Least Squares Problems
	(Algorithm AS 75)," Applied Statistics (1974) Vol 23; No. 3

=head2 INTERNALS

R=Rbar is an upperright triangular matrix, kept in normalized
form with implicit 1's on the diagonal.  D is a diagonal scaling
matrix.  These correspond to "standard Regression usage" as

                X' X  = R' D R

A backsubsitution routine (in thetacov) allows to invert the R
matrix (the inverse is upper-right triangular, too!). Call this
matrix H, that is H=R^(-1).

	  (X' X)^(-1) = [(R' D^(1/2)') (D^(1/2) R)]^(-1)
	  = [ R^-1 D^(-1/2) ] [ R^-1 D^(-1/2) ]'

=head2 Remarks

This algorithm is the statistical "standard." Insertion of a new observation
can be done one obs at any time (WITH A WEIGHT!), and still only takes a low
quadratic time.  The storage space requirement is of quadratic order (in the
indep variables). A practically infinite number of observations can easily be
processed!


=head1 METHODS

=cut

################################################################


#### let's start with handling of missing data ("nan" or "NaN")

my $nan= "NaN";
sub isNaN { 
  if ($_[0] !~ /[0-9nan]/) { die "definitely not a number in NaN: '$_[0]'"; }
  return ($_[0]=~ /NaN/i) || ($_[0] != $_[0]);
}


################################################################

=head2 new

 my $reg = Statistics::Regression->new($n, $name, \@var_names)

Receives the number of variables on each observations (i.e., an integer) and
returns the blessed data structure as a Statistics::Regression object. Also
takes an optional name for this regression to remember, as well as a reference
to a k*1 array of names for the X coefficients.

=cut

################################################################
sub new {
  my $classname= shift(@_);
  my $K= shift(@_); # the number of variables
  my $regname= shift(@_) || "with no name";

  if (!defined($K)) { die "Regression->new needs at least one argument for the number of variables"; }
  if ($K<=1) { die "Cannot run a regression without at least two variables."; }

  sub zerovec {
    my @rv;
    for (my $i=0; $i<=$_[0]; ++$i) { $rv[$i]=0; } 
    return \@rv;
  }

  bless {
	 k => $K,
	 regname => $regname,
	 xnames => shift(@_),

	 # constantly updated
	 n => 0,
	 sse => 0,
	 syy => 0,
	 sy => 0,
	 wghtn => 0,
	 d => zerovec($K),
	 thetabar => zerovec($K),
	 rbarsize => ($K+1)*$K/2+1,
	 rbar => zerovec(($K+1)*$K/2+1),

	 # other constants
	 neverabort => 0,

	 # computed on demand
	 theta => undef,
	 sigmasq => undef,
	 rsq => undef,
	 adjrsq => undef
	}, $classname;
}

################################################################

=head2 dump

  $reg->dump

Used for debugging.

=cut

################################################################
sub dump {
  my $this= $_[0];
  print "****************************************************************\n";
  print "Regression '$this->{regname}'\n";
  print "****************************************************************\n";
  sub print1val {
    no strict;
    print "$_[1]($_[2])=\t". ((defined($_[0]->{ $_[2] }) ? $_[0]->{ $_[2] } : "intentionally undef"));

    my $ref=$_[0]->{ $_[2] };

    if (ref($ref) eq 'ARRAY') {
      my $arrayref= $ref;
      print " $#$arrayref+1 elements:\n";
      if ($#$arrayref>30) {
	print "\t";
	for(my $i=0; $i<$#$arrayref+1; ++$i) { print "$i='$arrayref->[$i]';"; }
	print "\n";
      }
      else {
	for(my $i=0; $i<$#$arrayref+1; ++$i) { print "\t$i=\t'$arrayref->[$i]'\n"; }
      }
    }
    elsif (ref($ref) eq 'HASH') {
      my $hashref= $ref;
      print " ".scalar(keys(%$hashref))." elements\n";
      while (my ($key, $val) = each(%$hashref)) {
	print "\t'$key'=>'$val';\n";
      }
    }
    else {
      print " [was scalar]\n"; }
  }

  while (my ($key, $val) = each(%$this)) {
    $this->print1val($key, $key);
  }
  print "****************************************************************\n";
}

################################################################

=head2 print

  $reg->print

prints the estimated coefficients, and R^2 and N. For an example see the
SYNOPSIS.

=cut

################################################################
sub print {
  my $this= $_[0];
  print "****************************************************************\n";
  print "Regression '$this->{regname}'\n";
  print "****************************************************************\n";

  my $theta= $this->theta();

  for (my $i=0; $i< $this->k(); ++$i) {
    print "Theta[$i".(defined($this->{xnames}->[$i]) ? "='$this->{xnames}->[$i]'":"")."]= ".sprintf("%12.4f", $theta->[$i])."\n";
  }
  print "R^2= ".sprintf("%.3f", $this->rsq()).", N= ".$this->n()."\n";
  print "****************************************************************\n";
}


################################################################

=head2 include

  $n = $reg->include( $y, [ $x1, $x2, $x3 ... $xk ], $weight );

Add one new observation. The weight is optional. Note that inclusion with a
weight of -1 can be used to delete an observation.

Returns the number of observations so far included.

=cut

################################################################
sub include {
  my $this = shift();
  my $yelement= shift();
  my $xrow= shift();
  my $weight= shift() || 1.0;

  # omit observations with missing observations;
  if (!defined($yelement)) { die "Internal Error: yelement is undef"; }
  if (isNaN($yelement)) { return $this->{n}; }

  my @xcopy;
  for (my $i=1; $i<=$this->{k}; ++$i) { 
    if (!defined($xrow->[$i-1])) { die "Internal Error: xrow [ $i-1 ] is undef"; }
    if (isNaN($xrow->[$i-1])) { return $this->{n}; }
    $xcopy[$i]= $xrow->[$i-1];
  }

  $this->{syy}+= ($weight*($yelement*$yelement));
  $this->{sy}+= ($weight*($yelement));
  if ($weight>=0.0) { ++$this->{n}; } else { --$this->{n}; }

  $this->{wghtn}+= $weight;

  for (my $i=1; $i<=$this->{k};++$i) {
    if ($weight==0.0) { return $this->{n}; }
    if (abs($xcopy[$i])>(TINY)) {
      my $xi=$xcopy[$i];

      my $di=$this->{d}->[$i];
      my $dprimei=$di+$weight*($xi*$xi);
      my $cbar= $di/$dprimei;
      my $sbar= $weight*$xi/$dprimei;
      $weight*=($cbar);
      $this->{d}->[$i]=$dprimei;
      my $nextr=int( (($i-1)*( (2.0*$this->{k}-$i))/2.0+1) );
      if (!($nextr<=$this->{rbarsize}) ) { die "Internal Error 2"; }
      my $xk;
      for (my $kc=$i+1;$kc<=$this->{k};++$kc) {
	$xk=$xcopy[$kc]; $xcopy[$kc]=$xk-$xi*$this->{rbar}->[$nextr];
	$this->{rbar}->[$nextr]= $cbar * $this->{rbar}->[$nextr]+$sbar*$xk;
	++$nextr;
      }
      $xk=$yelement; $yelement-= $xi*$this->{thetabar}->[$i];
      $this->{thetabar}->[$i]= $cbar*$this->{thetabar}->[$i]+$sbar*$xk;
    }
  }
  $this->{sse}+=$weight*($yelement*$yelement);

  # indicate that Theta is garbage now
  $this->{theta}= undef;
  $this->{sigmasq}= undef; $this->{rsq}= undef; $this->{adjrsq}= undef;

  return $this->{n};
}



################################################################


=head2 theta

  $theta = $reg->theta
  @theta = $reg->theta

Estimates and returns the vector of coefficients. In scalar context returns an
array reference; in list context it returns the list of coefficients.

=cut

################################################################

sub theta {
  my $this= shift();

  if (defined($this->{theta})) { 
    return wantarray ? @{$this->{theta}} : $this->{theta}; 
  }

  if ($this->{n} < $this->{k}) { return; }
  for (my $i=($this->{k}); $i>=1; --$i) {
    $this->{theta}->[$i]= $this->{thetabar}->[$i];
    my $nextr= int (($i-1)*((2.0*$this->{k}-$i))/2.0+1);
    if (!($nextr<=$this->{rbarsize})) { die "Internal Error 3"; }
    for (my $kc=$i+1;$kc<=$this->{k};++$kc) {
      $this->{theta}->[$i]-=($this->{rbar}->[$nextr]*$this->{theta}->[$kc]);
      ++$nextr;
    }
  }

  my $ref = $this->{theta}; shift(@$ref); # we are counting from 0

  # if in a scalar context, otherwise please return the array directly
  wantarray ? @{$this->{theta}} : $this->{theta};
}

################################################################

=head2 rsq, adjrsq, sigmasq, ybar, sst, k, n

  $rsq = $reg->rsq; # etc...

These methods provide common auxiliary information.  rsq, adjrsq,
sigmasq, sst, and ybar have not been checked but are likely correct.
The results are stored for later usage, although this is somewhat
unnecessary because the computation is so simple anyway.

=cut

################################################################

sub rsq {
  my $this= shift();
  return $this->{rsq}= 1.0- $this->{sse} / $this->sst();
}

sub adjrsq {
  my $this= shift();
  return $this->{adjrsq}= 1.0- (1.0- $this->rsq())*($this->{n}-1)/($this->{n} - $this->{k});
}

sub sigmasq {
  my $this= shift();
  return $this->{sigmasq}= ($this->{n}<=$this->{k}) ? "Inf" : ($this->{sse}/($this->{n} - $this->{k}));
}

sub ybar {
  my $this= shift();
  return $this->{ybar}= $this->{sy}/$this->{wghtn};
}

sub sst {
  my $this= shift();
  return $this->{sst}= ($this->{syy} - $this->{wghtn}*($this->ybar())**2);
}

sub k {
  my $this= shift();
  return $this->{k};
}
sub n {
  my $this= shift();
  return $this->{n};
}


################################################################

=head1 BUGS/PROBLEMS

=over 4

=item Missing

This package lacks routines to compute the standard errors of
the coefficients.  This requires access to a matrix inversion
package, and I do not have one at my disposal.  If you want to
add one, please let me know.

=item Perl Problem

perl is unaware of IEEE number representations.  This makes it a
pain to test whether an observation contains any missing
variables (coded as 'NaN' in Regression.pm).

=back

=for comment
pod2html -noindex -title "perl weighted least squares regression package" Regression.pm > Regression.html

=head1 VERSION

0.15

=head1 AUTHOR

Naturally, Gentleman invented this algorithm.  Adaptation by ivo welch. Alan
Miller (alan@dmsmelb.mel.dms.CSIRO.AU) pointed out nicer ways to compute the
R^2. Ivan Tubert-Brohman helped wrap the module as as standard CPAN
distribution.

=head1 LICENSE

This module is released for free public use under a GPL license.

(C) Ivo Welch, 2001,2004.

=cut

1;

# vim: sw=2 sts=2
