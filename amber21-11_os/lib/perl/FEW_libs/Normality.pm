package Statistics::Normality;

use warnings;
use strict;
use Carp;
use Distributions;

=head1 NAME

Statistics::Normality - test whether an empirical distribution can be taken as being drawn from a normally-distributed population

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

use Statistics::Normality ':all';
use Statistics::Normality 'shapiro_wilk_test';
use Statistics::Normality 'dagostino_k_square_test';

=head1 DESCRIPTION

Various situations call for testing whether an empirical sample
can be presumed to have been drawn from a normally
(L<Gaussian|http://en.wikipedia.org/wiki/Normal_distribution>)
distributed population, especially because many downstream
significance tests depend upon the assumption of
normality.
This package implements some of the more
L<well-known tests|http://en.wikipedia.org/wiki/Normality_test>
from the mathematical statistics literature, though there
are also others that are not
included.
The tests here are all so-called I<omnibus> tests that
find departures from normality
on the basis of skewness and/or kurtosis [Dagostino71].
Note that, although the
L<Kolmogorov-Smirnov test|http://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test>
can also be used in this capacity, it is a
I<distance> test and therefore not advisable [Dagostino71].
This, and other distance tests (e.g. Chi-square) are not implemented
here.

=head1 TESTS

The subtleties and esoterica of various statistical tests for
normality require some familiarity with the mathematical statistics
literature.
We give rules-of-thumb for specific tests, where they exist, but it may be
advisable to try several different tests to check the consistency of the
conclusion.
It is probably also a good idea to check results graphically, either by
direct plotting or by a L<Q-Q plot|http://en.wikipedia.org/wiki/Q-Q_plot>.
In general, small samples will often pass a normality test suggesting
the I<possibility> that there is insufficient information to detect
departure from normal for such cases, should it
exist.

Each of the methods here is a frequentist test, i.e. one that tests against the
L<null-hypothesis|http://en.wikipedia.org/wiki/Null_hypothesis>
that the sample is
normal.
In other words, a low p-value recommends rejecting the
null.

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=cut

use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw/shapiro_wilk_test dagostino_k_square_test/;
%EXPORT_TAGS = (all => [qw/shapiro_wilk_test dagostino_k_square_test/]);
my $pkg = 'Statistics::Normality';


#  GENERAL IMPLEMENTATION NOTES FOR STATISTICAL TESTS
#
#  (1) Use standard Horner's Rule for polynomial evaluations, see e.g. Forsythe,
#      Malcolm, and Moler "Computer Methods for Mathematical Computations"
#      (1977) Prentice-Hall, pp 68.
#
#  (2) VAGARIES OF THE PERL Statistics::Distributions PACKAGE    symmetric
#      This package is implemented in the opposite way that one     :   |
#      finds tables of the standard normal function presented in   /:\  |
#      textbooks, where F(z) = A1 is the area from -infinity to : / : \ |
#      Z (or sometimes from 0 to Z). Instead, the Perl package  :/  :  \|
#      gives f(z) = A2 as the area from Z to +infinity, i.e.    /   :   \
#      in the *context of a significance test*. Note the       /:   :   |\
#      following implications for this package:               / :   A1  | \
#                                                            /  :   :   |A2\
#          f(Z) = F(-Z)       F(Z) + f(Z) + 1             __/___:___:___|___\___
#                                                              -Z   0   Z
#                         -1          -1              -1
#         udistr:    Z = f  [f(Z)] = f  [1 - F(Z)] = f  (1 - A1)
#
#       and the same appears to be true for other distributions in this package,
#       e.g. chi-square, student's T, etc.
#
#   (3) Tests that should perhaps be implemented in future versions:
#
#       * Anderson-Darling test
#       * Jarque-Bera test

#######################
#  SHAPIRO-WILK TEST  #
#######################

=head2 Shapiro-Wilk Test

The L<Shapiro-Wilk W-Statistic test|http://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test>
[Shapiro65] is considered to be among the most objective tests of
normality [Royston92] and also one of the most powerful ones for detecting
non-normality [Chen71].
Its statistic is essentially the roughly best unbiased estimator
of population standard deviation to the sample variance [Dagostino71].
The test is mathematically complex and most implementations use several
conventional approximations (as we do here), including Blom's formula
for the expected value of the order statistics [Harter61] and transformation
to standard normal distribution for evaluation, especially for large
samples [Royston92].

	$pval = shapiro_wilk_test ([0.34, -0.2, 0.8, ...]);
	($pval, $w_statistic) = shapiro_wilk_test ([0.34, -0.2, 0.8, ...]);

This test may not be the best if there are many repeated values in the test
distribution or when the number of points in the test distribution is very
large, e.g. more than 5000.
The routine will L<carp|Carp> about the latter, but not the
former.
This particular implementation of the test also requires at
least 6 data points in the sample distribution and will L<croak|Carp>
otherwise.

=cut

#  IMPLEMENTATION NOTES FOR SHAPIRO-WILK TEST
#
#  (1) THIS ROUTINE IS NON-TRIVIAL TO IMPLEMENT --- THE IMPLEMENTATION HERE
#      ROUGHLY FOLLOWS THE EXPLANATIONS GIVEN IN A PAPER BY GUNER AND
#      JOHNSON (GJ), EXCEPT WHERE THIS PAPER HAS ERRORS, AS NOTED IN THE
#      CODE. A DRAFT OF THIS PAPER IS AVAILABLE AT:
#
#             http://esl.eng.ohio-state.edu/~rstheory/iip/shapiro.pdf
#
#      BUT NOTE THAT WE CHANGE THEIR 1->N VECTOR INDEX NOTATION TO 0->N-1 TO BE
#      CONSISTENT WITH HOW PERL STORES LISTS
#
#  (2) Size limits for empirical distributions to be tested. The log-standard
#      normal transformation "Z-form" used here has been empirically shown to be
#      good up to 2000 data points according to GJ, but GJ also claims up to
#      5000 is OK. This version that uses the Blom approximation also requires
#      at least 6 points because the vector of m-values is anti-symmetric and
#      6 is the lowest number of members for which the numerator of epsilon
#      does not vanish:
#
#      num = M*M - 2[m_(n)^2 + m_(n-1)^2]
#          = [m_(1)^2 + m_(2)^2 + m_(3)^2 +...+ m_(n)^2 - 2[m_(n)^2 + m_(n-1)^2]
#
#        # points   notes
#        1          no distribution (a point)
#        2          no distribution (a line)
#        3          m1 & m3 cancelled by 2*m3 and m2=0
#        4          termwise cancellation: m1 & m4 by 2*m4 and m2 & m3 by 2*m3
#        5          similar except m3 identically 0
#
#  (3) test case from http://www.statsdirect.com/help/parametric_methods/swt.htm
#
#      the list qw/0.0987 0.0000 0.0533 -0.0026 0.0293 -0.0036 0.0246 -0.0042
#                  0.0200 -0.0114 0.0194 -0.0139 0.0191 -0.0222 0.0180 -0.0333
#                  0.0172 -0.0348 0.0132 -0.0363 0.0102 -0.0363 0.0084 -0.0402
#                  0.0077 -0.0583 0.0058 -0.1184 0.0016 -0.1420/
#
#      should give  $w_statistic \approx 0.89218 and $pval \approx 0.00544

sub shapiro_wilk_test {
   my ($sample_distribution) = @_;

#__SOME CONSTANTS
   my $three_eighths = 3/8;
   my $one_fourth = 1/4;

#__ARGUMENT MUST BE A REFERENCE TO A LIST OF NUMBERS (NOTE REGEXP::COMMON
#  DOES NOT SEEM TO PROVIDE A GENERAL REGEXP FOR FLOATING POINTS OF VARIOUS
#  FORMS (AS ITS TITLE WOULD SEEM TO INDICATE) --- FOUND A PRETTY GOOD REGEXP
#  AT http://perl.active-venture.com/pod/perlretut-regexp.html WHICH SEEMS TO
#  BE PART OF THE PERL REGULAR EXPRESSIONS TUTORIAL
   croak "arg must be list ref" unless ref $sample_distribution eq "ARRAY";
   foreach my $string (@{$sample_distribution}) {
      croak "'$string' not a number" unless
         $string =~ /^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$/;
   }
 
#__SORT AND SIZE
   @{$sample_distribution} = sort _numerical_ @{$sample_distribution};
   my $num_vals = scalar @{$sample_distribution};

#__LIST SIZE HAS CERTAIN LIMITS AS NOTED ABOVE
   carp "results not guaranteed for >5000 points" if $num_vals > 5000;
   croak "sample must have at least 6 points" unless $num_vals >= 6;

#__VECTOR OF M-VALUES (GJ EQN 7) I.E. EXPECTED VALUES OF THE ORDER STATISTICS
#  OF THE STANDARD NORMAL DISTRIBUTION USING BLOM'S APPROXIMATION -- SEE
#  [Harter61] PP 153 AND NOTE INDEXING CORRECTIONS AND IMPLEMENTATION NOTES
#  FOR "UDISTR"
#
#  CERTIFICATION: THIS BLOCK RETURNS VALUES THAT ARE CONSISTENT WITH
#                 [Harter61] TABLE 1 PP 158-165, THOUGH NOT EXACT BECAUSE
#                 OF THE APPROXIMATION, WHERE THEY LIST THE POSITIVE VALUES
#                 CORRESPONDING TO THE I-TH *LARGEST* NORMAL DEVIATE --- NOTE
#                 THAT THIS LOOP IS BASED ON ORDER FROM LEAST TO GREATEST
#                 SO WE CALCULATE THE I-TH *SMALLEST* NORMAL DEVIATES FIRST
#
#                 THE RESULT IS ANTI-SYMMETRIC
   my @mvals = ();
   my $mean = 0;
   for (my $i = 0; $i < $num_vals; $i++) {
      my $index = $i + 1;
      my $arg_p = ($index - $three_eighths)/($num_vals + $one_fourth);
      my $mval = Statistics::Distributions::udistr (1 - $arg_p);
      push @mvals, $mval;
      $mean += $sample_distribution->[$i];
   }
   $mean /= $num_vals;

#__DOT-PRODUCT OF VECTOR M WITH ITSELF
   my $m_dot_m = 0;
   foreach my $mval (@mvals) {
      $m_dot_m += $mval * $mval;
   }

#__VECTOR OF C-VALUES (GJ EQN 9) --- NO $INDEX CORRECTION NECESSARY HERE
   my $sqrt_m_dot_m = sqrt ($m_dot_m);
   my @cvals = ();
   for (my $i = 0; $i < $num_vals; $i++) {
      my $cval = $mvals[$i] / $sqrt_m_dot_m;
      push @cvals, $cval;
   }

#__CALCULATE THE LAST 2 MEMBERS OF THE VECTOR OF A-VALUES USING POLYNOMIAL
#  APPROXIMATION [Royston92] EVALUATED USING HORNER'S RULE (GJ EQNS 4 AND 5)
#__NOTE INDEXING CORRECTIONS
   my $position_n = $num_vals - 1;
   my $uval = 1 / sqrt($num_vals);
   my $a_n     = ((((-2.706056*$uval+4.434685)*$uval-2.07119)*$uval-0.147981)*$uval+0.221157)*$uval+$cvals[$position_n];
   my $a_n_m_1 = ((((-3.582633*$uval+5.682633)*$uval-1.752461)*$uval-0.293762)*$uval+0.042981)*$uval+$cvals[$position_n - 1];

#__EPSILON (GJ EQN 8) --- NOTE INDEXING CORRECTIONS
   my $epsilon = $m_dot_m - 2 * (
         $mvals[$position_n] * $mvals[$position_n]
         +
         $mvals[$position_n - 1] * $mvals[$position_n - 1]
      );
   $epsilon /= 1 - 2 * ($a_n * $a_n + $a_n_m_1 * $a_n_m_1);
   my $sqrt_epsilon = sqrt($epsilon);

#__BUILD THE VECTOR OF A-VALUES (GJ EQNS 4-6) --- NOTE INDEXING CORRECTIONS
   my @avals = ();
   push @avals, - $a_n;
   push @avals, - $a_n_m_1;
   for (my $i = 2; $i <= $num_vals - 3; $i++) {
      push @avals, $mvals[$i] / $sqrt_epsilon;
   }
   push @avals, $a_n_m_1;
   push @avals, $a_n;

#__SHAPIRO-WILK W-STATISTIC --- NOTE GJ EQN 1 FOR THIS ENTITY IS INCORRECT, AS
#  THE ACTUAL NUMERATOR IS SQUARED AND THE ACTUAL DIFFERENCE IN THE DENOMINATOR
#  IS BETWEEN THE I-TH ELEMENT OF THE TEST (SAMPLE) DISTRIBUTION AND ITS AVERAGE
   my ($sw_numerator, $sw_denominator) = (0, 0);
   for (my $i = 0; $i < $num_vals; $i++) {
      $sw_numerator += $avals[$i] * $sample_distribution->[$i];
      $sw_denominator += ($sample_distribution->[$i] - $mean)**2;
   }
   $sw_numerator *= $sw_numerator;
   my $w_statistic = $sw_numerator / $sw_denominator;

#__TRANSFORMATION OF W-STATISTIC INTO STANDARD Z-STATISTIC USING POLYNOMIAL
#  APPROXIMATION [Royston92] EVAL'D USING HORNER'S RULE (GJ EQNS 10-12)
#  SO THAT WE CAN TEST USING THE STANDARD NORMAL DISTRIBUTION
   my $log_n = log ($num_vals);
   my $mu_z = ((0.0038915*$log_n - 0.083751)*$log_n - 0.31082)*$log_n - 1.5861;
   my $log_sigma_z = (0.0030302 * $log_n - 0.082676) * $log_n - 0.4803;
   my $sigma_z = exp ($log_sigma_z);
   my $z_statistic = (log (1 - $w_statistic) - $mu_z) / $sigma_z;

#__PVALUE FROM STANDARD NORMAL DISTRIBUTION -- SEE IMPLEMENTATION NOTES ABOVE:
#  "UPROB" GIVES SIGNIFICANCE PVALUE (I.E. AREA TO THE RIGHT OF THE STATISTIC)
   my $pval = Statistics::Distributions::uprob ($z_statistic);

#__RETURN RESULTS
   if (wantarray) {
      return ($pval, $w_statistic);
   } else {
      return $pval;
   }
}

##############################
#  D'AGOSTINO K-SQUARE TEST  #
##############################

=head2 D'Agostino K-Squared Test

The L<D'Agostino K-Squared test|http://en.wikipedia.org/wiki/D%27Agostino%27s_K-squared_test> is a good test against non-normality arising from
L<kurtosis|http://en.wikipedia.org/wiki/Kurtosis> and/or
L<skewness|http://en.wikipedia.org/wiki/Skewness> [Dagostino90].

	$pval = dagostino_k_square_test ([0.34, -0.2, ...]);
	($pval, $ksq_statistic) = dagostino_k_square_test ([0.34, -0.2, ...]);

The test statistic depends upon both the sample kurtosis and skewness, as
well as the moments of these parameters from a normal population, as quantified
by Pearson's coefficients [Pearson31].
These are transformed [Dagostino70,Anscombe83] to expressions that sum
to the K-squared statistic, which is essentially chi-square-distributed
with 2 degrees of
freedom [Dagostino90].
The kurtosis transform, and thus the overall test, generally works best
when the sample distribution has at least 20 data points [Anscombe83] and the
routine will L<carp|Carp>
otherwise.

=cut

#  IMPLEMENTATION NOTES FOR D'AGOSTINO K-SQUARED TEST
#
#  (1) WE ROUGHLY FOLLOW THE VARIABLE NOTATION GIVEN BY THE WIKIPEDIA ARTICLE
#      AT http://en.wikipedia.org/wiki/D'Agostino's_K-squared_test
#
#  (2) test case from [Dagostino90] PP 318
#
#      the list qw/393 353 334 336 327 300 300 308 283 285 270 270 272 278 278
#                  263 264 267 267 267 268 254 254 254 256 256 258 240 243 246
#                  247 248 230 230 230 230 231 232 232 232 234 234 236 236 238
#                  220 225 225 226 210 211 212 215 216 217 218 200 202 192 198
#                  184 167/;
#
#      should give  $ksq_statistic \approx 14.752 and $pval \approx 0.00063

sub dagostino_k_square_test {
   my ($sample_distribution) = @_;

#__ARGUMENT MUST BE A REFERENCE TO A LIST OF NUMBERS (NOTE REGEXP::COMMON
#  DOES NOT SEEM TO PROVIDE A GENERAL REGEXP FOR FLOATING POINTS OF VARIOUS
#  FORMS (AS ITS TITLE WOULD SEEM TO INDICATE) --- FOUND A PRETTY GOOD REGEXP
#  AT http://perl.active-venture.com/pod/perlretut-regexp.html WHICH SEEMS TO
#  BE PART OF THE PERL REGULAR EXPRESSIONS TUTORIAL
   croak "arg must be list ref" unless ref $sample_distribution eq "ARRAY";
   foreach my $string (@{$sample_distribution}) {
      croak "'$string' not a number" unless
         $string =~ /^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$/;
   }

#__MOMENTS OF SAMPLE DISTRIBUTION
   my ($n, $mean, $g_1, $g_2) = _sample_mean_skew_kurt_ ($sample_distribution);
   carp "results not guaranteed for <20 points" if $n < 20;

#__SELECTED PEARSON'S COEFFICIENTS [Pearson31] WITH HORNER'S RULE FOR POLYNOMS
   my $mu_2_g_1 = 6*($n-2) / (($n+1)*($n+3));
   my $gamma_2_g_1 = 36*($n-7)*(($n+2)*$n-5) / (($n-2)*($n+5)*($n+7)*($n+9));

   my $mu_1_g_2 = - 6 / ($n+1);
   my $mu_2_g_2 = 24*$n*($n-2)*($n-3) / (($n+1)*($n+1)*($n+3)*($n+5));
   my $gamma_1_g_2 = 6*(($n-5)*$n+2) / (($n+7)*($n+9));
   $gamma_1_g_2 *= sqrt(6*($n+3)*($n+5));
   $gamma_1_g_2 /= sqrt($n*($n-2)*($n-3));

#__TRANSFORMED SKEWNESS PARAMETER [Dagostino70]
   my $w_squared = sqrt(2*$gamma_2_g_1+4) - 1;
   my $delta = 1 / sqrt(log(sqrt($w_squared)));
   my $alpha_squared = 2/($w_squared - 1);
   my $ratio_squared = $g_1 * $g_1 / ($alpha_squared * $mu_2_g_1);
   my $z1 = $delta*log(sqrt($ratio_squared) + sqrt($ratio_squared + 1));

#__EQUIVALENT-ALTERNATE DERIVATION OF SKEWNESS TRANSFORM [Dagostino90]
#  my $y = $g_1 * sqrt(($n + 1)*($n + 3)/(6*($n - 2)));
#  my $beta_2 = 3*(($n + 27)*$n - 70)*($n + 1)*($n + 3);
#  $beta_2 /= ($n - 2)*($n + 5)*($n + 7)*($n + 9);
#  my $w_squared = sqrt(2*$beta_2 - 2) - 1;
#  my $delta = 1 / sqrt(log(sqrt($w_squared)));
#  my $alpha = sqrt(2/($w_squared - 1));
#  my $z1 = $delta*log($y/$alpha + sqrt(($y/$alpha)*($y/$alpha) + 1));

#__TRANSFORMED KURTOSIS PARAMETER [Anscombe83]
   my $a = sqrt(1 + 4 / ($gamma_1_g_2 * $gamma_1_g_2)) + 2/$gamma_1_g_2;
   $a *= 8/$gamma_1_g_2;
   $a += 6;
   my $term = 1 - 2/$a;
   $term /= 1 + ($g_2 - $mu_1_g_2) * sqrt(2/($a-4)) / sqrt($mu_2_g_2);
   my $z2 = 1 - 2/(9*$a) - $term**(1/3);
   $z2 *= sqrt(9*$a/2);

#__OMNIBUS K-SQUARED STATISTIC WHICH IS DSITRBUTED ROUGHLY CHI-SQUARED
#  WITH 2 DEGREES OF FREEDOM [Dagostino90]
   my $k_squared_stat = $z1*$z1 + $z2*$z2;
   my $pval = Statistics::Distributions::chisqrprob (2, $k_squared_stat);

#__RETURN RESULTS
   if (wantarray) {
      return ($pval, $k_squared_stat);
   } else {
      return $pval;
   }
}

sub _numerical_ {$a <=> $b}

sub _sample_mean_skew_kurt_ {
   my ($sample_distribution) = @_;

#__MEAN
   my $mean = 0;
   my $num_vals = scalar @{$sample_distribution};
   foreach my $datum (@{$sample_distribution}) {
      $mean += $datum;
   }
   $mean /= $num_vals;

#__SAMPLE SKEWNESS AND KURTOSIS
   my ($sum_square_diffs, $sum_cube_diffs, $sum_quad_diffs) = (0, 0, 0);
   foreach my $datum (@{$sample_distribution}) {
      my $sq_diff = ($datum - $mean) * ($datum - $mean);
      $sum_square_diffs += $sq_diff;
      $sum_cube_diffs += $sq_diff * ($datum - $mean);
      $sum_quad_diffs += $sq_diff * $sq_diff;
   }
   my $sum_square_diffs_over_n = $sum_square_diffs / $num_vals;
   my $g_1 = $sum_cube_diffs / $num_vals;
   $g_1 /= $sum_square_diffs_over_n**1.5;
   my $g_2 = $sum_quad_diffs / $num_vals;
   $g_2 /= $sum_square_diffs_over_n**2;
   $g_2 -= 3;

#__RESULTS: NUMBER OF VALUES, MEAN, SKEWNESS, KURTOSIS
   return ($num_vals, $mean, $g_1, $g_2);
}

=head1 REFERENCES

=over

=item *

[Anscombe83] Anscombe, F. J. and Glynn, W. J. (1983)
I<Distribution of the Kurtosis Statistic B2 for Normal Samples>,
Biometrika B<70>(1), 227-234.

=item *

[Chen71] Chen, E. H. (1971)
I<The Power of the Shapiro-Wilk W Test for Normality in Samples from Contaminated Normal Distributions>,
Journal of the American Statistical Association B<66>(336), 760-762.

=item *

[Dagostino70] D'Agostino, R. B. (1970)
I<Transformation to Normality of the Null Distribution of G1>,
Biometrika B<57>(3), 679-681.

=item *

[Dagostino71] D'Agostino, R. B. (1971)
I<An Omnibus Test of Normality for Moderate and Large Size Samples>,
Biometrika B<58>(2), 341-348.

=item *

[Dagostino90] D'Agostino, R. B. et al. (1990)
I<A Suggestion for Using Powerful and Informative Tests of Normality>,
American Statistician B<44>(4), 316-321.

=item *

[Harter61] Harter, H. L. (1961)
I<Expected values of normal order statistics>,
Biometrika B<48>(1/2), 151-165.

=item *

[Pearson31] Pearson, E. S. (1931)
I<Notes on Tests for Normality>,
Biometrika B<22>(3/4), 423-424.

=item *

[Royston92] Royston, J. P. (1992)
I<Approximating the Shapiro-Wilk W-test for non-normality>,
Statistics and Computing B<2>(3) 117-119.

=item *

[Shapiro65] Shapiro, S. S. and Wilk, M. B. (1965)
I<An analysis of variance test for normality - complete samp1es>,
Biometrika B<52>(3/4), 591-611.

=back

=head1 AUTHOR

Mike Wendl, C<< <mwendl at genome.wustl.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-statistics-normality at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Statistics-Normality>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Statistics::Normality


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Statistics-Normality>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Statistics-Normality>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Statistics-Normality>

=item * Search CPAN

L<http://search.cpan.org/dist/Statistics-Normality/>

=back

=head1 COPYRIGHT & LICENSE

Copyright (C) 2011 Washington University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

=cut

1; # End of Statistics::Normality
