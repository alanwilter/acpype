#!/usr/bin/perl
#
# matrix_test
#
# To a test of the conversion functions between Math::VectorReal
# and Math::MatrixReal functions by defining a generalise rotation
# matrix generator function.
#
# For details see
#     http://www.sct.gu.edu.au/~anthony/info/graphics/matrix.hints
#
# location of my libraries (I normally set this via PERL5LIB)
# #use lib '/home/anthony/lib/perl5';  # location of my libraries
use FindBin;
use lib "$FindBin::Bin/blib/lib";

use strict;
use Math::MatrixReal;
use Math::VectorReal qw( :all );

sub rotation_matrix {  # create rotation matrix to rotate first vector to second 
                      # second vector, through plane common to both.
  # Input vectors (normalise the input)
  my $v = shift->norm;
  my $w = shift->norm;

  # workout the axis framework around $v (to become X)
  my $n = ( $v x $w )->norm;  # axis of rotation (to become Z)
  my $m = $n x $v;            # othogonal to $v and $n (to become Y)

  print "Axis of Rotation (n) = $n\n\n";

  # Rotation matrix from this framework to the XYZ framework (inverted)
  # $v mapping to Z axis
  #
  # Creating a Matrix using strings is a crazy way to to do things!
  # $Q = Math::MatrixReal->new_from_string("[ $v ]\n[ $m ]\n[ $n ]\n");
  #
  # Math::VectorReal has a vector_matrix() is MUCH easier and more exact.
  my $Q = vector_matrix( $v, $m, $n );

  print "Matrix to rotate Z axis to n => Q\n$Q\n";

  # Transfrom  $w  to the new coordinate system
  # Note ~$Q is the transpose of $Q which is also its inverse
  my $wx =  $w * ~$Q;    # produce w'
  my $wy =  Z x $wx;    # and the third orthogonal vector.

  # Rotate around Z axis so   v' ($x)  -->  w' ($wx)
  my $T = vector_matrix( $wx, $wy, Z );
  print "Rotate around Z axis => T\n$T\n";

  # multiply together the three resulting rotation matrix
  # to produce the general rotation matrix to map $v to $w
  return  (~$Q) * $T * $Q;
}


my $v = vector( 1, 1, 1 );        # rotate this vectior to Y
#my $v = vector( shift, shift, shift );        # rotate this vectior to Y

$Math::VectorReal::FORMAT = "[ %g %g %g ]";
print "Generate a rotation matix to rotate $v to ", Y, "\n",
      "Vectors are normalised and rotation is by the minimal amount.\n\n";

$Math::VectorReal::FORMAT = "[ %5f %5f %5f ]";
my $R = rotation_matrix( $v, Y ); # generate the general rotation matrix
print "Resultant Rotation Matrix => R = transpose(Q) * T * Q\n$R\n";

my( $nx, $ny, $nz );

print "Rotate some vectors with matrix =>  v' =  v * R\n";
print "v -> ", $v * $R, "\n";   # test out the matrix generated
print "x -> ", $nx = X * $R, "\n";
print "y -> ", $ny = Y * $R, "\n";
print "z -> ", $nz = Z * $R, "\n";
print "\n";

print "Check orthogonaly is preserved\n";
print "'x x 'y -> ", $nx x $ny, " = 'z\n";
print "'y x 'z -> ", $ny x $nz, " = 'x\n";
print "'z x 'x -> ", $nz x $nx, " = 'y\n";
print "\n";

exit(0);

