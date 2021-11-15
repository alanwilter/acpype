#!/usr/bin/perl
#
# Direct extraction of the VectorReal POD synopsis
#
# location of my libraries (I normally set this via PERL5LIB)
# # #use lib '/home/anthony/lib/perl5';  # location of my libraries
use FindBin;
use lib "$FindBin::Bin/blib/lib";

use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors

$a = vector( 1, 2, .5 );
print "Vector as string (default ouput format)\n\$a => ", $a;
print "\n";

print  "Specified output format\n";
print  $a->stringify("\$a => { %g, %g, %g }\n");
print "\n";

# I hate newline in the default output format (defined as MatrixReal)
$Math::VectorReal::FORMAT = "[ %.5f %.5f %.5f ]";
print "Modified default output format\n\$a => $a\n";
print "\n";

print "General Vector Mathematics\n";
print "length\n",     $a->length, "\n";
print "normalised\n", $a->norm,   "\n";
print 'string concat   $a."**" = ', $a."**", "\n";
print 'vector constant    X    = ',   X,    "\n";
print 'subtraction     $a - Z  = ', $a - Z, "\n";
print 'scalar divide   $a / 3  = ', $a / 3, "\n";
print 'dot product     $a . Y  = ', $a . Y, "\n";
print 'cross product   $a x Y  = ', $a x Y, "\n";
print "\n";


print "Plane containing points X, \$a, Z (in anti-clockwise order)\n";
($n,$d) = plane( X, $a, Z );  # return normal and disance from O
print '      normal      =    $n     = ', $n, "\n";
print '  disance from O  =    $d     = ', $d, "\n";
print ' Y axis intersect = $d/($n.Y) = ', $d/($n.Y), "\n";
print "\n";

use Math::MatrixReal;  # Not required for pure vector math as above

print "Conversions to MatrixReal objects\n";
$r = $a->vector2matrix_row;  # convert to MatrixReal Row Vector
$c = $a->vector2matrix_col;  # convert to MatrixReal Column Vector
print 'Vector as a MatrixReal Row $r (vector -> matrix) => ', "\n", $r;
print 'Vector as a MatrixReal Col $c (vector -> matrix) => ', "\n", $c;
print "\n";

print "Rotation Matrix from 3 Vectors\n";
$nx = $a->norm;   $ny = $nx x Z;  $nz = $nx x $ny; # orthogonal vectors
$R = vector_matrix( $nx, $ny, $nz );   # make the rotation matrix
print '$R   => ',"\n", $R, "\n";

print "Extract the Y row from the matrix as a VectorReal\n";
print '$R->matrix_row2vector(1) => ', $R->matrix_row2vector(1), "\n\n";

print "Rotate a vector with above rotation matrix\n";
print '$a * $R (vector -> vector)',"\n", $a * $R, "\n\n";

print "Rotate a MatrixReal column (post multiply)\n";
print "(NB: matrix must be transposed (~) to match column format)\n";
print '~$R * $c (col_matrix -> col_matrix) =>',"\n", ~$R * $c, "\n";

