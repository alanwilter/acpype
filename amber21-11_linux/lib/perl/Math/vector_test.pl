#!/usr/bin/perl
#
# Testing the Math::VectorReal module
#
# location of my libraries (I normally set this via PERL5LIB)
#use lib '/home/anthony/lib/perl5';  # location of my libraries
use FindBin;
use lib "$FindBin::Bin/blib/lib";

use Math::VectorReal qw( :all );
#$Math::VectorReal::TRACE = 1; # Tracing of overload operators on/off

print "Stringify output testing (MatrixReal default)\n";
print 'O->stringify => ', O->stringify;
print "\n";

print "Changing default vector to string format\n";
print '$Math::VectorReal::FORMAT = "[ %g %g %g ]";', "\n";
$Math::VectorReal::FORMAT = "[ %g %g %g ]";
print "\n";

print "Axis functions, assign to constants\n";
print ' $o = O => ', $o=O, "\n";
print ' $x = X => ', $x=X, "\n";
print ' $y = Y => ', $y=Y, "\n";
print ' $z = Z => ', $z=Z, "\n";
print "\n";

print "String conversion operation testing\n";
print "Note: this include some automatic stringify concat ('.') operations\n";
print ' "$o"  => ', "$o", "\n";
print '""$x   =>', " $x", "\n";
print '  $y"" => ', "$y\n";
print '  $z   => ',  $z, "\n";
print 'vector(1,2,3) => ', vector(1,2,3), "\n";
print "\n";

print "Addition\n";
$a = $x + Y;
print '$a = $x + Y => ', $a, "\n";
$a += $y ;
print '$a += $y    => ', $a, "\n";
print "\n";

print "Clone and Addition Tests\n";
$b = $y;
print '$b = $y  => ', $b, "\n";
$b += Z;
print '$b += Z  => ', $b, "\n";
print '   $y    => ', $y, "\n";
print "\n";

print "Subtraction\n";
$b -= $z ;
print '$b -= $z    => ', $b, "\n";
$b = $b - Z ;
print '$b = $b - Z => ', $b, "\n";
print "\n";


print "Scalar Multiply\n";
$a = $z * 2;
print '$a = $z * 2 => ', $a, "\n";
$a = 2 * Z;
print '$a = 2 * Z  => ', $a, "\n";
$a *= 2.5;
print '$a *= 2.5   => ', $a, "\n";
print "\n";

print "Scalar Divide\n";
$a = $b / 2;
print '$a = $b / 2 => ', $a, "\n";
$a /= 3e14;
print '$a /= 3e14  => ', $a, "\n";
print "\n";


print "Unary - and more subtraction\n";
$b = -$b;
print '$b = -$b       => ', $b, "\n";
$b -= Z;
print '$b -= Z        => ', $b, "\n";
$b -= $z - -$y;
print '$b -= $z - -$y => ', $b, "\n";
$b = $o - $b;
print '$b = $o - $b   => ', $b, "\n";
print "\n";


print "Cross Product\n";
$a = $b x X;
print '$a = $b x X   => ', $a, "\n";
$a = $b x $y;
print '$a = $b x $y  => ', $a, "\n";
$a = $b x $z;
print '$a = $b x $z  => ', $a, "\n";
print "\n";


print "Dot Product / String Concatenation\n";
$a = Z . $b;
print '$a = Z . $b   => ', $a, "\n";
$a = $b . -$y;
print '$a = $b . -$y => ', $a, "\n";
$s = $b . "!";
print '$s = $b . "!" => ', $s, "\n";
$s = "!" . $b;
print '$s = "!" . $b => ', $s, "\n";
$a .= $b;
print '$a .= $b      => ', $a, "\n";
print "\n";

print "Special Functions (length, norm, plane)\n";
print '$b->length    => ', $b->length, "\n";
print '$b->norm      => ', $b->norm, "\n";
@a = plane(X,Y,Z);
print '@a = plane(X,Y,Z) => ', "\n  @a\n";
print "check output from plane() function\n";
$a = X+Y+Z;
print 'normal   => ', $a->norm, "\n";
print 'distance => ', ($a/3)->length, "\n";
print "\n";

print "Are defined constants still OK\n";
print '$o => ', $o, "\n";
print '$x => ', $x, "\n";
print '$y => ', $y, "\n";
print '$z => ', $z, "\n";
print "\n";

