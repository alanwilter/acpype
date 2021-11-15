#
# Math::VectorReal     Vector Mathematics 
#
#
# Copyright (c) 2001 Anthony Thyssen. All rights reserved. This program
# is free software; you can redistribute it and/or modify it under the
# same terms as Perl itself.
#
package Math::VectorReal;

=head1 NAME

Math::VectorReal - Module to handle 3D Vector Mathematics

=head1 SYNOPSIS

    #!/usr/bin/perl
    use Math::VectorReal;

    $a = vector( 1, 2, .5 );
    print "Vector as string (MatrixReal default format)\n\$a => ", $a;

    print  $a->stringify("Formated Output   \$a => { %g, %g, %g }\n");

    # I hate newline in the default output format (defined as MatrixReal)
    $Math::VectorReal::FORMAT = "[ %.5f %.5f %.5f ]";
    print "Modified default output format   \$a => $a\n";

    print 'length     => ', $a->length, "\n";
    print 'normalised => ', $a->norm, "\n";

    use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
    print 'string concat   $a."**" = ', $a."**", "\n";
    print 'vector constant    X    = ',   X,    "\n";
    print 'subtraction     $a - Z  = ', $a - Z, "\n";
    print 'scalar divide   $a / 3  = ', $a / 3, "\n";
    print 'dot product     $a . Y  = ', $a . Y, "\n";
    print 'cross product   $a x Y  = ', $a x Y, "\n";

    print "Plane containing points X, \$a, Z (in anti-clockwise order)\n";
    ($n,$d) = plane( X, $a, Z ); # return normal and disance from O
    print '      normal      =    $n     = ', $n, "\n";
    print '  disance from O  =    $d     = ', $d, "\n";
    print ' Y axis intersect = $d/($n.Y) = ', $d/($n.Y), "\n";

    print "VectorReal and MatrixReal interaction\n\n";
    use Math::MatrixReal;  # Not required for pure vector math as above

    $r = $a->vector2matrix_row;  # convert to MatrixReal Row Vector
    $c = $a->vector2matrix_col;  # convert to MatrixReal Column Vector
    print 'Vector as a MatrixReal Row $r (vector -> matrix) => ', "\n", $r;
    print 'Vector as a MatrixReal Col $c (vector -> matrix) => ', "\n", $c;

    $nx = $a->norm;   $ny = $nx x Z;  $nz = $nx x $ny; # orthogonal vectors
    $R = vector_matrix( $nx, $ny, $nz );   # make the rotation matrix
    print 'Rotation Matrix from 3 Vectors   $R   => ',"\n", $R, "\n";

    print "Extract the Y row from the matrix as a VectorReal..\n";
    print '$R->matrix_row2vector(1) => ', $R->matrix_row2vector(1), "\n";

    print "Rotate a vector with above rotation matrix\n";
    print '$a * $R (vector -> vector)',"\n", $a * $R, "\n";

    print "Rotate a MatrixReal column (post multiply)...\n";
    print "(NB: matrix must be transposed (~) to match column format)\n";
    print '~$R * $c (col_matrix -> col_matrix) =>',"\n", ~$R * $c, "\n";

=head1 DESCRIPTION

The C<Math::VectorReal> package defines a 3D mathematical "vector", in a way
that is compatible with the previous CPAN module C<Math::MatrixReal>. However
it provides a more vector oriented set of mathematical functions and overload
operators, to the C<MatrixReal> package.  For example the normal perl string
functions "x" and "." have been overloaded to allow vector cross and dot
product operations. Vector math formula thus looks like vector math formula in
perl programs using this package.

While this package is compatible with Math::MatrixReal, you DO NOT need to
have that package to perform purely vector orientated calculations. You will
need it however if you wish to do matrix operations with these vectors. The
interface has been designed with this package flexibility in mind.

The vectors are defined in the same way as a "row" C<Math::MatrixReal> matrix,
instead of that packages choice of "column" definition for vector operations.
Such vectors are multiplied to matices with the vector on the left and the
matrix on the right. EG:   v * M -> 'v

Not only is this the way I prefer to handle vectors, but it is the way most
graphics books use vectors. As a bonus it results in no overload conflicts
between this package and that of Math::MatrixReal, (the left objects overload
operator is called to do the mathematics). It also is a lot simpler than
C<MatrixReal> column vector methods, which were designed for equation solving
rather than 3D geometry operations.

The  vector_matrix()  function provided, simplifies the creation a
C<MatrixReal> object from 3 (usually orthogonal) vectors. This with its vector
orientated math operators makes it very easy to define orthogonal rotation
matrices from C<Math::VectorReal> objects.  See a rough example in the
synopsis above, or in the file "matrix_test" in the packages source.

NOTE: the 6th element the C<Math::MatrixReal> array object is used to hold the
length of the vector so that it can be re-used without needing to be
re-calculated all the time. This means the expensive sqrt() function, need not
be called unless nessary.  This usage should not effect the direct use of
these objects in the C<Math::MatrixReal> functions.

=cut

use strict;
#require Math::MatrixReal;  # not required!

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw( vector plane vector_matrix );

@EXPORT_OK = qw( O X Y Z );

%EXPORT_TAGS = (
   axis => [ qw( O X Y Z ) ],     # Unix Axis Vector Constants
   all  => [@EXPORT, @EXPORT_OK]
);

$VERSION = '1.0';

use Carp;
use vars qw( $FORMAT $TRACE );
$TRACE = 0;
$FORMAT  = "[ %#19.12E %#19.12E %#19.12E ]\n"; # output format (as MatrixReal)

=head1 CONSTANTS

Four constant vectors are available for export (using an ":all" tag).
these are

    0 = [ 0 0 0 ]   the zero vector or origin
    X = [ 1 0 0 ]   |
    Y = [ 0 1 0 ]    > Unit axis vectors
    Z = [ 0 0 1 ]   |

=cut 

# Constant Vector Functions
# The format is as per a Math::MatrixReal object, with extra length item
sub O() { bless [ [[0,0,0]], 1,3, undef,undef,undef, 0 ], __PACKAGE__; }
sub X() { bless [ [[1,0,0]], 1,3, undef,undef,undef, 1 ], __PACKAGE__; }
sub Y() { bless [ [[0,1,0]], 1,3, undef,undef,undef, 1 ], __PACKAGE__; }
sub Z() { bless [ [[0,0,1]], 1,3, undef,undef,undef, 1 ], __PACKAGE__; }

=head1 CONSTRUCTORS

=over 4

=item new(x,y,z)

Create a new vector with the values of C<x>, C<y>, C<z> returning the
appropriate object.

=item vector(x,y,z)

As C<new> but is a exported function which does not require a package
reference to create a C<Math::VectorReal> object.

=item clone()

Return a completely new copy of the referring C<Math::VectorReal> object.

=cut

sub new {  # typical object creation (not many checks)
  croak "Usage: \$vector = ".__PACKAGE__."->new(x,y,z);\n" unless @_;
  my $ref = shift;
  return bless [ [[ @_ ]], 1,3 ], ref $ref || $ref;
}


sub vector {  # normal way to create a vector - Exported function
              # This works as both a Object Method or Exported Function
  croak "Usage: \$vector = ".__PACKAGE__."->vector(x,y,z);\n".
        "  or   \$vector = vector(x,y,z);\n"
                          unless @_ == 3  ||  @_ == 4 && ref $_[0];
  my $class = __PACKAGE__;
  $class = ref shift  if @_ == 4;
  return $class->new(@_);
}

sub clone {
  croak "Usage: \$vector_copy = \$vector->clone;\n" unless @_ == 1;
  my $v = shift;
  my $c = $v->new( $v->array );          # create a new vector using values
  $c->[6] = $v->[6] if defined $v->[6];  # also note its length (if known)
  return $c;
}

=head1 METHODS

=item array()

Return the x,y,z elements of the referring vector are an array of values.

=item x()

Return the x element of the referring vector.

=item y()

Return the y element of the referring vector.

=item z()

Return the z element of the referring vector.

=item stringify( [ FORMAT ] )

Return the referring verctor as a string. The C<FORMAT> if given is used
to sprintf format the vector. This is used for all VectorReal to String
conversions.

By default this format is the same as it would be for a C<Math::MatrixReal>
object, "[ %#19.12E %#19.12E %#19.12E ]\n".  Note that this includes a newline
character!.

However unlike C<Math::MatrixReal> you can assign a new default sprintf
format by assigning it to the packages C<$FORMAT> variable. For Example

   $Math::VectorReal::FORMAT = "{ %g, %g, %g }"

Which is a good format to output vectors for use by the POVray (Persistance of
Vision Raytracer) program.

=item length()

Return the length of the given vector. As a side effect the length is saved
into that vectors object to avoid the use of the expensive sqrt() function.

=item norm()

Normalise the Vector. That is scalar divide the vector by its length, so that
it becomes of length one.  Normal vectors are commonly use to define
directions, without scale, or orientation of a 3 dimensional plane.

=cut

sub array {   # return vector as an array of values
  my $v = shift;
  return @{$v->[0][0]};
}

sub x {
  my $v = shift;
  return ($v->array)[0];
}

sub y {
  my $v = shift;
  return ($v->array)[1];
}

sub z {
  my $v = shift;
  return ($v->array)[2];
}

sub stringify {   # convert a vector to a string (with optional format)
  my( $v, $fmt ) = @_;
  $fmt = $FORMAT   unless defined $fmt;  # if not given use current default
  return sprintf $fmt, $v->array;
}

sub length {   # convert a vector to a string
  my $v = shift;
  return $v->[6] if defined $v->[6];
  return $v->[6] = sqrt( $v.$v );
}

sub norm {   # scale vector to a length of one
  my $v = shift;
  return $v / $v->length;
}

=item plane( v1, v2, v3 )

Given three points defined counter clockwise on a plane, return an array in
which the first element is the planes normal unit vector, and the second its
distance from the origin, along that vector.  NOTE: the distance may be
negitive, in which case the origon is above the defined plane in 3d space.

=cut

sub plane { # Given three points on the plane (right-hand rule)
            # return a normal vector and distance from origin for a plane
  croak "Usage: (\$normal, \$distance) = plane(\$p1,\$p2,\$p3);\n"
                                                         unless @_ == 3;
  my ($a, $b, $c) = @_;
  my $normal = (($b - $a) x ($c - $b))->norm;
  return ( $normal, $a . $normal );
}

=item vector_matrix( nx, ny, nz )

Given the new location for the X, Y and Z vectors, concatanate them together
(row wise) to create a C<Math::MatrixReal> translation matrix. For example
if the 3 vectors are othogonal to each other, the matrix created will be
a rotation matrix to rotate the X, Y and Z axis to the given vectors. See
above for an example.

=cut

sub vector_matrix {
  my( $nx, $ny, $nz ) = @_;
  bless [ [[$nx->array],
           [$ny->array],
           [$nz->array]], 3, 3 ],  "Math::MatrixReal";
}

# ------------------------------------------------------------------
# Convertsions between Math::MatrixReal and Math::VectorReal packages

=back

=head1 VECTOR/MATRIX CONVERSION

The following functions provide links between the C<Math::VectorReal> and
C<Math::MatrixReal> packages.

NOTE: While this package is closely related to C<Math::MatrixReal>, it does
NOT require that that package to be installed unless you actually want to
perform matrix operations.

Also the overload operations will automatically handle vector/matrix
mathematics (See below).

=head2 Vector to Matrix Conversion

=over 4

=item vector2matrix_row( [CLASS] )

=item vector2matrix_col( [CLASS] )

Convert C<Math::VectorReal> objects to a C<Math::MatrixReal> objects.
Optional argument defines the object class to be returned (defaults to
C<Math::MatrixReal>).

Note that as a C<Math::VectorReal> is internally equivelent to a
C<Math::MatrixReal> row matrix, C<vector2matrix_row> is essentually just a
bless operation, which is NOT required to use with C<Math::MatrixReal>
functions.

The C<vector2matrix_col> performs the required transpose to convert the
C<Math::VectorReal> object into a C<Math::MatrixReal> version of a vector (a
column matrix).

=cut

sub vector2matrix_row {
  my( $v, $ref ) = @_;
  $ref ||= "Math::MatrixReal";
  bless $v->clone,  ref $ref || $ref;  # clone and bless (object unchanged)
}

sub vector2matrix_col {
  my( $v, $ref ) = @_;
  $ref ||= "Math::MatrixReal";
  my @v = $v->array;
  bless [ [[$v[0]],[$v[1]],[$v[2]]], 3, 1 ], ref $ref || $ref;
}

=head2 Matrix to Vector Conversion

=item matrix_row2vector( [ROW] )

=item matrix_col2vector( [COLUMN] )

When referred to by a C<Math::MatrixReal> object, extracts the vector
from the matrix. the optional argument defines which row or column of the
matrix is to be extracted as a C<Math::VectorReal> vector.

=cut

{  # Enclose MartixReal package in a block
package Math::MatrixReal; # Fake a change into the Math::MatrixReal package
use Carp;                 # import carp into this package

sub matrix_row2vector {
  my $m = shift;    my($rows,$cols) = ($m->[1],$m->[2]);
  my $r = shift;   # optional, which column from matrix
  croak "Error: matrix does not have 3D rows" unless ($cols == 3);
  if ( defined $r ) {
    croak "Error: matrix does not have that row" unless ( $r < $rows);
  }
  else {    # if no option, it must be a Math::MatrixReal Row Vector
    croak "Error: matrix given to matrix_row2vector is not a 3D row matrix"
           unless ($rows == 1);
    $r = 0;
  }
  return Math::VectorReal->new(@{$m->[0][$r]}); # same result, only cleaned up
}

sub matrix_col2vector {
  my $m = shift;    my($rows,$cols) = ($m->[1],$m->[2]);
  my $c = shift;   # optional, which column from matrix
  croak "Error: matrix does not have 3D rows" unless ($rows == 3);
  if ( defined $c ) {
    croak "Error: matrix does not have that column" unless ( $c < $cols);
  }
  else {    # if no option, it must be a Math::MatrixReal Column Vector
    croak "Error: matrix given to matrix_col2vector is not a 3D column matrix"
           unless ($cols == 1);
    $c = 0;
  }
  return Math::VectorReal->new($m->[0][0][$c], $m->[0][1][$c], $m->[0][2][$c]);
}

} # Return to the Math::VectorReal package we are really defining

# ------------------------------------------------------------------
# Overloaded Math functions

=back

=head1 OPERATOR OVERLOADING

Overload operations are provided to perform the usual string conversion,
addition, subtraction, unary minus, scalar multiplation & division.  On top of
this however the multiply have been expanded to look for and execute
C<MatrixReal> multiplation.

The Main purpose of this package however was to provide the special vector
product operations: dot product "." and cross product "x".  In perl these
operations are normally used for string operations, but if either argument
is a C<VectorReal> object, the operation will attempt the approprate
vector math operation instead.

Note however that if one side of the dot "." operator is already a string,
then the vector will be converted to a sting and a string concatantion will be
performed. The cross operator "x" will just croak() as it is non-sensical to
either repeat the string conversion of a vector, OR to repeat a string,
vector, times!

Overloaded operator summery...
    neg     unary minus - multiply vector by -1
     ""     automatic string conversion using stringify() function
      +     vector addition
      -     vector subtraction
      /     scalar division (left argument must be the vector)
      *     scalar multiplication OR MatrixReal multiplication
      x     vector/cross product of two vectors
      .     dot product of two vectors OR vector/string concatanation

Posible future addition   '~'  to transpose a C<VectorReal> into a
C<MatrixReal> column vector (as per that operator on C<MatrixReal> objects).
It was not added as it just did not seem to be needed.

=cut

use overload
     'neg' => \&_negate,
      '""' => \&_stringify,
       '+' => \&_addition,
       '-' => \&_subtract,
       '*' => \&_multiply,
       '/' => \&_scalar_divide,
       'x' => \&_cross_product, # Redefination of the string function
       '.' => \&_dot_product,   # These includes stingify/concatanation
'fallback' => undef;


sub _trace {
    return unless $TRACE;
    my($text,$object,$argument,$flip) = @_;
    unless (defined $object)   { $object   = 'undef'; };
    unless (defined $argument) { $argument = 'undef'; };
    unless (defined $flip)     { $flip     = 'undef'; };
    if (ref($object))   { $object   = ref($object);   }
    if (ref($argument)) { $argument = ref($argument); }
    $argument =~ s/\n/\\n/g;
    print "$text: \$obj='$object' \$arg='$argument' \$flip='$flip'\n";
}


sub _negate {
  my($object,$argument,$flip) = @_;
  _trace("'neg'",$object,$argument,$flip);
  my $v = $object->clone;
  for ( 0 .. 2 ) { $v->[0][0][$_] = -$v->[0][0][$_]; }
  # $v->[6]; does not change.
  return $v
}

sub _stringify {
  my($object,$argument,$flip) = @_;
  _trace("'\"\"'",$object,$argument,$flip);
  return $object->stringify;
}


sub _addition {
  # Operation on two vectors, as such $flip will be undefined or false
  # The operation is also communitive - order does not matter.
  my($object,$argument,$flip) = @_;
  _trace("'+'",$object,$argument,$flip);
  if ( (defined $argument) && ref($argument) &&
       (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/) ) {
    my $v = $object->clone;
    for ( 0 .. 2 ) { $v->[0][0][$_] += $argument->[0][0][$_]; }
    $#{$v} = 2;   # any cached vector length is now invalid
    return $v;
  }
  croak("non-vector argument for '+'");
}


sub _subtract {
  my($object,$argument,$flip) = @_;
  _trace("'-'",$object,$argument,$flip);
  # Operation on two vectors, as such $flip will be undefined or false
  # Note; however this is not communitive - order matters
  if ( (defined $argument) && ref($argument) &&
       (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/) ) {
    my $v = $object->clone;
    for ( 0 .. 2 ) { $v->[0][0][$_] -= $argument->[0][0][$_]; }
    $#{$v} = 2;    # any cached vector length is now invalid
    return $v;
  }
  croak("non-vector argument for '-'");
}


sub _multiply {
  my($object,$argument,$flip) = @_;
  _trace("'*'",$object,$argument,$flip);
  if ( ref($argument) ) {
    # Assume multiply by  Math::MatrixReal object  EG:  $v * $M --> $new_v
    # Order is communicative, but $flip should NOT be true
    if ( ! $flip ) {
      return ( $object->vector2matrix_row($argument)
                        * $argument )->matrix_row2vector;
    } else { # just in case flip is true..
      return ( $argument *
                $object->vector2matrix_row($argument) )->matrix_row2vector;
    }
  }
  elsif ( defined $argument ) {
    # defined $argument must be a scalar, so Scalar Multiply
    # Communitive - order does not matter, $flip can be ignored
    my $v = $object->clone;
    for ( 0 .. 2 ) { $v->[0][0][$_] *= $argument; }
    $v->[6] *= abs($argument) if defined $v->[6]; # multiply vector length
    return $v;
  }
  croak("undefined argument given for vector multiply");
}


sub _scalar_divide {
  my($object,$argument,$flip) = @_;
  _trace("'/'",$object,$argument,$flip);
  # The order is very important, you can NOT divide a scalar by a vector
  croak("You can not divide a scalar by a vector") if $flip;
  # The provided $argument must be a defined scalar
  if ( (defined $argument) && ! ref($argument)  ) {
    my $v = $object->clone;
    for ( 0 .. 2 ) { $v->[0][0][$_] /= $argument; }
    $v->[6] /= abs($argument) if defined $v->[6]; # do vector length
    return $v;
  }
  croak("non-scalar given for vector scalar divide");
}


sub _cross_product {
  my($object,$argument,$flip) = @_;
  # Operation on two vectors, as such $flip will be undefined or false
  # Note: however this is not communitive - order does matters
  _trace("'x'",$object,$argument,$flip);
  if ( (defined $argument) && ref($argument) &&
       (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/) ) {
    my $v = $object->new;
    my @o = $object->array;
    my @a = $argument->array;
    @{$v->[0][0]} = ( $o[1]*$a[2] - $o[2]*$a[1],
                      $o[2]*$a[0] - $o[0]*$a[2],
                      $o[0]*$a[1] - $o[1]*$a[0] );
    $#{$v} = 2;    # any cached vector length is now invalid
    return $v;
  }
  croak("string 'x' with a vector does not make sense!");
}


sub _dot_product {
  my($object,$argument,$flip) = @_;
  if ( (defined $argument) && ref($argument) &&
       (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/) ) {
    # Operation on two vectors, and communitive - order does not matter
    _trace("'.'",$object,$argument,$flip);
    my $v = 0;   # result is NOT an object, but a scalar
    for ( 0 .. 2 ) { $v +=  $object->[0][0][$_] * $argument->[0][0][$_]; }
    return $v;
  }
  # Argument is NOT a vector! Assume String concatenation wanted
  elsif ( defined $flip ) {
    if ( $flip ) {
      _trace("'.\"\"'",$object,$argument,$flip);
      return $argument . $object->stringify;
    } else {
      _trace("'\"\".'",$object,$argument,$flip);
      return $object->stringify . $argument;
    }
  }
  # concatenate a string to a vector
  _trace("'.='",$object,$argument,$flip);
  return $object->stringify . $argument;
  # Concatenate a vector to string is handled automatically with '""' operator
}

1;
# ------------------------------------------------------------------

=head1 SEE ALSO

The C<Math::MatrixReal> CPAN Module by   Steffen Beyer
and the C<Math::MatrixReal-Ext1> CPAN extension by  Mike South

=head1 AUTHOR

Anthony Thyssen E<lt>F<anthony@cit.gu.edu.au>E<gt>

=head1 COPYRIGHT

Copyright (c) 2001 Anthony Thyssen. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself. I would appreciate any suggestions however.

=cut

