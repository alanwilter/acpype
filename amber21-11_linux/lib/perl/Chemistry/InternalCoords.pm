package Chemistry::InternalCoords;

$VERSION = '0.18';
# $Id: InternalCoords.pm,v 1.6 2004/09/24 20:48:23 itubert Exp $

use 5.006;
use strict;
use warnings;
use Math::VectorReal qw(:all);
use Carp;
use overload '""' => "stringify", fallback => 1;
use Scalar::Util 'weaken';
use Math::Trig 'deg2rad';


=head1 NAME

Chemistry::InternalCoords - Represent the position of an atom using internal
coordinates and convert it to Cartesian coordinates. 

=head1 SYNOPSIS
    
    use Chemistry::InternalCoords;

    # ... have a molecule in $mol
    my $atom = $mol->new_atom;

    # create an internal coordinate object for $atom 
    # with respect to atoms with indices 4, 3, and 2.
    my $ic = Chemistry::InternalCoords->new(
        $atom, 4, 1.1, 3, 109.5, 2, 180.0
    );

    # can also use atom object references instead of indices
    ($atom4, $atom3, $atom2) = $mol->atoms(4,3,2);
    my $ic = Chemistry::InternalCoords->new(
        $atom, $atom4, 1.1, $atom3, 109.5, $atom2, 180.0
    );

    # calculate the Cartesian coordinates for
    # the atom from the internal coordinates
    my $vector = $ic->cartesians;

    # calculate and set permanently the Cartesian coordinates
    # for the atom from the internal coordinates
    my $vector = $ic->add_cartesians;
    # same as $atom->coords($ic->cartesians);

    # dump as string
    print $ic;
    # same as print $ic->stringify;

=head1 DESCRIPTION

This module implements an object class for representing internal coordinates 
and provides methods for converting them to Cartesian coordinates.

For generating an internal coordinate representation (aka a Z-matrix) of a
molecule from its Cartesian coordinates, see the
L<Chemistry::InternalCoords::Builder> module.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 METHODS

=over

=item my $ic = Chemistry::InternalCoords->new($atom, $len_ref, $len_val,
$ang_ref, $ang_val, $dih_ref, $dih_val)

Create a new internal coordinate object. $atom is the atom to which the
coordinates apply. $len_ref, $ang_ref, and $dih_ref are either atom references
or atom indices and are used to specify the distance, angle, and dihedral that
are used to define the current position. $len_val, $ang_val, and $dih_val are
the values of the distance, angle, and dihedral. The angle and the dihedral are
expected to be in degrees.

For example,

    my $ic = Chemistry::InternalCoords->new(
        $atom, 4, 1.1, 3, 109.5, 2, 180.0
    );

means that $atom is 1.1 distance units from atom 4, the angle $atom-4-3 is
109.5 degrees, and the dihedral $atom-4-3-2 is 180.0 degrees.

The first three atoms in the molecule don't need all the internal coordinates:
the first atom doesn't need anything (except for the atom reference $atom) 
because it will always be placed at the origin; the second atom only needs
a distance, and it will be placed on the X axis; the third atom needs a 
distance and an angle, and it will be placed on the XY plane.

=cut

sub new {
    my $class = shift;
    my $atom = shift;

    my $self = bless {
        atom => $atom,
        len_ref => shift, len_val => shift,
        ang_ref => shift, ang_val => shift,
        dih_ref => shift, dih_val => shift,
    }, ref $class || $class;

    weaken($self->{atom}); # to avoid memory leaks

    for (@$self{qw(len_ref ang_ref dih_ref)}) {
        if ($_ and not ref) {
            $_ = $atom->parent->atoms($_); 
        }
    }

    $self;
}

=item my ($atom, $distance) = $ic->distance

Return the atom reference and distance value contained in the
Chemistry::InternalCoords object.

=cut

sub distance {
    my ($self) = @_;
    ($self->{len_ref}, $self->{len_val});
}

=item my ($atom, $angle) = $ic->angle

Return the atom reference and angle value contained in the
Chemistry::InternalCoords object.

=cut

sub angle {
    my ($self) = @_;
    ($self->{ang_ref}, $self->{ang_val});
}

=item my ($atom, $dihedral) = $ic->dihedral

Return the atom reference and dihedral value contained in the
Chemistry::InternalCoords object.

=cut

sub dihedral {
    my ($self) = @_;
    ($self->{dih_ref}, $self->{dih_val});
}

=item my $vector = $ic->cartesians

Calculate the Cartesian coordinates from an internal coordinate object.
Returns a Math::VectorReal object. Note that the Cartesian coordinates of the
atoms referenced by the $ic object should already be calculated.

=cut

sub cartesians {
    my ($self) = @_;

    #print "cartesians\n";
    unless ($self->{len_ref}) { # origin
        return vector(0, 0, 0); 
    }

    unless ($self->{ang_ref}) { # second atom; place on X axis
        return vector($self->{len_val}, 0, 0);
    }

    if (!$self->{dih_ref}  # third atom; place on XY plane 
        or abs($self->{ang_val} % 180) < 0.005) { # linear angle
        my $len = $self->{len_val};
        my $ang = deg2rad(180 - $self->{ang_val});
        my $d1 = $self->{len_ref}->coords - $self->{ang_ref}->coords;
        $d1 = $d1->norm;
        my $v = $len * $d1 * cos($ang) + $len * Y * sin($ang);
        #vector($len * cos($ang), $len * sin($ang), 0) ;
        #print "len=$len; ang=$ang; v=$v; d1=$d1\n";
        return( $v + $self->{len_ref}->coords);
        
    }

    # the real thing...
    my $v1 = $self->{dih_ref}->coords; # 'oldest' point
    my $v2 = $self->{ang_ref}->coords;
    my $v3 = $self->{len_ref}->coords; # 'newest' point
    my $d1 = $v1 - $v2;
    my $d2 = $v3 - $v2;

    # $xp = normal to atoms 1 2 3 
    my $xp = $d1 x $d2;

    if ($xp->length == 0) { # ill-defined dihedral
        my $len = $self->{len_val};
        my $ang = deg2rad(180 - $self->{ang_val});
        my $d1 = $self->{len_ref}->coords - $self->{ang_ref}->coords;
        $d1 = $d1->norm;
        my $v = $len * $d1 * cos($ang) + $len * Y * sin($ang);
        return( $v + $self->{len_ref}->coords);
    }

    # $yp = normal to xp and atoms 2 3
    my $yp = $d2 x $xp;

    my $ang1 = deg2rad($self->{dih_val});   # dihedral
    # $r = normal to atoms 2 3 4 (where 4 is the new atom)
    #      obtained by rotating $xp through $d2
    my $r = $xp->norm * cos($ang1) + $yp->norm * sin($ang1);

    my $ypp = $d2 x $r; # complete new frame of reference
    my $ang2 = deg2rad($self->{ang_val});   # angle
    my $d3 = -$d2->norm * cos($ang2) + $ypp->norm * sin($ang2);

    $d3 = $d3 * $self->{len_val}; # mult by distance to $v3
    my $v4 = $v3 + $d3; # define new point

    return $v4;
}

=item my $vector = $ic->add_cartesians

Same as $ic->cartesians, but also adds the newly calculated Cartesian
coordinates to the atom. It is just shorthand for the following:

    $atom->coords($ic->cartesians);

The best way of calculating the Cartesian coordinates for an entire molecule,
assuming that every atom is defined only in terms of previous atoms (as it 
should be), is the following:

    # we have all the internal coords in @ics
    for my $ic (@ics) {
        $ic->add_cartesians;
    }

=cut

sub add_cartesians {
    my ($self) = @_;
    my $v = $self->cartesians;
    $self->{atom}->coords($v);
    $v;
}

=item my $string = $ic->stringify

Dump the object as a string representation. May be useful for debugging.
This method overloads the "" operator.

=cut

sub stringify {
    my ($self) = shift;
    no warnings 'uninitialized';
    my $ret;
    for my $key (qw(len_ref len_val ang_ref ang_val dih_ref dih_val)) {
        $ret .= "$key=($self->{$key}), ";
    }
    "$ret\n";
}

=back 

=head1 VERSION

0.18

=head1 SEE ALSO

L<Chemistry::InternalCoords::Builder>,
L<Chemistry::Mol>, L<Chemistry::Atom>, 
L<Math::VectorReal>, L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

1;

