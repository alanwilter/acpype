package Chemistry::Ring;
$VERSION = '0.19';
#$Id: Ring.pm,v 1.1.1.1 2005/03/29 23:57:36 itubert Exp $

=head1 NAME

Chemistry::Ring - Represent a ring as a substructure of a molecule

=head1 SYNOPSIS

    use Chemistry::Ring;
    
    # already have a molecule in $mol...
    # create a ring with the first six atoms in $mol
    my $ring = Chemistry::Ring->new;
    $ring->add_atom($_) for $mol->atoms(1 .. 6);

    # find the centroid
    my $vector = $ring->centroid;

    # find the plane that fits the ring
    my ($normal, $distance) = $ring->plane;

    # is the ring aromatic?
    print "is aromatic!\n" if $ring->is_aromatic;

    # "aromatize" a molecule
    Chemistry::Ring::aromatize_mol($mol);

    # get the rings involving an atom (after aromatizing)
    my $rings = $mol->atoms(3)->attr('ring/rings');

=head1 DESCRIPTION

This module provides some basic methods for representing a ring. A ring is
a subclass of molecule, because it has atoms and bonds. Besides that, it
has some useful geometric methods for finding the centroid and the ring plane,
and methods for aromaticity detection.

This module does not detect the rings by itself; for that, look at 
L<Chemistry::Ring::Find>.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=cut

use strict;
use warnings;
use Math::VectorReal qw(:axis vector);
use Statistics::Regression;
use Chemistry::Mol;
use base 'Chemistry::Mol', 'Exporter';
use Scalar::Util 'weaken';

our @EXPORT_OK = qw(aromatize_mol);
our %EXPORT_TAGS = ( all => \@EXPORT_OK ); 

our $N = 0;
our $DEBUG = 0;

=head1 METHODS

=over 4

=item Chemistry::Ring->new(name => value, ...)

Create a new Ring object with the specified attributes. Same as
C<< Chemistry::Mol->new >>.

=cut

sub nextID {
    "ring".++$N; 
}


# make sure we don't become parent of the atoms added to us
sub add_atom { shift->SUPER::add_atom_np(@_) }
sub add_bond { shift->SUPER::add_bond_np(@_) }

sub print {
    my $self = shift;
    return <<EOF;
    ring:
        id: $self->{id}
        atoms: @{$self->{atoms}}
        bonds: @{$self->{bonds}}
EOF
}

=item $ring->centroid

Returs a vector with the centroid, defined as the average of the coordinates
of all the atoms in the ring. The vecotr is a L<Math::VectorReal> object.

=cut

sub centroid {
    my $self = shift;
    my $c = O; # origin
    my $n = 0;
    for my $a ($self->atoms) {
	$c += $a->coords;
	++$n;
    }
    $c = $c / $n; 
}

=item my ($norm, $d) = $ring->plane

Returns the normal and distance to the origin that define the plane that best
fits the atoms in the ring, by using multivariate regression. The normal 
vector is a L<Math::VectorReal> object.

=cut

sub plane {
    my $self = shift;
    my $reg = Statistics::Regression->new(3, "plane for $self", [qw(b mx my)]);
    for my $atom ($self->atoms) {
        my ($x, $y, $z) = $atom->coords->array;
        $reg->include($z, [1.0, $x, $y]);
    }
    $reg->print if $DEBUG;

    # convert the theta vector (z = a + bx + cy) to a normal vector and
    # distance to the origin
    my @coef = (@{$reg->theta}, -1.0);   # -1 is d in a + bx + cx + dz = 0
    my $d = shift @coef;                 # distance (not normalized)
    my $sum_sq = 0;                      # normalization constant
    $sum_sq += $_*$_ for @coef;   
    $sum_sq ||= 1;
    ($d, @coef) = map { $_ / $sum_sq } ($d, @coef); # normalize
    return (vector(@coef)->norm, $d);
}

=item $ring->is_aromatic

Naively guess whether ring is aromatic from the molecular graph, with a method
based on Hückel's rule. This method is not very accurate, but works for simple
molecules. Returns true or false.

=cut

sub is_aromatic {
    my ($self) = @_;
    my $n_pi = 0;

    for my $atom ($self->atoms) {
        no warnings 'uninitialized';
        return 0 if ($atom->bonds + $atom->hydrogens > 3);        

        # build bond order histogram
        my @order_freq = (0,0,0,0);
        for my $bond ($atom->bonds) {
            $order_freq[$bond->order]++;
        }
        
        return 0 if ($order_freq[3] or $order_freq[2] > 1);
        if ($order_freq[2] == 1) {
            $n_pi += 1;
        } elsif ($atom->symbol =~ /^[NOS]$/) {
            $n_pi += 2;
        }
    }
    #print "n_pi = $n_pi\n";
    return ($n_pi % 4 == 2) ? 1 : 0;
}

1;

=back

=head1 EXPORTABLE SUBROUTINES

Nothing is exported by default, but you can export these subroutines
explicitly, or all of them by using the ':all' tag.

=over

=item aromatize_mol($mol)

Finds all the aromatic rings in the molecule and marks all the atoms and bonds
in those rings as aromatic. 

It also adds the 'ring/rings' attribute to the molecule and to all ring atoms
and bonds; this attribute is an array reference containing the list of rings
that involve that atom or bond (or all the rings in the case of the molecule).
NOTE (the ring/rings attribute is experimental and might change in future
versions).

=cut

sub aromatize_mol {
    my ($mol) = @_;

    require Chemistry::Ring::Find;

    $_->aromatic(0) for ($mol->atoms, $mol->bonds);

    my @rings = Chemistry::Ring::Find::find_rings($mol);
    $mol->attr("ring/rings", \@rings);
    $_->attr("ring/rings", []) for ($mol->atoms, $mol->bonds);

    for my $ring (@rings) {
        if ($ring->is_aromatic) {
            $_->aromatic(1) for ($ring->atoms, $ring->bonds);
        }
        for ($ring->atoms, $ring->bonds) {
            my $ringlist = $_->attr("ring/rings") || [];
            push @$ringlist, $ring;
            weaken($ringlist->[-1]);
            $_->attr("ring/rings", $ringlist);
        }
    }
    @rings;
}

=back

=head1 VERSION

0.19

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Ring::Find>, 
L<Math::VectorReal>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

