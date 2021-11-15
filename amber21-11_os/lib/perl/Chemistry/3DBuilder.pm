package Chemistry::3DBuilder;
# $Id: 3DBuilder.pm,v 1.2 2005/05/06 21:34:46 itubert Exp $

$VERSION = '0.10';

use strict;
use warnings;
use Chemistry::Ring;
use Chemistry::InternalCoords;
use List::Util qw(first);

use base qw(Exporter);

our @EXPORT_OK = qw(build_3d);
our %EXPORT_TAGS = ( all => \@EXPORT_OK );

our $DEBUG = 0;

=head1 NAME

Chemistry::3DBuilder - Generate 3D coordinates from a connection table

=head1 SYNOPSIS

    # example: convert SMILES to MDL molfile
    use Chemistry::3DBuilder qw(build_3d);
    use Chemistry::File::SMILES;
    use Chemistry::File::MDLMol;

    my $s = '[O-]C(=O)C(N)C(C)CC';
    my $mol = Chemistry::Mol->parse($s, format => 'smiles');

    build_3d($mol);

    print $mol->print(format => 'mdl');

=head1 DESCRIPTION

This module generates a three-dimensional molecular structure from a connection
table, such as that obtained by a 2D representation of the molecule or from a
SMILES string.

B<NOTE>: this module is still at a very early stage of development so it has
important limitations. 1) It doesn't handle rings or stereochemistry yet!  2)
The bond lengths and atoms are very approximate as they don't really account
for different elements. 3) Only the sp3, sp2, and sp hybridizations are
supported.

=head1 SUBROUTINES

These subroutines may be exported; to export all, use the ':all' tag.

=over

=item build_3d($mol)

Add internal and cartesian coordinates to the molecule C<$mol>.

=cut

# POSSIBLE PROBLEMS:
# allenes
# stereochemistry
# (E)-esters and acids

sub build_3d {
    my ($mol) = @_;

    # prepare molecule
    $mol->collapse_hydrogens;
    $_->{internal_coords} = undef for $mol->atoms;
    Chemistry::Ring::aromatize_mol($mol);

    # check for rings
    my $rings = $mol->attr('ring/rings');
    if (@$rings) {
        warn "mol has " . @$rings . " rings\n";
        return;
    }

    my $user = { n => 0 };
    build_chain($mol, $user);

    #print "sprout_hydrogens\n" if $DEBUG;
    # add the hydrogens
    $mol->sprout_hydrogens;
    for my $h (
        grep { $_->symbol eq 'H' and ! $_->internal_coords } $mol->atoms
    ) {
        my ($bond) = $h->bonds;
        my ($nei)  = $h->neighbors;
        add_atom(atom => $h, bond => $bond, from => $nei, user => $user);
    }
    #warn "IC($_): " . $_->internal_coords for $mol->atoms;

    # calculate the cartesian coordinates
    $_->internal_coords->add_cartesians for $mol->atoms;
}

sub build_chain {
    my ($mol, $user) = @_;
    dfs($mol, on_atom => \&add_atom, user => $user);
}

sub bond_length {
    my ($bond) = @_;
    my (@atoms) = $bond->atoms;
    if (grep { $_->symbol eq 'H' } @atoms) {
        return 1.1;
    } else {
        return 1.5 - ($bond->order - 1) * 0.15;
    }
}

sub choose_bond {
    my ($atom, $bond, $from) = @_;
    ($from, bond_length($bond));
}

my %angle_table = (
    sp => 180,
    sp2 => 120,
    sp3 => 109.47,
);

sub bond_angle {
    my ($atom1, $atom2, $atom3) = @_;
    my $hyb = atom_hybridization($atom2);
    #warn "hyb angle($atom)=$hyb\n";
    $angle_table{$hyb} || 109.47;
}

sub atom_hybridization {
    no warnings qw(uninitialized);
    my ($atom) = @_;
    my %ord_hist;
    $ord_hist{$_}++ for map { $_->order } $atom->bonds;
    my $deg = $atom->bonds + $atom->implicit_hydrogens;
    if ($deg == 2 and ($ord_hist{3} or $ord_hist{2} == 2)) {
        return "sp";
    } elsif ($deg == 3 and $ord_hist{2} == 1) {
        return "sp2";
    } else {
        return "sp3";
    }
}

sub bond_dihedral {
    my ($atom, $from, $ang_ref) = @_;
    180;
}

sub choose_angle {
    my ($atom, $bond, $from) = @_;
    my ($len_ref, $len_val) = choose_bond($atom, $bond, $from);
    my ($ang_ref) = first { $_->internal_coords } $from->neighbors($atom);
    ($len_ref, $len_val, $ang_ref, bond_angle($atom, $from, $ang_ref));
}

sub choose_dihedral {
    my ($atom, $bond, $from) = @_;
    my ($len_ref, $len_val, $ang_ref, $ang_val) =
        choose_angle($atom, $bond, $from);
    my ($dih_ref, $dih_val);
    my $hyb = atom_hybridization($len_ref);
    my $fd;
   
    #warn "hyb($len_ref)=$hyb\n";
    if ($fd = $len_ref->attr("3dbuilder/first_dihedral")) {
        $fd = $atom->parent->by_id($fd); # XXX
        $dih_ref = $fd;
        ($ang_ref, $ang_val) = $fd->internal_coords->angle;
        $dih_val = $hyb eq 'sp2' ? 180 : 
            $len_ref->attr("3dbuilder/second_dihedral") ?
            240 : 120;
        $len_ref->attr("3dbuilder/second_dihedral", $atom->id);
    } elsif ($fd = $len_ref->attr("3dbuilder/first_dihedral_ang")) {
        $fd = $atom->parent->by_id($fd); # XXX
        ($dih_ref) = $fd->internal_coords->dihedral;
        ($ang_ref) = $fd->internal_coords->distance;
        $dih_val = $hyb eq 'sp2' ? 180 : 
            $len_ref->attr("3dbuilder/second_dihedral_ang") ?
            240 : 120;
        $len_ref->attr("3dbuilder/second_dihedral_ang", $atom->id);
    } else {
        # for fd(first_dihedralang): ang_ref = ok; dih_ref = dih_ref(fd)
        $dih_ref = first { $_->internal_coords } $ang_ref->neighbors($len_ref);
        if ($dih_ref) {
            $dih_val = 180;
            $ang_ref->attr("3dbuilder/first_dihedral_ang", $atom->id);
        } else {
            $dih_ref = first { $_->internal_coords and $_ ne $ang_ref } 
                $len_ref->neighbors($atom);
            $dih_val = $hyb eq 'sp2' ? 180 : 120;
        }
        $len_ref->attr("3dbuilder/first_dihedral", $atom->id);
    }
    
    ($len_ref, $len_val, $ang_ref, $ang_val, $dih_ref, $dih_val);
}

sub add_atom {
    my (%opts) = @_;
    my ($atom, $bond, $from, $user) = @opts{qw(atom bond from user)};
    my $stack = $user->{stack};
    my $n = ++$user->{n};

    #print "adding atom $n: $atom!\n" if $DEBUG;
    my $ic;
    if ($n == 1) {                             
        # first atom; place at origin
        $ic = Chemistry::InternalCoords->new($atom);
    } elsif ($n == 2) {
        # second atom only needs a bond
        $ic = Chemistry::InternalCoords->new(
            $atom, choose_bond($atom, $bond, $from));
    } elsif ($n == 3) {
        # third atom, needs an angle
        $ic = Chemistry::InternalCoords->new(
            $atom, choose_angle($atom, $bond, $from));
    } else {
        # other atoms need dihedral
        $ic = Chemistry::InternalCoords->new(
            $atom, choose_dihedral($atom, $bond, $from));
    }
    $atom->internal_coords($ic);
}

# do a depth-first traversal of a molecule. This should be added to 
# Chemistry::Mol

sub dfs {
    my ($mol, %opts) = @_;
    my ($on_atom, $on_bond, $user) = @opts{qw(on_atom on_bond user)};
    my $dfs;
    my %visited;
    $dfs = sub {
        my ($atom) = @_;
        $visited{$atom} = 1;
        for my $bn ($atom->bonds_neighbors) {
            my $nei  = $bn->{to};
            my $bond = $bn->{bond};

            next if $visited{$bond}++;
            $on_bond->(mol => $mol, from => $atom, to => $nei, bond => $bond,
                user => $user)
                if $on_bond;

            next if $visited{$nei};
            $on_atom->(mol => $mol, from => $atom, atom => $nei, 
                bond => $bond, user => $user) if $on_atom;
            $dfs->($nei);
        }
    };
    return unless $mol->atoms;
    $on_atom->(mol => $mol, atom => $mol->atoms(1), from => undef, 
        bond => undef, user => $user) if $on_atom;
    $dfs->($mol->atoms(1));
}

1;

=back

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::InternalCoords>.

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

