package Chemistry::Ring::Find;

$VERSION = 0.19;
# $Id: Find.pm,v 1.1.1.1 2005/03/29 23:57:36 itubert Exp $

=head1 NAME

Chemistry::Ring::Find - Find the rings (cycles) in a molecule

=head1 SYNOPSIS

    use Chemistry::Ring::Find ':all';

    # find the smallest ring containing $atom
    my $ring = find_ring($atom);

    # find all the rings containing $bond
    my @rings = find_ring($bond, all => 1);

    # see below for more options

    # find the six 4-atom rings in cubane
    @rings = find_rings($cubane);

    # find a cubane SSSR with five rings
    @rings = find_rings($cubane, sssr => 1);

=head1 DESCRIPTION

The Chemistry::Ring::Find module implements a breadth-first ring finding
algorithm, and it can find all rings that contain a given atom or bond, the
Smallest Set of Smallest Rings (SSSR), or the "almost SSSR", which is an
unambiguous set of rings for cases such as cubane.The algorithms are  based on
ideas from:

1) Leach, A. R.; Dolata, D. P.; Prout, P. Automated Conformational Analysis and
Structure Generation: Algorithms for Molecular Perception J. Chem. Inf. Comput.
Sci. 1990, 30, 316-324

2) Figueras, J. Ring perception using breadth-first search. J. Chem. Inf.
Comput.  Sci. 1996, 36, 986-991.

Ref. 2 is only used for find_ring, not for find_rings, because it has been
shown that the overall SSSR method in ref 2 has bugs. Ref 1 inspired
find_rings, which depends on find_ring.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 FUNCTIONS

These functions may be exported explicitly, or all by using the :all tag, but
nothing is exported by default.

=over 

=cut


use strict;
use warnings;
use Chemistry::Ring;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(find_ring find_rings);
our %EXPORT_TAGS = ( all => \@EXPORT_OK ); 

our $DEBUG = 0;

=item find_ring($origin, %opts)

Find the smallest ring containg $origin, which may be either an atom or a bond.
Returns a Chemistry::Ring object. Options:

=over

=item all

If true, find all the rings containing $origin. If false, return the first ring
found. Defaults to false. "All" is supposed to include only "simple" rings,
that is, rings that are not a combination of smaller rings. 

=item min

Only find rings with a the given minimum size. Defaults to zero.

=item max

Only find rings up to the given maximium size. Defaults to unlimited size.

=item size

Only find rings with this size. Same as setting min and max to the same size.
Default: unspecified.

=item exclude

An array reference containing a list of atoms that must NOT be present in the
ring. Defaults to the empty list.

=item mirror

If true, find each ring twice (forwards and backwards). Defaults to false.

=back

=cut

# $origin is an atom
# options: min, max, size, all, mirror, exclude
sub find_ring {
    my ($origin, %opts) = @_;
    my $min_size = $opts{min} || $opts{size} || 0;
    my $max_size = $opts{max} || $opts{size};
    my %paths;
    my %bond_paths;
    my @q;
    my @rings;
    my %used_end_nodes;
    my $required_bond;
    my %exclude;
    @exclude{ @{$opts{exclude} || []} } = ();

    if ($origin->isa("Chemistry::Bond")) {
        $required_bond = $origin;
        ($origin) = $origin->atoms;
    }
    @q = ($origin);
    $paths{$origin} = [$origin];
    $bond_paths{$origin} = [];
    # $path{$atom} means how to get to $atom from $origin

    my $a;
    while ($a = shift @q) {
        my $from = $paths{$a}[-2];
        print "at $a from $from\n" if $DEBUG;
        for my $bn ($a->bonds_neighbors($from)) {
            my $nei  = $bn->{to};
            my $bond = $bn->{bond};
            next if exists $exclude{$nei};
            print "  -> $nei\n" if $DEBUG;
            if ($paths{$nei}) {
                print "found a path collision... " if $DEBUG;
                # check to make sure that the ring really started at $origin
                # and the size is what was requested
                my $size = @{$paths{$nei}} + @{$paths{$a}} - 1;
                if($paths{$nei}[1] != $paths{$a}[1]
                    and $size >= $min_size
                    and !$max_size || $size <= $max_size)
                {
                    print "VALID\n" if $DEBUG;
                    my @atoms = (@{$paths{$a}}, reverse @{$paths{$nei}});
                    print "RING = ", print_path(\@atoms) if $DEBUG;
                    pop @atoms;
                    if ($used_end_nodes{$atoms[1]} and !$opts{mirror}) {
                        print "skipping redundant ring\n" if $DEBUG;
                        next; # don't want to find rings twice
                    }
                    my @bonds = (@{$bond_paths{$a}}, $bond,
                        reverse @{$bond_paths{$nei}});
                    if ($required_bond 
                        and not grep {$_ eq $required_bond} @bonds) {
                        print "does not include required bond\n" if $DEBUG;
                        next;
                    }
                    if (contains_ring(\@atoms, \@rings)) {
                        print "contains another ring\n" if $DEBUG;
                        next;
                    }
                    my $r = Chemistry::Ring->new;
                    $r->add_atom(@atoms);
                    $r->add_bond(@bonds);
                    return $r unless $opts{all};  # FOUND VALID RING
                    push @rings, $r;
                    $used_end_nodes{$atoms[-1]} = 1;
                    #@used_nodes{@atoms} = ();
                } else {
                    print "NOT VALID", 
                        print_path( [@{$paths{$a}}, 
                            reverse @{$paths{$nei->id}}]) if $DEBUG;
                }
            } else {
                if (!$max_size || @{$paths{$a}} < ($max_size / 2) + 0.1) {
                    push @q, $nei;
                    print "    pushing path\n" if $DEBUG;        
                    $paths{$nei} = [@{$paths{$a}}, $nei];
                    $bond_paths{$nei} = [@{$bond_paths{$a}}, $bond];
                    print print_path($paths{$nei}) if $DEBUG;
                } else {
                    print "path too long; " if $DEBUG;
                    print print_path($paths{$a}) if $DEBUG;
                    #path too long
                }
            }
        }
    }
    @rings;
}


sub print_path {
    my $p = shift;
    my $ret = "    PATH: ";
    for my $a (@$p) {
        $ret .=  "$a - ";
    }
    $ret .= "\n";
}

# contains_ring($atoms, $rings)
# returns true if one of the rings in the array ref $rings is a proper subset
# of the atom list in the array ref $atom.
sub contains_ring {
    my ($atoms, $rings) = @_;
    my %seen;
    @seen{@$atoms} = ();
    for my $ring (@$rings) {
        my $unique_atoms = $ring->atoms; 
        next if $unique_atoms >= @$atoms; # ring is same size or bigger
        # make sure that $ring has at least one atom not in $atoms
        for my $atom ($ring->atoms) {
            if (exists $seen{$atom}) {
                $unique_atoms--;
            } else {
                last; # there's at least one unique atom!
            }
        }
        return 1 unless $unique_atoms;
    }
    0;
}

=item @rings = find_rings($mol, %options)

Find "all" the rings in the molecule. In general it return the Smallest Set of
Smallest Rings (SSSR). However, since it is well known that the SSSR is not
unique for molecules such as cubane (where the SSSR consists of five
unspecified four-member rings, even if the symmetry of the molecule would
suggest that the six faces of the cube are equivalent), in such cases
find_rings will return a non-ambiguous "non-smallest" set of smallest rings,
unless the "sssr" option is given. For example,

    @rings = find_rings($cubane);
    # returns SIX four-member rings

    @rings = find_rings($cubane, sssr => 1);
    # returns FIVE four-member rings (an unspecified subset of
    # the six rings above.)

=cut

sub find_rings {
    my ($mol, %opts) = @_;
    my $visited = {};
    my @ring_bonds;
    for my $atom ($mol->atoms) {
        next if $visited->{$atom};
        push @ring_bonds, find_ring_bonds($mol, \%opts, $atom, $visited);
    }
    my @rings;
    my $n_rings = @ring_bonds;
    #print "cyclomatic number=$n_rings\n";
    for my $ring_bond (@ring_bonds) {
        push @rings, find_ring($ring_bond, all => 1)
    }
    my %seen;
    my @ring_keys = map { join " ", sort $_->atoms } @rings;
    @seen{@ring_keys} = @rings;
    @rings = sort { $a->atoms <=> $b->atoms } values %seen;
    if (!$opts{sssr} and @rings > $n_rings) {
        $n_rings++ if $rings[$n_rings]->atoms == $rings[$n_rings-1]->atoms;
    }
    splice @rings, $n_rings;
    @rings;
}

sub find_ring_bonds {
    my ($mol, $opts, $atom, $visited) = @_;

    my @ring_bonds;
    $visited->{$atom}  = 1;
    for my $bn ($atom->bonds_neighbors) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond}  = 1;
        if ($visited->{$nei}) { # closed ring
            #print "closing ring\n";
            push @ring_bonds, $bond;
        } else {
            push @ring_bonds, 
                find_ring_bonds($mol, $opts, $nei, $visited);
        }
    }
    @ring_bonds;
}

1;

=back

=head1 BUGS

The "all" option in find_ring doesn't quite work as expected. It finds all
simple rings and some bridged rings. It never finds fused rings (which is good).

=head1 VERSION

0.19

=head1 SEE ALSO

L<Chemistry::Ring>, L<http://www.perlmol.org>.

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

