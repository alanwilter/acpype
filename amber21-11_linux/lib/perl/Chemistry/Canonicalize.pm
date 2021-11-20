package Chemistry::Canonicalize;

# $Id: Canonicalize.pm,v 1.1.1.1 2004/06/16 19:40:37 ivan Exp $
$VERSION = '0.10';

use strict;
use warnings;

use Carp;
use base 'Exporter';
our @EXPORT_OK = qw(canonicalize);
our %EXPORT_TAGS = ( all => \@EXPORT_OK ); 


=head1 NAME

Chemistry::Canonicalize - Number the atoms in a molecule in a unique way

=head1 SYNOPSIS

    use Chemistry::Canonicalize ':all';

    # $mol is a Chemistry::Mol object
    canonicalize($mol);
    print "The canonical number for atom 1 is: ", 
        $mol->atoms(1)->attr("canon/class");
    print "The symmetry class for for atom 1 is: ", 
        $mol->atoms(1)->attr("canon/symmetry_class");

=head1 DESCRIPTION

This module provides functions for "canonicalizing" a molecular structure; that
is, to number the atoms in a unique way regardless of the input order.

The canonicalization algorithm is based on: Weininger, et. al., J. Chem. Inf.
Comp. Sci. 29[2], 97-101 (1989)

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 ATOM ATTRIBUTES

During the canonicalization process, the following attributes are set on each
atom:

=over

=item canon/class

The unique canonical number; it is an integer going from 1 to the number of
atoms.

=item canon/symmetry_class

The symmetry class number. Atoms that have the same symmetry class are
considered to be topologicaly equivalent. For example, the two methyl carbons
on 2-propanol would have the same symmetry class.

=back

=head1 FUNCTIONS

These functions may be exported, although nothing is exported by default.

=over

=cut

my @PRIMES = qw( 1
      2      3      5      7     11     13     17     19     23     29 
     31     37     41     43     47     53     59     61     67     71 
     73     79     83     89     97    101    103    107    109    113 
    127    131    137    139    149    151    157    163    167    173 
    179    181    191    193    197    199    211    223    227    229 
    233    239    241    251    257    263    269    271    277    281 
    283    293    307    311    313    317    331    337    347    349 
    353    359    367    373    379    383    389    397    401    409 
    419    421    431    433    439    443    449    457    461    463 
    467    479    487    491    499    503    509    521    523    541 
    547    557    563    569    571    577    587    593    599    601 
    607    613    617    619    631    641    643    647    653    659 
    661    673    677    683    691    701    709    719    727    733 
    739    743    751    757    761    769    773    787    797    809 
    811    821    823    827    829    839    853    857    859    863 
    877    881    883    887    907    911    919    929    937    941 
    947    953    967    971    977    983    991    997   1009   1013 
   1019   1021   1031   1033   1039   1049   1051   1061   1063   1069 
   1087   1091   1093   1097   1103   1109   1117   1123   1129   1151 
   1153   1163   1171   1181   1187   1193   1201   1213   1217   1223 
);

=item canonicalize($mol, %opts)

Canonicalizes the molecule. It adds the canon/class and canon/symmetry class to
every atom, as discussed above. This function may take the following options:

=over

=item sort

If true, sort the atoms in the molecule in ascending canonical number order.

=item invariants

This should be a subroutine reference that takes an atom and returns a number.
These number should be based on the topological invariant properties of the
atom, such as symbol, charge, number of bonds, etc.

=back

=cut

sub canonicalize {
    my ($mol, %opts) = @_;

    if ($mol->atoms > @PRIMES - 1) {
        croak "maximum number of atoms exceeded for canonicalization\n";
    }

    my $invariants_sub = $opts{invariants} || \&atom_invariants;

    # set up initial classes
    for my $atom ($mol->atoms) {
        $atom->attr("canon/class", $invariants_sub->($atom));
        $atom->attr("canon/prev_class", 1);
    }
    #printf "$_: %s\n", $_->attr("canon/class") for $mol->atoms;

    # run one canonicalization step
    my $atoms;
    my $n_classes;
    ($atoms, $n_classes) = rank_classes($mol);
    ($atoms, $n_classes) = canon($mol, $n_classes);
    my $n_atom = $mol->atoms;

    # atoms with the same class are topologically symmetric
    for my $atom ($mol->atoms) {
        $atom->attr("canon/symmetry_class", $atom->attr("canon/class"));
    }
    #printf "$_: %s\n", $_->attr("canon/class") for $mol->atoms;

    # break symmetry to get a canonical numbering
    while ($n_classes < $n_atom) {

        # multiply all classes by 2
        for my $atom (@$atoms) {
            my $class = $atom->attr("canon/class");
            $atom->attr("canon/class", $class * 2); 
        }

        # break first tie
        my $last_class = -1;
        my $last_atom;
        for my $atom (@$atoms) {
            my $class = $atom->attr("canon/class");
            if ($class == $last_class) { # tie
                #print "breaking tie for $last_atom\n";
                $last_atom->attr("canon/class", $class - 1);
                last;
            }
            $last_class = $class;
            $last_atom  = $atom;
        }
        #printf "$_: %s\n", $_->attr("canon/class") for $mol->atoms;
        #print "---\n";

        # run another canonicalization step
        ($atoms, $n_classes) = canon($mol, $n_classes);
        #printf "$_: %s\n", $_->attr("canon/class") for $mol->atoms;
    }
    if ($opts{'sort'}) {
        $mol->sort_atoms( 
            sub { $_[0]->attr("canon/class") <=> $_[1]->attr("canon/class") } 
        );
    }
    # clean up temporary classes
    $_->del_attr("canon/new_class") for $mol->atoms;
    $n_classes;
}

sub atom_invariants {
    no warnings 'uninitialized';
    my ($atom) = @_;
    my $n_bonds = $atom->bonds;
    my $valence = 0;
    $valence += $_->order for $atom->bonds;
    my $Z = $atom->Z;
    my $q = $atom->formal_charge + 5;
    return $n_bonds*10_000 + $valence*1000 + $q*100 + $Z;
}

# atom class comparison function. Only compare the class if the 
# previous classes are equal
sub _cmp {
    $a->attr("canon/prev_class") <=> $b->attr("canon/prev_class")
    or $a->attr("canon/class") <=> $b->attr("canon/class")
}

sub rank_classes {
    my ($mol) = @_;
    my @atoms = sort _cmp $mol->atoms; # consider Schwartzian transform?
    my $n = 0;
    local ($a, $b);
    for $b (@atoms) {
        $n++ if (!$a || _cmp);
        $a = $b;
        $b->attr("canon/new_class", $n);
    }
    #use diagnostics;
    for my $atom (@atoms) {
        $atom->attr("canon/class", $atom->attr("canon/new_class"));
    }
    (\@atoms, $n);
}

sub canon {
    my ($mol, $n) = @_;

    my $old_classes = 0;
    my $n_atom = $mol->atoms;
    my $atoms;
    while ($n > $old_classes and $n < $n_atom) {
        $old_classes = $n;
        # save current classes
        for my $atom ($mol->atoms) {
            $atom->attr("canon/prev_class", $atom->attr("canon/class"));
        }

        # set new class to product of neighbor's primes
        for my $atom ($mol->atoms) {
            my $class = 1;
            for my $neighbor ($atom->neighbors) {
                $class *= $PRIMES[$neighbor->attr("canon/prev_class")];
            }
            $atom->attr("canon/class", $class);
        }
        ($atoms, $n) = rank_classes($mol);
    }
    ($atoms, $n);
}

1;

=back

=head1 VERSION

0.10

=head1 TO DO

Add some tests.

=head1 CAVEATS

Currently there is an atom limit of 200 atoms.

These algorithm is known to fail to discriminate between non-equivalent atoms
for some complicated cases. These are usually highly bridged structures
explicitly designed to break canonicalization algorithms; I don't know of any
"real-looking structure" (meaning something that someone would actually
synthesize or find in nature) that fails, but don't say I didn't warn you!

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Obj>,
L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

