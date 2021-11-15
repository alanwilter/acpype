package Chemistry::Pattern::Atom;
$VERSION = '0.26';
# $Id: Atom.pm,v 1.11 2005/05/16 22:08:15 itubert Exp $

=head1 NAME

Chemistry::Pattern::Atom - An atom that knows how to match

=head1 SYNOPSIS

    my $patt_atom = Chemistry::Pattern::Atom->new(symbol => C);
    $patt_atom->test_sub( sub {
        my ($what, $where) = @_; 
        $where->bonds == 3 ? 1 : 0; # only match atoms with three bonds
    });

=head1 DESCRIPTION

Objects of this class represent atoms in a pattern. This is a subclass of
Chemistry::Atom. In addition to the properties of regular atoms, 
pattern atoms have a method for testing if they match an atom in a molecule.
By default, a pattern atom matches an atom if they have the same symbol.
It is possible to substitute this by an arbitrary criterion by providing
a custom test subroutine.

=cut

use 5.006;
use strict;
use Carp;
use base qw(Chemistry::Atom);

=head1 METHODS

=over 4

=cut

our $Debug = 0;

=item $patt_atom->test($atom)

Tests if the pattern atom matches the atom given by $atom. Returns true or
false.

=cut

sub test {
    my ($what, $where) = @_;
    #print "\t\ttesting $where against $what\n";
    #print $where->print, $what->print;
    if ($what->test_sub) {
        #print "\t\thave a test sub\n";
        return $what->test_sub->($what, $where);
    } else {
        #print "\t\tdon't have a test sub\n";
        return $what->symbol eq $where->symbol;
    }
}

=item $patt_atom->test_sub(\&my_test_sub)

Specify an arbitrary test subroutine to be used instead of the default one.
&my_test_sub must take two parameters; the first one is the pattern atom 
and the second is the atom to match. It must return true if there is a match.

=cut

Chemistry::Obj::accessor('test_sub');


=item $patt_atom->map_to([$atom])

Returns or sets the atom that is considered to be matched by $patt_atom.

=cut

Chemistry::Obj::accessor('map_to');

1;

=back

=head1 VERSION

0.26

=head1 SEE ALSO

L<Chemistry::Pattern>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

