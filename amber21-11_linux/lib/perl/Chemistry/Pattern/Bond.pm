package Chemistry::Pattern::Bond;
$VERSION = '0.26';
# $Id: Bond.pm,v 1.10 2005/05/16 22:08:15 itubert Exp $

=head1 NAME

Chemistry::Pattern::Bond - A bond that knows how to match

=head1 SYNOPSIS

    my $patt_bond = Chemistry::Pattern::Bond->new(order => 2);
    $patt_bond->test_sub( sub {
        my ($what, $where) = @_; 
        $where->type eq 'purple' ? 1 : 0; # only match purple bonds
    });

=head1 DESCRIPTION

Objects of this class represent bonds in a pattern. This is a subclass of
Chemistry::Bond. In addition to the properties of regular bonds, 
pattern bonds have a method for testing if they match an bond in a molecule.
By default, a pattern bond matches an bond if they have the same bond order.
It is possible to substitute this by an arbitrary criterion by providing
a custom test subroutine.

=cut

use 5.006;
use strict;
use base qw(Chemistry::Bond);

=head1 METHODS

=over 4

=cut

=item $patt_bond->test($bond)

Tests if the pattern bond matches the bond given by $bond. Returns true or
false.

=cut

sub test {
    my ($what, $where) = @_;
    if ($what->test_sub) {
         return $what->test_sub->($what, $where);
    } else {
         return $what->order eq $where->order;
    }
}



=item $patt_bond->test_sub(\&my_test_sub)

Specify an arbitrary test subroutine to be used instead of the default one.
&my_test_sub must take two parameters; the first one is the pattern bond 
and the second is the bond to match. It must return true if there is a match.

=cut

Chemistry::Obj::accessor('test_sub');

=item $patt_bond->map_to([$bond])

Returns or sets the bond that is considered to be matched by $patt_bond.

=cut

#Chemistry::Obj::accessor('map_to');
sub map_to {
    my $self = shift;
    if (@_) {
        #print "\t\tmapping $self to '@_'\n";
        ($self->{map_to}) = @_;
        $self;
    } else {
        #print "\t\t$self is mapped to '$self->{map_to}'\n";
        $self->{map_to};
    }
}

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

