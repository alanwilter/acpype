package Chemistry::Domain;

$VERSION = '0.06';
# $Id: Domain.pm,v 1.6 2004/07/03 19:19:44 itubert Exp $


use 5.006;
use strict;
use warnings;
use base qw(Chemistry::Mol);
use Carp;
use Scalar::Util 'weaken';

=head1 NAME

Chemistry::Domain - Class for domains in macromolecules

=head1 SYNOPSIS

  use Chemistry::Domain;
  my $domain = Chemistry::Domain->new(parent => $bigmol);

=head1 DESCRIPTION

A domain is a substructure of a larger molecule. It is typically used to
represent aminoacid residues within a protein, or bases within a nucleic acid,
but you could use it for any arbitrary substructure such as functional groups
and rings. A domain has all the properties of a molecule, plus a "parent". The
parent is the larger molecule that contains the domain. In other words, the
Chemistry::Domain class inherits from Chemistry::Mol.

=head1 METHODS

Note: the methods that are inherited from Chemistry::Mol are not repeated here.

=over 4

=item Chemistry::Domain->new(parent => $mol, name => value, ...)

Create a new Domain object with the specified attributes. You can use the same 
attributes as for Chemistry::Mol->new, plus the parent attribute, which is 
required.

=cut

sub new {
    my $class = shift;
    my %args = @_;
    my $self = bless $class->SUPER::new(), $class;
    $self->$_($args{$_}) for (keys %args);
    croak "Must specify parent to Chemistry::Domain->new" unless $self->parent;
    return $self;
}


=item $domain->parent

Returns the parent of the domain.

=cut

sub parent {
    my $self = shift;
    if (@_) {
        $self->{parent} = shift;
        weaken($self->{parent});
        return $self;
    } else {
        return $self->{parent};
    }
}

=item $domain->add_atom($atom, ...)

Add one or more Atom objects to the domain. Returns the last atom added. It 
also automatically adds the atoms to the atom table of the parent molecule.

=cut

sub add_atom {
    my $self = shift;
    $self->add_atom_np(@_); # add atom to self (domain)
    $self->parent->add_atom(@_); # add atom to parent (macromol)
}


=item $domain->add_bond($bond, ...)

Add one or more Bond objects to the domain. Returns the last bond added. It 
also automatically adds the bond to the bond table of the parent molecule.

=cut

sub add_bond {
    my $self = shift;
    $self->add_bond_np(@_); # add bond to self (domain)
    $self->parent->add_bond(@_); # add bond to parent (macromol)
}


1;

=back

=head1 VERSION

0.06

=head1 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>, L<Chemistry::Atom>, 
L<Chemistry::Bond>

=head1 AUTHOR

Ivan Tubert, E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2004 by Ivan Tubert

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut

