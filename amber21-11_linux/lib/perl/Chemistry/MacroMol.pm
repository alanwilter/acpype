package Chemistry::MacroMol;

$VERSION = '0.06';
# $Id: MacroMol.pm,v 1.6 2004/07/03 19:19:44 itubert Exp $

use 5.006;
use strict;
use warnings;
use base qw(Chemistry::Mol);

=head1 NAME

Chemistry::MacroMol - Perl module for macromolecules

=head1 SYNOPSIS

    use Chemistry::MacroMol;

    my $mol = Chemistry::MacroMol->new(name => 'my big molecule');
    $mol->new_domain(name => "ASP"); # see Chemistry::Domain for details
    my @domains = $mol->domains;

=head1 DESCRIPTION

For the purposes of this module, a macromolecule is just a molecule that 
consists of several "domains". For example, a protein consists of aminoacid
residues, or a nucleic acid consists of bases. Therefore Chemistry::MacroMol 
is derived from Chemistry::Mol, with additional methods to handle the domains.

The way things are currently structured, an atom in a macromolecule "belong"
both to the MacroMol object and to a Domain object. This way you can get all the
atoms in $protein via $protein->atoms, or to the atoms in residue 123 via
$protein->domain(123)->atoms.

=head1 METHODS

Remember that this class inherits all the methods from Chemistry::Mol. They
won't be repeated here.

=over 4

=item Chemistry::MacroMol->new(name => value, ...)

Create a new MacroMol object with the specified attributes. You can use the
same attributes as for Chemistry::Mol->new.

=cut

sub new {
    my $class = shift;
    my %args = @_;
    my $self = bless $class->SUPER::new(), $class;
    $self->{domains} = [];
    $self->$_($args{$_}) for (keys %args);
    return $self;
}

=item $mol->add_domain($domain, ...)

Add one or more Domain objects to the molecule. Returns the last domain added.

=cut

sub add_domain {
    my $self = shift;
    for my $b (@_){
        push @{$self->{domains}}, $b;
	$self->{byId}{$b->{id}} = $b;
    }
    $_[-1];
}

=item $mol->domain_class

Returns the domain class that a macromolecule class expects to use by default.
Chemistry::MacroMol objects return "Chemistry::Domain", but subclasses will
likely override this method.

=cut

sub domain_class { "Chemistry::Domain" }

=item $mol->new_domain(name => value, ...)

Shorthand for $mol->add_domain($mol->domain_class->new(parent => $mol, name => value, ...));

=cut

sub new_domain {
    my $self = shift;
    $self->add_domain(Chemistry::Domain->new(parent => $self, @_));
}

=item $mol->domains($n1, ...)

Returns the domains with the given indices, or all by default. NOTE:
the indices start from one (1), not from zero.

=cut

sub domains {
    my $self = shift;
    if (@_) {
        my @doms = map {$_ - 1} @_;
        wantarray ? @{$self->{domains}}[@doms] : $self->{domains}[$doms[-1]];
    } else {
        @{$self->{domains}}; # return all
    }
}


1;

=back

=head1 VERSION

0.06

=head1 SEE ALSO

L<Chemistry::Domain>, L<Chemistry::Mol>

=head1 AUTHOR

Ivan Tubert, E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2004 by Ivan Tubert

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut

