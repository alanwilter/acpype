package Chemistry::File::MidasPattern; 

$VERSION = '0.11';

=head1 NAME

Chemistry::File::MidasPattern - Wrapper Chemistry::File class for Midas patterns

=head1 SYNOPSIS

    use Chemistry::File::MidasPattern;
    use Chemistry::File::PDB;

    # read a molecule
    my $mol = Chemistry::MacroMol->read("test.pdb");

    # define a pattern matching carbons alpha and beta
    # in all valine residues
    my $str  = ':VAL@CA,CB';
    my $patt = Chemistry::MidasPattern->parse($str, format => 'midas');
    # Chemistry::Mol->parse($str, format => 'midas') also works

    # apply the pattern to the molecule
    $patt->match($mol);

    # extract the results
    for my $atom ($patt->atom_map) {
        printf "%s\t%s\n",  $atom->attr("pdb/residue_name"), $atom->name;
    }
    printf "FOUND %d atoms\n", scalar($patt->atom_map);

=head1 DESCRIPTION

This is a wrapper class for reading Midas Patterns using the standard
Chemistry::Mol I/O interface. This allows Midas patterns to be used
interchangeably with other pattern languages such as SMARTS in the context of
programs such as L<mok>. All of the real work is done by
L<Chemistry::MidasPattern>.

This module register the 'midas' format with Chemistry::Mol.

=cut

use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::MidasPattern;

Chemistry::Mol->register_format(midas => __PACKAGE__);

sub parse_string {
    my ($self, $s, %opts) = @_;
    my $patt = Chemistry::MidasPattern->new($s);
    $patt->options(%opts);
    $patt;
}

1;

=head1 VERSION

0.11

=head1 SEE ALSO

L<Chemistry::MidasPattern>, L<Chemistry::File>, L<Chemistry::Mol>,
L<Chemistry::MacroMol>, L<mok>.

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

