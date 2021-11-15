package Chemistry::File::FormulaPattern; 

$VERSION = '0.10';

=head1 NAME

Chemistry::File::FormulaPattern - Wrapper Chemistry::File class for Formula patterns

=head1 SYNOPSIS

    use Chemistry::File::FormulaPattern;

    # somehow get a bunch of molecules...
    use Chemistry::File::SDF;
    my @mols = Chemistry::Mol->read("file.sdf");

    # we want molecules with six carbons and 8 or more hydrogens
    my $patt = Chemistry::Pattern->new("C6H8-", format => "formula_pattern");

    for my $mol (@mols) {
        if ($patt->match($mol)) {
            print $mol->name, " has a nice formula!\n";
        }
    }

    # a concise way of selecting molecules with grep
    my @matches = grep { $patt->match($mol) } @mols;

=head1 DESCRIPTION

This is a wrapper class for reading Formula Patterns using the standard
Chemistry::Mol I/O interface. This allows Formula patterns to be used
interchangeably with other pattern languages such as SMARTS in the context of
programs such as L<mok>. All of the real work is done by
L<Chemistry::FormulaPattern>.

This module register the 'formula_pattern' format with L<Chemistry::Mol>.

=cut

use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::FormulaPattern;

Chemistry::Mol->register_format(formula_pattern => __PACKAGE__);

sub parse_string {
    my ($self, $s, %opts) = @_;
    my $patt = Chemistry::FormulaPattern->new($s);
    #$patt->options(%opts);
    $patt;
}

1;

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Pattern>, L<Chemistry::File>, L<Chemistry::Mol>,
L<Chemistry::MacroMol>, L<mok>.

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

