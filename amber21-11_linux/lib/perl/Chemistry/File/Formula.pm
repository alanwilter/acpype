package Chemistry::File::Formula;
$VERSION = '0.35';
# $Id: Formula.pm,v 1.13 2005/05/20 19:01:05 itubert Exp $

use strict;
use base "Chemistry::File";
use Chemistry::Mol;
use Carp;
use Text::Balanced qw(extract_bracketed);

=head1 NAME

Chemistry::File::Formula - Molecular formula reader/formatter

=head1 SYNOPSIS

    use Chemistry::File::Formula;

    my $mol = Chemistry::Mol->parse("H2O");
    print $mol->print(format => formula);
    print $mol->formula;    # this is a shorthand for the above 
    print $mol->print(format => formula, 
        formula_format => "%s%d{<sub>%d</sub>});

=cut

Chemistry::Mol->register_format('formula');

=head1 DESCRIPTION

This module converts a molecule object to a string with the formula. It
registers the 'formula' format with Chemistry::Mol.
Besides its obvious use, it is included in the Chemistry::Mol distribution
because it is a very simple example of a Chemistry::File derived I/O module.

The format can be specified as a printf-like string with the following control
sequences, which are specified with the formula_format parameter to $mol->print
or $mol->write.

=over

=item %s  symbol

=item %D  number of atoms

=item %d  number of atoms, included only when it is greater than one

=item %d{substr}  substr is only included when number of atoms is greater than
one

=item %j{substr}  substr is inserted between the formatted string for each
element. (The 'j' stands for 'joiner'). The format should have only one joiner,
and its location in the string doesn't matter.

=item %% a percent sign

=back

If no format is specified, the default is "%s%d". Some examples follow. Let's
assume that the formula is C2H6O, as it would be formatted by default.

=over

=item C<< %s%D >>

Like the default, but include explicit indices for all atoms. 
The formula would be formatted as "C2H6O1"

=item C<< %s%d{E<lt>subE<gt>%dE<lt>/subE<gt>} >>

HTML format. The output would be
"CE<lt>subE<gt>2E<lt>/subE<gt>HE<lt>subE<gt>6E<lt>/subE<gt>O".

=item C<< %D %s%j{, } >>

Use a comma followed by a space as a joiner. The output would be 
"2 C, 6 H, 1 O".

=back

=head2 Parsing Formulas

Formulas can also be parsed back into Chemistry::Mol objects.
The formula may have parentheses and square or triangular brackets, and 
it may have the following abbreviations:

    Me => '(CH3)',
    Et => '(CH3CH2)',
    Bu => '(C4H9)',
    Bn => '(C6H5CH2)',
    Cp => '(C5H5)',
    Ph => '(C6H5)',
    Bz => '(C6H5CO)',

The formula may also be preceded by a number, which multiplies the whole
formula. Some examples of valid formulas:

=over

=item CH3(CH2)3CH3. Equivalent to C5H12.

=item C6H3Me3. Equivalent to C9H12.

=item 2Cu[NH3]4(NO3)2. Equivalent to Cu2H24N12O12.

=item 2C(C[C<C>5]4)3. Equivalent to C152 (kind of silly example...).

=item 2C(C(C(C)5)4)3. Equivalent to C152.

=item C 1 0 H 2 2. Equivalent to C10H22. Note that whitespace is completely
ignored.

=back

When a formula is parsed, a molecule object is created which consists of
the set of the atoms in the formula (no bonds or coordinates, of course).
The atoms are created in alphabetical order, so the molecule object for C2H5Br
would have the atoms in the following sequence: Br, C, C, H, H, H, H, H. 

=cut

sub parse_string {
    my ($self, $string, %opts) = @_;
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || "Chemistry::Atom";
    my $bond_class = $opts{bond_class} || "Chemistry::Bond";

    my $mol = $mol_class->new;
    #$string =~ /^(?:([A-Z][a-z]*)(\d*))+$/ or croak("invalid formula $string\n");
    #my (%formula) = $string =~ m/([A-Z][a-z]*)(\d*)/g;
    #for (values %formula) {
        #$_ = 1 unless length; # Add implicit indices
    #}
    my %formula = parse_formula($string);
    for my $sym (sort keys %formula) {
        for (my $i = 0; $i < $formula{$sym}; ++$i) {
            $mol->add_atom($atom_class->new(symbol => $sym));
        }
    }
    return $mol;
}

sub write_string {
    my ($self, $mol, %opts) = @_;
    my @formula_parts;
    my $format = $opts{formula_format} || "%s%d";   # default format
    my $fh = $mol->formula_hash;
    $format =~ s/%%/\\%/g;                          # escape %% with a \
    my $joiner = "";
    $joiner = $1 if $format =~ s/(?<!\\)%j{(.*?)}//;        # joiner %j{}
    for my $sym (sort keys %$fh) {
        my $s = $format;
        my $n = $fh->{$sym};
        $s =~ s/(?<!\\)%s/$sym/g;                           # %s
        $s =~ s/(?<!\\)%D/$n/g;                             # %D
        $s =~ s/(?<!\\)%d\{(.*?)\}/$n > 1 ? $1 : ''/eg;     # %d{}
        $s =~ s/(?<!\\)%d/$n > 1 ? $n : ''/eg;              # %d
        $s =~ s/\\(.)/$1/g;                                 # other \ escapes
        push @formula_parts, $s;
    }
    return join($joiner, @formula_parts);
}


sub file_is {
    return 0; # no files are identified automatically as having this format
}

### Code derived from formula.pl by Brent Gregersen follows

my %macros = (
    Me => '(CH3)',
    Et => '(CH3CH2)',
    Bu => '(C4H9)',
    Bn => '(C6H5CH2)',
    Cp => '(C5H5)',
    Ph => '(C6H5)',
    Bz => '(C6H5CO)',
    # Ac is an element
    # Pr is an element
);


sub parse_formula {
    my ($formula) = @_;
    my (%elements);

    #check balancing
    return %elements if (!ParensBalanced($formula));

    # replace other grouping with normal parens
    $formula =~ tr/<>{}[]/()()()/;

    # get rid of any spaces
    $formula =~ s/\s+//g;

    # perform macro expansion
    foreach (keys(%macros)) {
        $formula =~ s/$_/$macros{$_}/g;
    }

    # determine initial compound coeficent
    my $coef = ($formula =~ s/^(\d+)//) ? $1 : 1.0;

    # recursively process rest of formula
    return internal_formula_parser($formula, $coef, %elements);
}

sub internal_formula_parser {
    my ($formula, $coef, %form) = @_;
    my $tmp_coef;

    my ($extract, $remainder, $prefix) =
      extract_bracketed($formula, '()', '[^(]*');

    if (defined($extract) and $extract ne '') {
        $extract =~ s/^\((.*)\)$/$1/;
        if ($remainder =~ s/^(\d+)(.*)$/$2/) {
            $tmp_coef = $1 * $coef;
        } else {
            $tmp_coef = $coef;
        }

        # get formula of prefix ( it has no parens)
        %form = add_formula_strings($prefix, $coef, %form) if ($prefix ne '');

        # check remainder for more parens
        %form = internal_formula_parser($remainder, $coef, %form)
          if ($remainder ne '');

        # check extract for more parens
        %form =
          internal_formula_parser($extract, $tmp_coef, %form);    
          ## we already know this is ne ''
    } else { # get formula of complete string
        %form = add_formula_strings($remainder, $coef, %form)
          if ($remainder ne '');
    }
    return %form;
}

sub add_formula_strings {
    my ($formula, $coef, %elements) = @_;

#  print "Getting Formula of $formula\n";
    $formula =~ /^(?:([A-Z][a-z]*)(\d+)?)+$/o
        or croak "Invalid Portion of Formula $formula";
    while ($formula =~ m/([A-Z][a-z]*)(\d+)?/go) {
        my ($elm, $count) = ($1, $2);
        $count = 1 unless defined $count;
        if (defined $elements{$elm}) {
            $elements{$elm} += $count * $coef;
        } else {
            $elements{$elm} = $count * $coef;
        }
    }
    return %elements;
}

sub ParensBalanced {
    my ($form) = @_;
    my @stack  = ();
    my %pairs  = (
        '<' => '>',
        '{' => '}',
        '[' => ']',
        '(' => ')'
    );

    while ($form =~ m/([<>(){}\]\[])/go) { 
        my $current = $1;
        if ($current =~ /[<({\[]/) {
            push(@stack, $current);
            next;
        }
        return 0 if (scalar(@stack) == 0);
        return 0 if ($current ne $pairs{ pop @stack});
    }
    return @stack ? 0 : 1;
}

1;

=head1 VERSION

0.35

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>.

Formula parsing code contributed by Brent Gregersen.

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

