package Chemistry::FormulaPattern;

$VERSION = '0.10';

# $Id: FormulaPattern.pm,v 1.1.1.1 2004/08/11 23:03:04 itubert Exp $

=head1 NAME

Chemistry::FormulaPattern - Match molecule by formula

=head1 SYNOPSIS

    use Chemistry::FormulaPattern;

    # somehow get a bunch of molecules...
    use Chemistry::File::SDF;
    my @mols = Chemistry::Mol->read("file.sdf");

    # we want molecules with six carbons and 8 or more hydrogens
    my $patt = Chemistry::FormulaPattern->new("C6H8-");

    for my $mol (@mols) {
        if ($patt->match($mol)) {
            print $mol->name, " has a nice formula!\n";
        }
    }

    # a concise way of selecting molecules with grep
    my @matches = grep { $patt->match($mol) } @mols;

=head1 DESCRIPTION

This module implements a simple language for describing a range of
molecular formulas and allows one to find out whether a molecule matches
the formula specification. It can be used for searching for molecules by
formula, in a way similar to the NIST WebBook formula search
(L<http://webbook.nist.gov/chemistry/form-ser.html>). Note however that the
language used by this module is different from the one used by the WebBook!

Chemistry::FormulaPattern shares the same interface as L<Chemistry::Pattern>.
To perform a pattern matching operation on a molecule, follow these steps.

1) Create a pattern object, by parsing a string. Let's assume that the pattern
object is stored in $patt and that the molecule is $mol.

2) Execute the pattern on the molecule by calling $patt->match($mol).

If $patt->match returns true, there was a match. If $patt->match is called two
consecutive times with the same molecule, it returns false; then true (if there
is a match), then false, etc. This is because the Chemistry::Pattern interface
is designed to allow multiple matches for a given molecule, and then returns
false when there are no further matches; in the case of a formula pattern,
there is only one possible match.

    $patt->match($mol); # may return true
    $patt->match($mol); # always false
    $patt->match($mol); # may return true
    $patt->match($mol); # always false
    # ...

This allows one two use the standard while loop for all kinds of patterns
without having to worry about endless loops:

    # $patt might be a Chemistry::Pattern, Chemistry::FormulaPattern,
    # or Chemistry::MidasPattern object
    while ($patt->match($mol)) {
        # do something
    }

Also note that formula patterns don't really have the concept of an atom map,
so $patt->atom_map and $patt->bond_map always return the empty list.

=head1 FORMULA PATTERN LANGUAGE

In the simplest case, a formula pattern may be just a regular formula, as
used by the L<Chemistry::File::Formula> module. For example, the pattern
"C6H6" will only match molecules with six carbons, six hydrogens, and no other
atoms. 

The interesting thing is that one can also specify ranges for the elements, 
as two hyphen-separated numbers. "C6H8-10" will match molecules with six
carbons and eight to ten hydrogens.

Ranges may also be open, by omitting the upper part of the range. "C6H0-" will
match molecules with six carbons and any number of hydrogens (i.e., zero or
more).

A formula pattern may also allow for unspecified elements by means of the
asterisk special character, which can be placed anywhere in the formula
pattern. For example, "C2H6*" (or "C2*H6, etc.) will match C2H6, and also
C2H6O, C2H6S, C2H6SO, etc.

Ranges can also be used after a subformula in parentheses: "(CH2)1-2" will
match molecules with one or two carbons and two to four hydrogens. Note,
however, that the "structure" of the bracketed part of the formula is
forgotten, i.e., the multiplier applies to each element individually and does
not have to be an integer. That is, the above pattern will match CH2, CH3, CH4,
C2H2, C2H3, and C2H4.

=cut


use strict;
use warnings;
use Chemistry::Mol;
use Carp;
use base "Chemistry::Pattern";
use Text::Balanced qw(extract_bracketed);

our $DEBUG = 0;

sub new {
    my ($class, $str) = @_;
    my $self = bless { 
        atom_map => [],
        options  => {},
        matched  => 0,
    }, ref $class || $class;
    $self->parse_formula_pattern($str);
    $self;
}

sub atom_map { () }

sub bond_map { () }

sub options {
    my $self = shift;
    $self->{options} = {%{$self->{options}}, @_};    
    $self;
}

sub parse_formula_pattern {
    my ($self, $string) = @_;
    if ($string =~ s/\*//g) {
        $self->{options}{allow_others} = 1;
    }
    my %formula = parse_formula($string);
    $self->{formula_pattern} = \%formula;
}

sub match {
    my ($self, $mol) = @_;
    my $tree = $self->{tree};
    my @ret;
    my $first = 1;

    if ($self->{matched} and $self->{matched} eq $mol) {
        $self->{matched} = 0;
        return 0;
    }

    my $want_formula = $self->{formula_pattern};
    my $have_formula = $mol->formula_hash;
    my $good = 0;
    my @want_keys = keys %$want_formula;
    for my $symbol (@want_keys) {
        no warnings 'uninitialized';
        $have_formula->{$symbol} ||= 0;
        if ($have_formula->{$symbol} >= $want_formula->{$symbol}{low}
            and $have_formula->{$symbol} <= $want_formula->{$symbol}{high}) {
            $good++;
        } else {
            last;
        }
    }
    my @have_keys = keys %$have_formula;
    if ($good == @want_keys 
        and ($good == @have_keys or $self->{options}{allow_others})) {
        $self->{matched} = $mol;
        return 1;
    }
    0;
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
    my ($low, $high, $is_range) = (1,1);
    if ($formula =~ s/^(\d+)(-?)(\d*)//) {
        ($low, $is_range, $high)  = ($1, $2, $3);
        if ($is_range) {
            $high = 1E10 unless length $high;
        } else {
            $high = $low;
        }
    }

    # recursively process rest of formula
    return internal_formula_parser($formula, $low, $high, %elements);
}

sub internal_formula_parser {
    my ($formula, $low, $high, %form) = @_;
    my ($tmp_low, $tmp_high, $tmp_is_range);

    my ($extract, $remainder, $prefix) =
      extract_bracketed($formula, '()', '[^(]*');

    if (defined($extract) and $extract ne '') {
        $extract =~ s/^\((.*)\)$/$1/;
        if ($remainder =~ s/^(\d+)(-?)(\d*)//) {
            ($tmp_low, $tmp_is_range, $tmp_high)  = ($1, $2, $3);
            if ($tmp_is_range) {
                $tmp_high = 1E10 unless length $tmp_high;
            } else {
                $tmp_high = $tmp_low;
            }
            $tmp_low  = $tmp_low * $low;
            $tmp_high = $tmp_high * $high;
        } else {
            $tmp_low  = $low;
            $tmp_high = $high;
        }

        # get formula of prefix ( it has no parens)
        %form = add_formula_strings($prefix, $low, $high, %form) if ($prefix ne '');

        # check remainder for more parens
        %form = internal_formula_parser($remainder, $low, $high, %form)
          if ($remainder ne '');

        # check extract for more parens
        %form =
          internal_formula_parser($extract, $tmp_low, $tmp_high, %form);    
          ## we already know this is ne ''
    } else { # get formula of complete string
        %form = add_formula_strings($remainder, $low, $high, %form)
          if ($remainder ne '');
    }
    return %form;
}

sub add_formula_strings {
    my ($formula, $low_coef, $high_coef, %elements) = @_;

#  print "Getting Formula of $formula\n";
    $formula =~ /^(?:([A-Z][a-z]*)(\d+-?\d*)?)+$/o
        or croak "Invalid Portion of Formula $formula";
    while ($formula =~ m/([A-Z][a-z]*)(?:(\d+)(-?)(\d*))?/go) {
        my ($elm, $low, $is_range, $high) = ($1, $2, $3, $4);
        $low = 1 unless defined $low;
        if ($is_range) {
            $high = 1E10 unless length $high;
        } else {
            $high = $low;
        }
        $elements{$elm}{low} += $low_coef * $low;
        $elements{$elm}{high} += $high_coef * $high;
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

0.10

=head1 SEE ALSO

L<Chemistry::Pattern>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut


