package Chemistry::File::SMARTS;

$VERSION = "0.22";
# $Id: SMARTS.pm,v 1.8 2005/05/16 22:20:24 itubert Exp $

use 5.006;
use strict;
use warnings;
use Chemistry::Pattern;
use base "Chemistry::File";
use Carp;
#use Data::Dumper;
use List::Util 'sum', 'first';

=head1 NAME

Chemistry::File::SMARTS - SMARTS chemical substructure pattern linear notation parser

=head1 SYNOPSYS

    #!/usr/bin/perl
    use Chemistry::File::SMARTS;

    # this string matches an oxygen next to an atom with three 
    # neighbors, one of which is a hydrogen, and a positive charge
    my $smarts = 'O[D3H+]'; 

    # parse a SMARTS string and compile it into a
    # Chemistry::Pattern object
    my $patt = Chemistry::Pattern->parse("$smarts", format => 'smarts');

    # find matches of the pattern in a Chemistry::Mol object $mol
    my $mol = Chemistry::Mol->read("myfile.mol");
    while ($patt->match($mol)) {
        print "pattern matches atoms: ", $patt->atom_map, "\n"
    }

    # NOTE: if the SMARTS pattern relies on aromaticity or ring
    # properties, you have to make sure that the target 
    # molecule is "aromatized" first:
    my $smarts = 'c:a';
    my $patt = Chemistry::Pattern->parse("$smarts", format => 'smarts');
    use Chemistry::Ring 'aromatize_mol';
    aromatize_mol($mol);  # <--- AROMATIZE!!!
    while ($patt->match($mol)) {
        print "pattern matches atoms: ", $patt->atom_map, "\n"
    }

    # Note that "atom mapping numbers" end up as $atom->name
    my $patt = Chemistry::Pattern->parse("[C:7][C:8]", format => 'smarts');
    print $patt->atoms(1)->name;    # prints 7


=head1 DESCRIPTION

This module parse a SMARTS (SMiles ARbitrary Target Specification) string,
generating a L<Chemistry::Pattern> object.  It is a file I/O driver for the
PerlMol toolkit; it's not called directly but by means of the
Chemistry::Pattern->parse class method.

For a detailed description of the SMARTS language, see
L<http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>. Note that
this module doesn't implement the full language, as detailed under CAVEATS.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=cut

# Initialization
Chemistry::Mol->register_format(smarts => __PACKAGE__);
our $DEBUG = 0;

# Chemistry::File interface

sub parse_string {
    my ($self, $s, %opts) = @_;
    my $patt = parse_smarts($s, \%opts);
    $patt->name($s);
    $patt;
}


############## SMARTS PARSER ##############

sub parse_smarts {
    my ($s, $opts) = @_;
    my %digits;

    my $toks = tokenize($s);
    my $tok = shift(@$toks);
    print "tok: $tok\n" if $DEBUG;
    my $mol_class = $opts->{mol_class} || "Chemistry::Pattern";
    my $patt = $mol_class->new();
    $patt->options($opts->{pattern_options}||{});
    my @atom_stack;
    my $current_atom = parse_atom($patt, $tok, $toks);
    while (defined ($tok = shift @$toks)) {
        print "tok: $tok\n" if $DEBUG;
        if ($tok eq '(') {
            push @atom_stack, $current_atom;
        } elsif ($tok eq ')') {
            $current_atom = pop @atom_stack;
        } else {  # bond, atom
            my $next_tok = shift @$toks;
            if ($next_tok =~ /^\d+$/) {  # digit
                if ($digits{$next_tok}) {  # close ring
                    parse_bond($patt, $tok, $current_atom, 
                        $digits{$next_tok});
                    $digits{$next_tok} = undef;
                } else { # open ring
                    $digits{$next_tok} = $current_atom;
                }
            } else {
                my $next_atom = parse_atom($patt, $next_tok, $toks);
                parse_bond($patt, $tok, $current_atom, $next_atom);
                $current_atom = $next_atom;
            }
        }
    }
    $patt;
}

sub parse_atom {
    my ($patt, $s, $toks) = @_;
    
    my $n_rec = 0;
    my $name;

    if ($s =~ s/:(\d+)$//) {
        $name = $1;
    }
    my $expr = 
        join " and ", map { 
            join " || ", map { 
                join ' && ', map {
                    parse_atomic_primitive($_, \$n_rec);
                } split '&', $_;
            } split ',', $_;
        } split ';', $s;

    
    my @recs;
    for (1 .. $n_rec) {
        my $rec_smarts = shift @$toks;
        my $rec = Chemistry::Pattern->parse($rec_smarts, 
            pattern_options => {overlap=>0, permute=>0},
            format => 'smarts',
        );
        push @recs, $rec;
    }

    print "atom expr: $expr\n" if $DEBUG;
    my $sub = eval <<SUB;
        sub {
            no warnings;
            my (\$patt, \$atom) = \@_;
            $expr;
        };
SUB
    my $atom = $patt->atom_class->new(
        test_sub => $sub,
        name => $name,
    );
    $patt->add_atom($atom);
    $atom;
}


# missing primitives: @, @@
sub parse_atomic_primitive {
    my ($s, $n_rec) = @_;
    local $_ = $s;
    my @terms;
    no warnings 'uninitialized';

    s/^(!?)H// &&           # Hydrogen atom
        push @terms, "$1(\$atom->symbol eq 'H')";

    s/(!?)
        (Zr|Zn|Yb|Y|Xe|W|V|U|Tm|Tl|Ti|Th|
        Te|Tc|Tb|Ta|Sr|Sn|Sm|Si|Sg|Se|Sc|Sb|S|Ru|Rn|Rh|Rf|Re|Rb|Ra|
        Pu|Pt|Pr|Po|Pm|Pd|Pb|Pa|P|Os|O|Np|No|Ni|Ne|Nd|Nb|Na|N|Mt|Mt|
        Mo|Mn|Mg|Md|Lu|Lr|Li|La|Kr|K|Ir|In|I|Hs|Hs|Ho|Hg|Hf|He|Ge|
        Gd|Ga|Fr|Fm|Fe|F|Eu|Es|Er|Dy|Ds|Db|Cu|Cs|Cr|Co|Cm|Cl|Cf|Ce|
        Cd|Ca|C|Br|Bk|Bi|Bh|Be|Ba|B|Au|At|As|Ar|Am|Al|Ag|Ac)
    //x # Order is reverse alphabetical to ensure longest match
        && push @terms, "$1(\$atom->symbol eq '$2' && ! \$atom->aromatic)";

    s/(!?)\*// &&           # wildcard
        push @terms, "${1}1";

    s/(!?)D(\d?)// &&       # explicit connections (shouldn't count h)
        push @terms, "$1(\$atom->bonds == " . (length $2 ? $2 : 1) . ')';

    s/(!?)a// &&            # aromatic
        push @terms, "$1\$atom->aromatic";

    s/(!?)A// &&            # aliphatic
        push @terms, "$1(!\$atom->aromatic)";

    s/(!?)X(\d?)// &&       # total connections (should add implicit H)
        push @terms, "$1(\$atom->bonds + \$atom->hydrogens == " 
            . (length $2 ? $2 : 1) . ')';

    s/(!?)v(\d?)// &&       # valence
        push @terms, "$1(\$atom->valence == " 
            . (length $2 ? $2 : 1) . ')';

    s/(!?)H(\d?)// &&    # total H-count
        push @terms, "$1(sum(map {\$_->symbol eq 'H'} \$atom->neighbors) +
            \$atom->hydrogens == " . (length $2 ? $2 : 1) . ')';

    s/(!?)h(\d?)// &&    # implicit H-count
        push @terms, "$1(\$atom->hydrogens == " . (length $2 ? $2 : 1) . ')';

    s/(!?)R(\d?)// &&    # number of rings
        push @terms, "$1(\@{\$atom->attr('ring/rings')||[]} " 
            . (length $2 ? "== $2" : "") . ')';
    
    s/(!?)r(\d?)// &&    # ring size
        push @terms, "$1(first { " . (length $2 ? "\$_->atoms == $2" : "1") 
            ." } \@{\$atom->attr('ring/rings')||[]} )";
    
    s/(!?)#(\d+)// &&       # atomic number
        push @terms, "$1(\$atom->Z == $2)";

    s/(!?)([+-]\d+)// &&    # numerical charge 
        push @terms, "$1(\$atom->formal_charge == $2)";

    s/(!?)(\++)// &&        # positive charge
        push @terms, "$1(\$atom->formal_charge == " . length($2) . ')';

    s/(!?)(-+)// &&         # negative charge 
        push @terms, "$1(\$atom->formal_charge == -" . length($2) . ')';

    s/(!?)(\d+)// &&        # mass
        push @terms, "$1(\$atom->mass == $2)";

    s/(!?)([cnosp])// &&    # aromatic symbol
        push @terms, "$1(\$atom->symbol eq '@{[uc $2]}' && \$atom->aromatic)";

    #s/(!?)([A-Z][a-z]?)// &&    # aliphatic symbol
        #push @terms, "$1(\$atom->symbol eq '$2' && ! \$atom->aromatic)";

    while (s/(!?)\$//) {    #  recursive SMARTS
        push @terms, qq{$1(\$recs[$$n_rec]->match(\$atom->parent,atom=>\$atom))};
        $$n_rec++;
    }

    join ' && ', @terms;
}

sub parse_bond {
    my ($patt, $s, @atoms) = @_;

    return if $s eq '.'; # the disconnected non-bond

    my $expr;

    if ($s) {
        $expr = 
            join " and ", map { 
                join " || ", map { 
                    join ' && ', map {
                        parse_bond_primitive($_);
                    } split '&', $_;
                } split ',', $_;
            } split ';', $s;
    } else {
        $expr = '($bond->order == 1 || $bond->aromatic)';
    }

    print "bond expr: $expr\n" if $DEBUG;
    my $sub = eval <<SUB;
        sub {
            no warnings;
            my (\$patt, \$bond) = \@_;
            $expr;
        };
SUB
    my $bond = $patt->bond_class->new(test_sub => $sub, 
        atoms => \@atoms);
    $patt->add_bond($bond);
    $bond;
}

sub parse_bond_primitive {
    local ($_) = @_;
    my @terms;

    s/(!?)~// &&        # wildcard
        push @terms, "${1}1";

    s/(!?)-// &&        # single
        push @terms, "$1(\$bond->order == 1 && !\$bond->aromatic)";

    s/(!?)=// &&        # double
        push @terms, "$1(\$bond->order == 2)";

    s/(!?)#// &&        # triple
        push @terms, "$1(\$bond->order == 3)";

    s/(!?):// &&        # triple
        push @terms, "$1\$bond->aromatic";

    s/(!?)\@// &&       # ring bond
        push @terms, "$1(\@{\$bond->attr('ring/rings')||[]} )";
    
    join ' && ', @terms;
}

my %ORGANIC_ELEMS = (
    Br => 1, Cl => 1, B => 1, C => 1, N => 1, O => 1, P => 1, S => 1, 
    F => 1, I => 1, s => 1, p => 1, o => 1, n => 1, c => 1, b => 1,
);

sub tokenize {
    my ($s) = @_;
    my $state = 3;
    my $paren_depth = 0;
    my $atom;
    my $digit;
    my $rec_smart;
    my @rec_smarts;
    my $bond = '';
    my $symbol;
    my @toks;
    
    # 0 expects atom or branch
    # after atom: atom or bond or branch
    # after bond: atom or branch or /,;/
    # token types: atom, bond, (, )
    my @chars = split '', $s;
    my $char;
    while (defined ($char = shift @chars)) {
        print "char: $char\n" if $DEBUG;
        if ($state == 0) { # expect atom or branch (not used!)
            push(@toks,  $char), $state = 2, next if ($char =~ /[()]/);
            $state = 3, redo;
        } elsif ($state == 1) { # in complex atom
            if ($char eq ']') {
                push @toks, $atom, @rec_smarts;
                @rec_smarts = ();
                $state = 4, next;
            }
            $atom .= $char;
            $state = 5 if ($char eq '$'); # recursive smarts
            next;
        } elsif ($state == 2) { # expect bond
            if ($char =~ /[-=:~\\\/#@?,;&.!]/) {
                $bond .= $char;
                next;
            } else {
                push @toks, $bond;
                $bond = '';
                $state = 3, redo;
            }
        } elsif ($state == 3) {  # expect atom
            $state = 1, $atom = '', next if $char eq '[';
            if ($char eq '%') {
                $state = 7, next;
            }
            $symbol = $char;
            push(@toks, $symbol), last unless @chars;
            $char = shift @chars;
            if ($ORGANIC_ELEMS{$symbol.$char}) {
                push @toks, $symbol.$char;
                $state = 4, next;
            } else {
                push @toks, $symbol;
                $state = 4, redo;
            }
        } elsif ($state == 4) {  # expect atom or bond or branch
            push(@toks, $char), next if ($char =~ /[()]/); # branch
            $state = 2, redo; # bond
        } elsif ($state == 5) {  # expect left paren
            croak "expected (" unless $char eq '(';
            $rec_smart = '';
            $paren_depth++, $state = 6, next;
        } elsif ($state == 6) {  # within recursive smarts
            #croak "Recursive SMARTS not implemented yet\n";
            $paren_depth++ if $char eq '(';
            $paren_depth-- if $char eq ')';
            unless ($paren_depth) {
                push @rec_smarts, $rec_smart;
                $state = 1, next;
            }
            $rec_smart .= $char;
        } elsif ($state == 7) {  # double digit
            $digit = $char . (shift @chars || die "expected second digit");
            push @toks, $digit;
            $state = 2;
        } else {
            die "shouldn't be here";
        }
    }
    #print Dumper \@toks if $DEBUG;
    \@toks;
}


1;

=head1 CAVEATS

The following features are not implemented yet:

=over

=item chirality: @, @@

=item component-level gruouping

That is, the difference between these three cases:

    (SMARTS)
    (SMARTS).(SMARTS)
    (SMARTS).SMARTS

=back

The so-called parser is very lenient, so if you give it something that's not
quite reasonable it will ignore it or interpret it in a strange way without
warning.

As shown in the synopsis, you have to make sure that the molecule is
"aromatized" if you want to apply to it a pattern that relies on aromaticity
or ring properties.

=head1 VERSION

0.22

=head1 SEE ALSO

L<Chemistry::Pattern>, L<Chemistry::Mol>, L<Chemistry::File>,
L<Chemistry::File::SMILES>.

For more information about SMARTS, see the SMARTS Theory Manual at
L<http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

