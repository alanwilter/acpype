package Chemistry::Mok;

$VERSION = '0.25';
# $Id: Mok.pm,v 1.10 2005/05/16 21:54:21 itubert Exp $

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File ':auto';
use Chemistry::Pattern;
use Chemistry::Bond::Find qw(find_bonds assign_bond_orders);
use Chemistry::Ring 'aromatize_mol';
use Chemistry::3DBuilder 'build_3d';
use Text::Balanced ':ALL';
use Scalar::Util 'blessed';
use Data::Dumper;
use Carp;

our $DEBUG = 0;

=head1 NAME

Chemistry::Mok - molecular awk interpreter

=head1 SYNOPSIS

    use Chemistry::Mok;
    $code = '/CS/g{ $n++; $l += $match->bond_map(0)->length }
        END { printf "Average C-S bond length: %.3f\n", $l/$n; }';

    my $mok = Chemistry::Mok->new($code);
    $mok->run({ format => mdlmol }, glob("*.mol"));

=head1 DESCRIPTION

This module is the engine behind the mok program. See mok(1) for a detailed
description of the language. Mok is part of the PerlMol project,
L<http://www.perlmol.org>.

=head1 METHODS

=over

=cut

sub tokenize {
    my ($self, $code) = @_;

    $code =~ s/\s*$//; # Text::Balanced complains about trailing whitespace
    #$code =~ s/^\s*#.*//g; # remove comments at the top of the file
    #unless($code =~ /^\s*([\/{#]|sub|BEGIN|END)/) {
    unless($code =~ /^(\s*#.*)*\s*([\/{]|sub|BEGIN|END|\w+:\s*\/)/) {
        print "MOK: adding implicit braces\n" if $DEBUG;
        $code = "{$code}"; # add implicit brackets for simple one-liners
    }
    #print "code = '$code'\n";
    # (patt opt?)? code | sub code
    my @toks = extract_multiple(my $c = $code,
        [
            { 'Chemistry::Mok::Comment' => 
                qr/\s*#.*\s*/ },
            { 'Chemistry::Mok::Patt' => 
                sub { scalar extract_delimited($_[0],'/') } },
            { 'Chemistry::Mok::Sub'  => 
                qr/\s*(?:END|BEGIN|sub\s+\w+)\s*/ },
            { 'Chemistry::Mok::Block' => 
                sub { scalar extract_codeblock($_[0],'{') } },
            { 'Chemistry::Mok::PattLang' => 
                qr/(\s*\w+):(?=\s*\/)/ },
            { 'Chemistry::Mok::Opts' => 
                qr/[gopGOP]+/ },
        ],
    );
    die "Mok: error extracting: $@" if $@;
    print "MOK: TOKENS:\n", Dumper(\@toks), "\nCODE:<<<<$code>>>>\n\n"
        if $DEBUG;
    @toks;
}

sub parse {
    my ($self, @toks)  = @_;

    my (@subs, @blocks);
    for my $tok (@toks) {
        blessed $tok or die "unparsable token '$tok'\n";
    }

###  new parser

    my $st = 1;
    my ($patt, $opts, $block, $sub, $pattlang) = ('') x 5;
    my ($save) = 0;
    my $line;
    my $next_line = 1;
    while (my $tok = shift @toks) {
        $line = $next_line; 
        $next_line += $$tok =~ y/\n//;
        print "MOK: LINE=$line;\nTOK=<<<<$$tok>>>>;\nNEXT_LINE=$next_line\n\n" 
            if $DEBUG;
        next if $tok->isa("Chemistry::Mok::Comment");
        if ($st == 1) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            } elsif ($tok->isa("Chemistry::Mok::Sub")) {
                $sub = $$tok,       $st = 5,    next;
            } elsif ($tok->isa("Chemistry::Mok::PattLang")) {
                $pattlang = $$tok,  $st = 4,    next;
            } elsif ($tok->isa("Chemistry::Mok::Patt")) {
                $patt = $$tok,      $st = 2,    next;
            }
        } elsif ($st == 2) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            } elsif ($tok->isa("Chemistry::Mok::Opts")){
                $opts = $$tok,      $st = 3,    next;
            }
        } elsif ($st == 3) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            }
        } elsif ($st == 4) {
            if ($tok->isa("Chemistry::Mok::Patt")){
                $patt = $$tok,      $st = 2,    next;
            }
        } elsif ($st == 5) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            }
        } else {
            confess "unknown state '$st'";
        }
        if ($save) { # save block and go back to state 1
            if ($sub) {
                push @subs, { block => "$sub $$tok", line => $line };
            } else {
                push @blocks, { patt => $patt, opts => $opts, 
                    pattlang => $pattlang, block => $$tok,
                    line => $line};
            }
            $patt = $opts = $pattlang = $block = $sub = '';
            $st = 1,    $save = 0,  next;
        } else {
            die "unexpected token '$$tok' (type '" . ref($tok) . "'\n";
        }
    }
    print "MOK: BLOCKS\n", Dumper(\@blocks), "\nSUBS:\n", Dumper(\@subs), "\n"
        if $DEBUG;

    \@subs, \@blocks;
}

sub compile_subs {
    my ($self, @subs) = @_;
    my $pack = $self->{package};

    for my $sub (@subs) {
        my $code = <<END;
            package Chemistry::Mok::UserCode::$pack;
            no strict;
            no warnings;
#line $sub->{line} "mok code"
            $sub->{block}
END
        print "MOK: COMPILING SUB: <<<<$code>>>>\n\n" if $DEBUG;
        eval $code;
        die "Mok: error compiling sub: $@" if $@;
    }
}

sub compile_blocks {
    my ($self, @blocks) = @_;
    my $pack = $self->{package};
    my $format = $self->{pattern_format};
    my @compiled_blocks;

    for my $block (@blocks) {
        #use Data::Dumper; print Dumper $block;
        my $code = <<END;
            package Chemistry::Mok::UserCode::$pack;
            no strict;
            no warnings;
            sub {
                my (\$mol, \$file, \$match, \$patt) = \@_;
                my (\$MOL, \$FILE, \$MATCH, \$PATT, \$FH) = \@_;
                my (\@A) = \$MATCH ? \$MATCH->atom_map : \$MOL->atoms;
                my (\@B) = \$MATCH ? \$MATCH->bond_map : \$MOL->bonds;
#line $block->{line} "mok code"
                $block->{block};
            }
END
        print "MOK: COMPILING BLOCK: <<<<$code>>>>\n\n" if $DEBUG;
        my $sub = eval $code;
        die "Mol: Error compiling block: $@" if $@;

        my ($patt, $patt_str);
        if ($block->{patt}) {
            $block->{patt} =~ m#^/(.*)/$#;
            $patt_str = $1;
            $patt = Chemistry::Pattern->parse($patt_str, 
                format => $block->{pattlang} || $format);
            $patt->attr(global => 1) if $block->{opts} =~ /g/;
            $patt->options(overlap => 0) if $block->{opts} =~ /O/;
            $patt->options(permute => 1) if $block->{opts} =~ /p/;
        } 
        push @compiled_blocks, {'sub' => $sub, 
            patt => $patt, patt_str => $patt_str};
    }
    \@compiled_blocks;
}

=item Chemistry::Mok->new($code, %options)

Compile the code and return a Chemistry::Mok object. Available options:

=over

=item C<package>

If the C<package> option is given, the code runs in the
Chemistry::Mok::UserCode::$options{package} package instead of the
Chemistry::Mok::UserCode::Default package. Specifying a package name is
recommended if you have more than one mok object and you are using global
varaibles, in order to avoid namespace clashes.

=item C<pattern_format>

The name of the format which will be used for parsing slash-delimited patterns
that don't define an explicit format. Mok versions until 0.16 only used the
'smiles' format, but newer versions can use other formats such as 'smarts',
'midas', 'formula_pattern', and 'sln', if available. The default is 'smarts'.

=back

=cut

sub new {
    my ($class, $code, @a) = @_;
    my %opts;

    # for backwards compatibility with Chemistry::Mok->new($code, $package)
    unshift @a, "package" if (@a == 1);
    %opts = @a;
        
    my $self = bless {
        'package'      => $opts{package} || "Default",
        pattern_format => $opts{pattern_format} || "smarts",
    }, $class;

    $self->setup_package;
    my @toks = $self->tokenize($code);
    my ($subs, $blocks) = $self->parse(@toks);
    $self->compile_subs(@$subs);
    $self->{blocks} = $self->compile_blocks(@$blocks);
    
    return $self;
}

sub setup_package {
    my ($self) = @_;
    my $usr_pack = $self->{package};
    # import convenience functions into the user's namespace
    eval <<EVAL;
          package Chemistry::Mok::UserCode::$usr_pack;
          use Chemistry::Atom ':all';
          use Chemistry::Ring ':all';
          use Chemistry::Ring::Find ':all';
          use Chemistry::Bond::Find ':all';
          use Chemistry::Canonicalize ':all';
          use Chemistry::InternalCoords::Builder ':all';
          use Chemistry::Isotope ':all';
          use Math::VectorReal ':all';
          use Chemistry::3DBuilder ':all';
          sub println { print "\@_", "\n" }
EVAL
    die "Mok: error setting up 'Chemistry::Mok::UserCode::$usr_pack' $@" if $@;
}

=item $mok->run($options, @args)

Run the code on the filenames contained in @args. $options is a hash reference
with runtime options. Available options:

=over

=item build_3d

Generate 3D coordinates using Chemistry::3DBuilder.

=item aromatize          

"Aromatize" each molecule as it is read. This is needed for example for
matching SMARTS patterns that use aromaticity or ring primitives.

=item delete_dummies

Delete dummy atoms after reading each molecule. A dummy atom is defined as an
atom with an unknown symbol (i.e., it doesn't appear on the periodic table), or
an atomic number of zero.

=item find_bonds

If set to a true value, find bonds. Use it when reading files with no bond
information but 3D coordinates to detect the bonds if needed (for example, if
you want to do match a pattern that includes bonds). If the file has explicit
bonds, mok will not try to find the bonds, but it will reassign the bond orders
from scratch.

=item format

The format used when calling $mol_class->read. If not given, $mol_class->read
tries to identify the format automatically.

=item mol_class

The molecule class used for reading the files. Defaults to Chemistry::Mol.

=back

=cut

sub run {
    my ($self, $opt, @args) = @_;
    # MAIN LOOP
    my $mol_class = $opt->{mol_class} || "Chemistry::Mol";
    FILE: for my $file (@args) {
        #my (@mols) = $mol_class->read(
        my %reader_opts = (
            format      => $opt->{format},
            mol_class   => $opt->{mol_class},
        );
        my $reader = $mol_class->file(
            $file, 
            %reader_opts,
        );
        $reader->open('<');
        $reader->read_header;
        while (my @mols  = $reader->read_mol($reader->fh, %reader_opts)) {
            MOL: for my $mol (@mols) {
                if ($opt->{delete_dummies}) {
                    $_->delete for grep { ! $_->Z } $mol->atoms;
                }
                if ($opt->{find_bonds}) {
                    find_bonds($mol) unless $mol->bonds;
                    assign_bond_orders($mol);
                }
                if ($opt->{aromatize}) {
                    aromatize_mol($mol);
                }
                if ($opt->{build_3d}) {
                    build_3d($mol);
                }
                BLOCK: for my $block (@{$self->{blocks}}) {
                    my ($code_block, $patt, $patt_str) = 
                        @{$block}{qw(sub patt patt_str)};
                    if ($patt) {
                        MATCH: while ($patt->match($mol)) {
                            $code_block->($mol, $file, $patt, 
                                $patt_str, $reader->fh);
                            last unless $patt->attr('global');
                        }
                    } else {
                        $code_block->($mol, $file, $patt, 
                            $patt_str, $reader->fh);
                    }
                }
            }
        }
    }
}

1;

__END__

=back

=head1 VERSION

0.25

=head1 SEE ALSO

L<mok>, L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

