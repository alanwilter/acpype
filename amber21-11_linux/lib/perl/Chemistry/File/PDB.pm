package Chemistry::File::PDB;

$VERSION = '0.21';
# $Id: PDB.pm,v 1.11 2005/05/16 22:33:48 itubert Exp $

use base qw(Chemistry::File);
use Chemistry::MacroMol;
use Chemistry::Domain;
use Scalar::Util qw(weaken);
use Carp;
use strict;
use warnings;

=head1 NAME

Chemistry::File::PDB - Protein Data Bank file format reader/writer

=head1 SYNOPSIS

    use Chemistry::File::PDB;

    # read a PDB file
    my $macro_mol = Chemistry::MacroMol->read("myfile.pdb");

    # write a PDB file
    $macro_mol->write("out.pdb");

    # read all models in a multi-model file
    my @mols = Chemistry::MacroMol->read("models.pdb");

    # read one model at a time
    my $file = Chemistry::MacroMol->file("models.pdb");
    $file->open;
    while (my $mol = $file->read_mol($file->fh)) {
        # do something with $mol
    }

=cut

=head1 DESCRIPTION

This module reads and writes PDB files. The PDB file format is commonly used to
describe proteins, particularly those stored in the Protein Data Bank
(L<http://www.rcsb.org/pdb/>). The current version of this module only uses the
following record types, ignoring everything else:

    ATOM
    HETATM
    ENDMDL 

This module automatically registers the 'pdb' format with Chemistry::Mol,
so that PDB files may be identified and read by Chemistry::Mol->read(). For 
autodetection purpuses, it assumes that files ending in .pdb or having 
a line matching /^(ATOM  |HETATM)/ are PDB files.

The PDB reader and writer is designed for dealing with Chemistry::MacroMol
objects, but it can also create and use Chemistry::Mol objects by throwing some
information away.

=head2 Properties

When reading and writing files, this module stores or gets some of the
information in the following places:

=over

=item $domain->type

The residue type, such as "ARG".

=item $domain->name

The type and sequence number, such as "ARG114".

=item $domain->attr("pdb/sequence_number")

The residue sequence number as given in the PDB file.

=item $domain->attr("pdb/chain_id")

The chain to which this residue belongs (one character).

=item $domain->attr("pdb/insertion_code")

The residue insertion code (see the PDB specification for details).

=item $atom->name

The PDB atom name, such as "CA".

=item $atom->attr("pdb/residue_name")

The name of the residue, as discussed above.

=item $atom->attr("pdb/serial_number")

The serial number for the atom, as given in the PDB file.

=back

If some of this information is not available when writing a PDB file, this 
module tries to make it up (by counting the atoms or residues, for example).
The default residue name for writing is UNK (unknown). Atom names are just the
atomic symbols.

=head2 Multi-model files

If a PDB file has multiple models (separated by ENDMDL records), each call to
read_mol will return one model.

=cut

Chemistry::Mol->register_format(pdb => __PACKAGE__);

sub read_mol {
    my ($self, $fh, %options) = @_;
    return if $fh->eof;

    my $domain;
    my $mol_class = $options{mol_class} || "Chemistry::MacroMol";
    my $mol = $mol_class->new;
    my $is_macro = $mol->isa('Chemistry::MacroMol');
    $domain = $mol unless $is_macro;

    local $_;
    my $curr_residue = 0;
    while (<$fh>) {
	if (/^TER/) {
            # TODO read separated molecules?
	} elsif (/^(?:HETATM|ATOM)/) {
	    my ($atom_n, $symbol, $suff, $res_name, $ins_code, $chain_id,
                $seq_n, $x, $y, $z
            ) = unpack "x6A5x1A2A2x1A3A1A1A4x4A8A8A8", $_;
	    #print "S:$symbol; N:$name; x:$x; y:$y; z:$z\n";
            $seq_n *= 1;
            if (!$domain || $seq_n != $curr_residue) {
                # new residue
                if ($is_macro) {
                    $domain = Chemistry::Domain->new(
                        parent => $mol,      name => "$res_name$seq_n",
                        type   => $res_name);
                    $mol->add_domain($domain);
                    weaken $domain;
                    $domain->attr('pdb/sequence_number', $seq_n);
                    $domain->attr('pdb/chain_id',        $chain_id);
                    $domain->attr('pdb/insertion_code',  $ins_code);
                }
                $curr_residue = $seq_n;
            }
            my $atom_name = $symbol.$suff;
            $atom_name =~ s/ //g;
            $symbol    =~ s/\d//g;
	    my $a = $domain->new_atom(
		symbol => $symbol, 
		coords => [$x, $y, $z], 
                name => $atom_name,
	    );
            $a->attr('pdb/residue_name', "$res_name$seq_n");
            $a->attr('pdb/serial_number', $atom_n*1);
            $a->attr('pdb/residue', $domain);
	} elsif (/^ENDMDL/) {
            return $mol;
	}
    }
    return $mol;
}


sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.pdb$/i;
}

sub file_is {
    my ($class, $fname) = @_;
    
    return 1 if $fname =~ /\.pdb$/i;

    open F, $fname or croak "Could not open file $fname";
    
    while (<F>){
	if (/^ATOM  / or /^HETATM/) {
	    close F;
	    return 1;
	}
    }

    return 0;
}

sub write_string {
    my ($class, $mol, %opts) = @_;
    my $ret = '';
    if ($mol->isa("Chemistry::MacroMol")) {
        my $i = 1;
        my $j = 1;
        for my $res ($mol->domains) {
            my $seq_n = $res->attr("pdb/sequence_number");
            $seq_n = $j++ unless defined $seq_n;
            my $res_name = substr($res->name, 0, 3) || "UNK";

            for my $atom ($res->atoms) {
                my $serial_n  = $res->attr("pdb/serial_number");
                $serial_n     = $i++ unless defined $serial_n;
                my @coords    = $atom->coords->array;

                # make atom name
                no warnings 'uninitialized';
                my $symbol    = substr $atom->symbol, 0, 2;
                my $atom_name = $atom->name || $symbol;
                $atom_name =~ /^(\dH|$symbol)(.{0,2})/;
                #print "NAME: '$atom_name' ($1,$2); SYMBOL: '$symbol'\n";
                $atom_name = sprintf "%2s%-2s", $1, $2;
                $ret .= sprintf "ATOM  %5d %4s %-3s  %4d    %8.3f%8.3f%8.3f\n",
                    $serial_n, $atom_name, $res_name, $seq_n, @coords;
            }
        }
    } else {
        my $i = 1;
        for my $atom ($mol->atoms) {
            my @coords = $atom->coords->array;
            $ret .= sprintf "HETATM%5d %2s   UNK     1    %8.3f%8.3f%8.3f\n",
                $i++, $atom->symbol, @coords;
        }
    }
    $ret .= "TER\nEND\n";
    $ret;
}

1;

=head1 VERSION

0.21

=head1 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>, L<Chemistry::File>,
L<http://www.perlmol.org/>.

The PDB format description at 
L<http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html>

There is another PDB reader in Perl, as part of the BioPerl project:
L<Bio::Structure::IO::pdb>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

