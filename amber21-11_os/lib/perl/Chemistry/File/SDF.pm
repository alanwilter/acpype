package Chemistry::File::SDF;
$VERSION = '0.20';
# $Id: SDF.pm,v 1.9 2005/05/20 19:08:45 itubert Exp $

use base "Chemistry::File";
use Chemistry::Mol;
use Chemistry::File::MDLMol;
use strict;
use warnings;

=head1 NAME

Chemistry::File::SDF - MDL Structure Data File reader/writer

=head1 SYNOPSIS

    use Chemistry::File::SDF;

    # Simple interface (all at once)
    # read all the molecules in the file
    my @mols = Chemistry::Mol->read('myfile.sdf');

    # assuming that the file includes a <PKA> data item...
    print $mols[0]->attr("sdf/data")->{PKA}; 

    # write a bunch of molecules to an SDF file
    Chemistry::Mol->write('myfile.sdf', mols => \@mols);

    # or write just one molecule
    $mol->write('myfile.sdf');


    # Low level interface (one at a time)
    # create reader
    my $reader = Chemistry::Mol->file('myfile.sdf');
    $reader->open('<');
    while (my $mol = $reader->read_mol($reader->fh)) {
        # do something with $mol
    }

=cut

Chemistry::Mol->register_format(sdf => __PACKAGE__);

=head1 DESCRIPTION

MDL SDF (V2000) reader.

This module automatically registers the 'sdf' format with Chemistry::Mol.

The parser returns a list of Chemistry::Mol objects.  SDF data can be accessed
by the $mol->attr method. Attribute names are stored as a hash ref at the
"sdf/data" attribute, as shown in the synopsis. When a data item has a single
line in the SDF file, the attribute is stored as a string; when there's more
than one line, they are stored as an array reference. The rest of the
information on the line that holds the field name is ignored.

This module is part of the PerlMol project, L<http://www.perlmol.org>.

=cut

sub slurp_mol {
    my ($self, $fh, %opts) = @_;
    return if $fh->eof;
    my $s;
    while (<$fh>) {
        last if /^\$\$\$\$/;
        $s .= $_;
    }
    $s =~ s/\r\n?/\n/g; # normalize EOL
    $s;
}

sub skip_mol {
    my ($self, $fh, %opts) = @_;
    return if $fh->eof;
    while (<$fh>) {
        return 1 if /^\$\$\$\$/;
    }
    return 0;
}

sub read_mol {
    my ($self, $fh, %opts) = @_;
    my $s = $self->slurp_mol($fh, %opts) or return;
    my $mol = Chemistry::File::MDLMol->parse_string($s, %opts);
    $self->parse_data($mol, $s);
    $mol;
}

sub parse_data {
    my ($self, $mol, $mol_string) = @_;
    my (@items) = split /\n>/, $mol_string; 
    shift @items; # drop everything until first datum
    my %data_block;
    for my $item (@items) {
        my ($header, @data) = split /\n/, $item;
        my ($field_name) = $header =~ /<(.*?)>/g;
        warn "SDF: no field name\n", next unless $field_name;
        #$mol->attr("sdf/$field_name", @data == 1 ? $data[0] : \@data);
        $data_block{$field_name} = @data == 1 ? $data[0] : \@data;
        
    }
    $mol->attr("sdf/data", \%data_block);
}

sub write_string {
    my ($self, $mol_ref, %opts) = @_;
    my @mols;
    my $ret = '';

    if ($opts{mols}) {
        @mols = @{$opts{mols}};
    } else {
        @mols = $mol_ref; 
    }

    for my $mol (@mols) {
        $ret .= $mol->print(format => 'mdl');    
        $ret .= format_data($mol->attr('sdf/data')) . '$$$$'."\n";
    }
    $ret;
}

sub format_data {
    my ($data) = @_;
    my $ret = '';
    return $ret unless $data;
    for my $field_name (sort keys %$data) {
        $ret .= ">  <$field_name>\n";
        my $value = $data->{$field_name};
        if (ref $value) {
            $ret .= join "\n", @$value;
        } else {
            $ret .= "$value\n";
        }
        $ret .= "\n";
    }
    $ret;
}

sub file_is {
    my ($self, $fname) = @_;
    
    return 1 if $fname =~ /\.sdf?$/i;
    return 0;
}

sub name_is {
    my ($self, $fname) = @_;
    $fname =~ /\.sdf?$/i;
}

sub string_is {
    my ($self, $s) = @_;
    /\$\$\$\$/ ? 1 : 0;
}
1;

=head1 CAVEATS

Note that by storing the SDF data as a hash, there can be only one field with
a given name. The SDF format description is not entirely clear in this regard.
Also note that SDF data field names are considered to be case-sensitive.

=head1 VERSION

0.20

=head1 SEE ALSO

L<Chemistry::Mol>

The MDL file format specification.
L<http://www.mdl.com/downloads/public/ctfile/ctfile.pdf> or
Arthur Dalby et al., J. Chem. Inf. Comput. Sci, 1992, 32, 244-255.

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

