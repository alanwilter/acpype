package MOL2;

use base qw(Chemistry::File);
use Chemistry::MacroMol;
use Carp;
use strict;
use warnings;

=head1 NAME

MOL2 - 'mol2' file format reader

=head1 SYNOPSIS

    use MOL2.pm;

    # read a 'mol2' file
    my $macro_mol = Chemistry::Mol->read("myfile.mol2");

=cut

=head1 DESCRIPTION

This module was created in analogy to the PDB file reader provided in the
Chemistry perl module (see http://www.perlmol.org). 
Date of last change: 06/11/2013
The current version of this module only reads 'mol2' files and exclusively 
uses the information provided in the subsection.
    
 @<TRIPOS>ATOM

The module automatically detects 'mol2' files, so that 'mol2' files may be 
identified and read by Chemistry::Mol->read(). For autodetection purposes, 
it assumes that files ending in .mol2 or containing a line 

@<TRIPOS>MOLECULE  

are 'mol2' files.

=head2 Properties

Upon reading files, the module stores information in the following places:

=over

=item $atom->name

The atom name given in the mol2 file, e.g. "HE1".

=item $atom->attr("mol2/serial_number")

The serial number for the atom, that is provided in the 'mol2' file.

=item $atom->attr("mol2/residue_name")

The name of the residue given in the mol2 file.

=back

=cut


Chemistry::Mol->register_format(mol2 => __PACKAGE__);

sub read_mol {
    my ($self, $fh, %options) = @_;
    return if $fh->eof;

    my $mol_class = $options{mol_class} || "Chemistry::MacroMol";
    my $mol = $mol_class->new;

    local $_;
	
	my $start_atom_reading = 0;

    while (<$fh>) {
		if(/^@<TRIPOS>ATOM/){
			$start_atom_reading = 1;
			next;
		}
		elsif(/^@<TRIPOS>BOND/) {
			return $mol;
		} 
		elsif($start_atom_reading == 1){
			my ($atom_n, $atom_name, $x, $y, $z, $type, $res_id, $res_name, $charge);
			chomp($_);
			$_ =~ s/^\s+//g;
            my @line = split(/\s+/, $_);
			$atom_n = $line[0];
			$atom_name = $line[1];
			$x = $line[2];
			$y = $line[3];
			$z = $line[4];
			$type = $line[5];
			$res_id = $line[6];
			$res_name = $line[7];
			$charge = $line[8];
			
            $type =~ s/ //g;
			my @split_type = split(/\./, $type);
			my $symbol = $split_type[0]; 
            $atom_name =~ s/ //g;
	    	my $a = $mol->new_atom(
				symbol => $symbol, 
				coords => [$x, $y, $z], 
                name => $atom_name,
				formal_charge => $charge,
				type => $type,
	    	);
			if(length $res_name > 3){
				$res_name = substr($res_name, 0, 3);
			}
            $a->attr('mol2/residue_name', "$res_name");
            $a->attr('mol2/serial_number', $atom_n*1);
		}
    }
    return $mol;
}


sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.mol2$/i;
}

sub file_is {
    my ($class, $fname) = @_;
    
    return 1 if $fname =~ /\.mol2$/i;

    open F, $fname or croak "Could not open file $fname";
    
    while (<F>){
		if (/^@<TRIPOS>MOLECULE/) {
	    	close F;
	    	return 1;
		}
    }

    return 0;
}

1;


=head4 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>, L<Chemistry::File>,
L<http://www.perlmol.org/>.

The MOL2 format description at 
L<http://www.tripos.com/data/support/mol2.pdf>

=cut

