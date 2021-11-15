package Chemistry::Mol;
$VERSION = '0.35';
# $Id: Mol.pm,v 1.48 2005/05/20 19:01:04 itubert Exp $

=head1 NAME

Chemistry::Mol - Molecule object toolkit

=head1 SYNOPSIS

    use Chemistry::Mol;

    $mol = Chemistry::Mol->new(id => "mol_id", name => "my molecule");
    $c = $mol->new_atom(symbol => "C", coords => [0,0,0]); 
    $o = $mol->new_atom(symbol => "O", coords => [0,0,1.23]); 
    $mol->new_bond(atoms => [$c, $o], order => 3);

    print $mol->print;

=head1 DESCRIPTION

This package, along with Chemistry::Atom and Chemistry::Bond, includes basic
objects and methods to describe molecules. 

The core methods try not to enforce a particular convention.  This means that
only a minimal set of attributes is provided by default, and some attributes
have very loosely defined meaning. This is because each program and file type
has different idea of what each concept (such as bond and atom type) means.
Bonds are defined as a list of atoms (typically two) with an arbitrary type.
Atoms are defined by a symbol and a Z, and may have 3D and internal coordinates
(2D coming soon).

=cut

use 5.006;
use strict;
use warnings;
use Chemistry::Atom;
use Chemistry::Bond;
use Carp;
use base qw(Chemistry::Obj Exporter);
use Storable 'dclone';

our @EXPORT_OK = qw(read_mol);
our @EXPORT = ();
our %EXPORT_TAGS = (
      all  => [@EXPORT, @EXPORT_OK],
);



my %FILE_FORMATS = ();

=head1 METHODS

See also L<Chemistry::Obj> for generic attributes.

=over 4

=item Chemistry::Mol->new(name => value, ...)

Create a new Mol object with the specified attributes. 

    $mol = Chemistry::Mol->new(id => 'm123', name => 'my mol')

is the same as

    Chemistry::Mol->new()
    $mol->id('m123')
    $mol->name('my mol')

=cut

sub new {
    my $class = shift;
    my %args = @_;
    my $self = bless {
        id => $class->nextID,
        byId => {}, 
        atoms => [], 
        bonds => [], 
        name => "",
    }, ref $class || $class;
    $self->$_($args{$_}) for (keys %args);
    return $self;
}

my $N = 0; # atom ID counter
sub nextID { "mol".++$N; }
sub reset_id { $N = 0; }

=item $mol->add_atom($atom, ...)

Add one or more Atom objects to the molecule. Returns the last atom added.

=cut

sub add_atom {
    my $self = shift;
    for my $atom (@_){
        #if ($self->by_id($atom->id)) {
            #croak "Duplicate ID when adding atom '$atom' to mol '$self'";
        #}
        push @{$self->{atoms}}, $atom;
        $self->{byId}{$atom->id} = $atom;
        $atom->parent($self);
    }
    $_[-1];
}

sub add_atom_np {
    my $self = shift;
    for my $atom (@_){
        push @{$self->{atoms}}, $atom;
        $self->{byId}{$atom->id} = $atom;
    }
    $_[-1];
}

=item $mol->atom_class

Returns the atom class that a molecule or molecule class expects to use by
default. L<Chemistry::Mol> objects return "Chemistry::Atom", but subclasses
will likely override this method.

=cut

sub atom_class {
    "Chemistry::Atom";
}

=item $mol->new_atom(name => value, ...)

Shorthand for C<< $mol->add_atom($mol->atom_class->new(name => value, ...)) >>.

=cut

sub new_atom {
    my $self = shift;
    $self->add_atom($self->atom_class->new(@_));
}

=item $mol->delete_atom($atom, ...)

Deletes an atom from the molecule. It automatically deletes all the bonds in
which the atom participates as well. $atom should be a Chemistry::Atom
reference. This method also accepts the atom index, but this use is deprecated
(and buggy if multiple indices are given, unless they are in descending order).

=cut

sub delete_atom {
    my $self = shift;
    for my $i (@_) {
        my ($atom);
        if (ref $i) {
            $atom = $i;
        } else {
            $atom = $self->atoms($i)
                or croak "$self->delete_atom: no such atom $i\n";
        }
        $atom->delete($i);
    }
}

# takes an atom ref to delete and optionally the atom index
# 1) deletes bonds that belonged to atom
# 2) deletes atom
sub _delete_atom {
    my ($self, $atom) = @_;
    my $index = $self->get_atom_index($atom)    
        or croak "$self->delete_atom: no such atom $atom\n";
    my $id = $atom->id;
    $self->delete_bond($atom->bonds);
    delete $self->{byId}{$id};
    splice @{$self->{atoms}}, $index - 1, 1;
}

=item $mol->add_bond($bond, ...)

Add one or more Bond objects to the molecule. Returns the last bond added.

=cut

sub add_bond {
    my $self = shift;
    for my $bond (@_){
        #if ($self->by_id($bond->id)) {
            #croak "Duplicate ID when adding bond '$bond' to mol '$self'";
        #}
        push @{$self->{bonds}}, $bond;
        $self->{byId}{$bond->id} = $bond;
        if ($bond->{deleted}) {
            $_->add_bond($bond) for $bond->atoms;
            $bond->{deleted} = 0;
        }
        $bond->parent($self);
    }
    $_[-1];
}

sub add_bond_np {
    my $self = shift;
    for my $bond (@_){
        push @{$self->{bonds}}, $bond;
        $self->{byId}{$bond->id} = $bond;
    }
    $_[-1];
}

=item $mol->bond_class

Returns the bond class that a molecule or molecule class expects to use by
default. L<Chemistry::Mol> objects return "Chemistry::Bond", but subclasses
will likely override this method.

=cut

sub bond_class {
    "Chemistry::Bond";
}

=item $mol->new_bond(name => value, ...)

Shorthand for C<< $mol->add_bond($mol->bond_class->new(name => value, ...)) >>.

=cut

sub new_bond {
    my $self = shift;
    $self->add_bond($self->bond_class->new(@_));
}

sub get_bond_index {
    my ($self, $bond) = @_;
    my $i;
    for ($self->bonds) {
        ++$i;
        return $i if ($_ eq $bond);
    }
    undef;
}

sub get_atom_index {
    my ($self, $atom) = @_;
    my $i;
    for ($self->atoms) {
        ++$i;
        return $i if ($_ eq $atom);
    }
    undef;
}

=item $mol->delete_bond($bond, ...)

Deletes a bond from the molecule. $bond should be a L<Chemistry::Bond> object.

=cut

# mol deletes bond
# bond tells atoms involved to forget about it

sub delete_bond {
    my $self = shift;
    for my $i (@_){
        my ($bond);
        if (ref $i) {
            $bond = $i;
        } else {
            $bond = $self->bonds($i)
                or croak "$self->delete_bond: no such bond $i\n";
        }
        $bond->delete;
    }
}

sub _delete_bond {
    my ($self, $bond) = @_;
    my $index = $self->get_bond_index($bond)    
        or croak "$self->delete_bond: no such bond $bond\n";
    my $id = $bond->id;
    delete $self->{byId}{$id};
    splice @{$self->{bonds}}, $index - 1, 1;
    $bond->delete_atoms;
}

=item $mol->by_id($id)

Return the atom or bond object with the corresponding id.

=cut

sub by_id {
    my $self = shift;
    my ($id) = @_;
    $self->{byId}{$id};
}

sub _change_id {
    my ($self, $old_id, $new_id) = @_;
    my $ref = $self->{byId}{$old_id};
    $self->{byId}{$new_id} = $ref;
    delete $self->{byId}{$old_id};
}

=item $mol->atoms($n1, ...)

Returns the atoms with the given indices, or all by default. 
Indices start from one, not from zero.

=cut

sub atoms {
    my $self = shift;
    if (@_) {
        my @ats = map {$_ - 1} @_;
        @{$self->{atoms}}[@ats];
    } else {
        @{$self->{atoms}};
    }
}

=item $mol->atoms_by_name($name)

Returns the atoms with the given name (treated as an anchored regular
expression).

=cut

sub atoms_by_name {
    my $self = shift;
    my $re = qr/^$_[0]$/;
    no warnings;
    my @ret = grep {$_->name =~ $re} $self->atoms;
    wantarray ? @ret : $ret[0];
}

=item $mol->sort_atoms($sub_ref)

Sort the atoms in the molecule by using the comparison function given in
$sub_ref. This function should take two atoms as parameters and return -1, 0,
or 1 depending on whether the first atom should go before, same, or after the
second atom. For example, to sort by atomic number, you could use the
following:

    $mol->sort_atoms( sub { $_[0]->Z <=> $_[1]->Z } );

Note that the atoms are passed as parameters and not as the package variables
$a and $b like the core sort function does. This is because $mol->sort will
likely be called from another package and we don't want to play with another
package's symbol table.

=cut

sub sort_atoms {
    my ($self, $sub) = @_;
    my @a = $self->atoms;
    @a = sort { $sub->($a,$b) } @a;
    $self->{atoms} = \@a;
    $self;
}

=item $mol->bonds($n1, ...)

Returns the bonds with the given indices, or all by default.
Indices start from one, not from zero.

=cut

sub bonds {
    my $self = shift;
    if (@_) {
        my @bonds = map {$_ - 1} @_;
        @{$self->{bonds}}[@bonds];
    } else {
        @{$self->{bonds}};
    }
}

=item $mol->print(option => value...)

Convert the molecule to a string representation. If no options are given, 
a default YAML-like format is used (this may change in the future). Otherwise,
the format should be specified by using the C<format> option.

=cut

sub print {
    my $self = shift;
    my (%opts) = @_;
    my $ret;
    local $" = ""; #"

    if ($opts{format}) {
        return $self->formats($opts{format})->write_string($self, %opts);
    }
    # else use default printout 
    $ret = <<END;
$self->{id}:
    name: $self->{name}
END
    $ret .= "    attr:\n";
    $ret .= $self->print_attr(2);
    $ret .= "    atoms:\n";
    for my $a (@{$self->{atoms}}) { $ret .= $a->print(2) }
    $ret .= "    bonds:\n";
    for my $b (@{$self->{bonds}}) { $ret .= $b->print(2) }
    $ret;
}

=item $s = $mol->sprintf($format)

Format interesting molecular information in a concise way, as specified by
a printf-like format.

    %n - name
    %f - formula 
    %f{formula with format} - (note: right braces within
        the format should be escaped with a backslash)
    %s - SMILES representation
    %S - canonical SMILES representation
    %m - mass
    %8.3m - mass, formatted as %8.3f with core sprintf
    %q - formal charge
    %a - atom count
    %b - bond count
    %t - type
    %i - id
    %% - %

For example, if you want just about everything:

    $mol->sprintf("%s - %n (%f). %a atoms, %b bonds; "
        . "mass=%m; charge =%q; type=%t; id=%i");

Note that you have to C<use Chemistry::File::SMILES> before using C<%s> or
C<%S> on C<< $mol->sprintf >>.

=cut

sub sprintf {
    my ($mol, $format) = @_;
    no warnings 'uninitialized'; # don't care if some properties are undefined
    $format ||= "%f";
    $format =~ s/%%/\\%/g;              # escape %% with a \
    $format =~ s/(?<!\\)%f\{(.*?)(?<!\\)\}/$mol->formula($1)/eg; # %f{}
    $format =~ s/(?<!\\)%f/$mol->formula/eg;                    # %f
    $format =~ s/(?<!\\)%s/$mol->print(format=>'smiles')/eg;    # %s
    $format =~ s/(?<!\\)%S/$mol->print(format=>'smiles', unique => 1)/eg;    # %s
    $format =~ s/(?<!\\)%n/$mol->name/eg;                       # %n
    $format =~ s/(?<!\\)%(\d*\.?\d*)m/
        $1 ? sprintf "%$1f", $mol->mass : $mol->mass/eg;        # %m
    $format =~ s/(?<!\\)%q/$mol->charge/eg;                     # %q
    $format =~ s/(?<!\\)%a/$mol->atoms/eg;                      # %a
    $format =~ s/(?<!\\)%b/$mol->bonds/eg;                      # %b
    $format =~ s/(?<!\\)%t/$mol->type/eg;                       # %t
    $format =~ s/(?<!\\)%i/$mol->id/eg;                         # %i
    $format =~ s/\\(.)/$1/g;                             # other \ escapes
    $format;
}

=item $mol->printf($format)

Same as C<< $mol->sprintf >>, but prints to standard output automatically.
Used for quick and dirty molecular information dumping.

=cut

sub printf {
    my ($mol, $format) = @_;
    print $mol->sprintf($format);
}

=item Chemistry::Mol->parse($string, option => value...)

Parse the molecule encoded in C<$string>. The format should be specified
with the the C<format> option; otherwise, it will be guessed.

=cut

sub parse {
    my $self = shift;
    my $s = shift;
    my %opts = (mol_class => $self, @_);

    if ($opts{format}) {
        return $self->formats($opts{format})->parse_string($s, %opts);
    } else {
        croak "Parse does not support autodetection yet.",
            "Please specify a format.";
    }
    return;
}

=item Chemistry::Mol->read($fname, option => value ...)

Read a file and return a list of Mol objects, or croaks if there was a problem.
The type of file will be guessed if not specified via the C<format> option.

Note that only registered file readers will be used. Readers may be registered
using C<register_type()>; modules that include readers (such as
L<Chemistry::File::PDB>) usually register them automatically when they are
loaded.

Automatic decompression of gzipped files is supported if the L<Compress::Zlib>
module is installed. Files ending in .gz are assumed to be compressed;
otherwise it is possible to force decompression by passing the gzip => 1
option (or no decompression with gzip => 0).

=cut

sub read_mol { # for backwards compatibility
    my ($fname, $type) = shift;
    __PACKAGE__->read($fname, format => $type);
}

sub read {
    my $self = shift;
    my $fname = shift;
    my %opts = (mol_class => $self, @_);

    if ($opts{format}) {
        return $self->formats($opts{format})->parse_file($fname, %opts);
    } else { # guess format
        for my $type ($self->formats) {
            if ($self->formats($type)->file_is($fname)) {
                return $self->formats($type)->parse_file($fname, %opts);
            }
        }
    }
    croak "Couldn't guess format of file '$fname'";
}

=item $mol->write($fname, option => value ...)

Write a molecule file, or croak if there was a problem. The type of file will
be guessed if not specified via the C<format> option.

Note that only registered file formats will be used. 

Automatic gzip compression is supported if the IO::Zlib module is installed.
Files ending in .gz are assumed to be compressed; otherwise it is possible to
force compression by passing the gzip => 1 option (or no compression with gzip
=> 0). Specific compression levels between 2 (fastest) and 9 (most compressed)
may also be used (e.g., gzip => 9).

=cut

sub write {
    my ($self, $fname, %opts) = (@_);

    if ($opts{format}) {
        return $self->formats($opts{format})->write_file(@_);
    } else { # guess format
        for my $type ($self->formats) {
            if ($self->formats($type)->name_is($fname)) {
                return $self->formats($type)->write_file(@_);
            }
        }
    }
    croak "Couldn't guess format for writing file '$fname'";
}

=item Chemistry::Mol->file($file, option => value ...)

Create a L<Chemistry::File>-derived object for reading or writing to a file.
The object can then be used to read the molecules or other information in the
file.

This has more flexibility than calling C<< Chemistry::Mol->read >> when
dealing with multi-molecule files or files that have higher structure or that
have information that does not belong to the molecules themselves. For
example, a reaction file may have a list of molecules, but also general
information like the reaction name, yield, etc. as well as the classification
of the molecules as reactants or products. The exact information that is
available will depend on the file reader class that is being used. The
following is a hypothetical example for reading MDL rxnfiles.

    # assuming this module existed...
    use Chemistry::File::Rxn;

    my $rxn = Chemistry::Mol->file('test.rxn');
    $rxn->read;
    $name      = $rxn->name;
    @reactants = $rxn->reactants; # mol objects
    @products  = $rxn->products;
    $yield     = $rxn->yield;     # a number

Note that only registered file readers will be used. Readers may be registered
using register_type(); modules that include readers (such as
Chemistry::File::PDB) usually register them automatically.

=cut

sub file {
    my ($self,  $file, %opts) = @_;
    %opts = (mol_class => $self, %opts);

    if ($opts{format}) {
        return $self->formats($opts{format})->new(file => $file, 
            opts => \%opts);
    } else { # guess format
        for my $type ($self->formats) {
            if ($self->formats($type)->file_is($file)) {
                return $self->formats($type)->new(file => $file, 
                    opts => \%opts);
            }
        }
    }
    croak "Couldn't guess format of file '$file'";
}

=item Chemistry::Mol->register_format($name, $ref)

Register a file type. The identifier $name must be unique.  $ref is either a
class name (a package) or an object that complies with the L<Chemistry::File>
interface (e.g., a subclass of Chemistry::File).  If $ref is omitted, the
calling package is used automatically. More than one format can be registered
at a time, but then $ref must be included for each format (e.g.,
Chemistry::Mol->register_format(format1 => "package1", format2 => package2).

The typical user doesn't have to care about this function. It is used
automatically by molecule file I/O modules.

=cut

sub register_format {
    my $class = shift;
    if (@_ == 1) {
        $FILE_FORMATS{$_[0]} = caller;
        return;
    }
    my %opts = @_;
    $FILE_FORMATS{$_} = $opts{$_} for keys %opts;
}

=item Chemistry::Mol->formats

Returns a list of the file formats that have been installed by
register_type()

=cut

sub formats {
    my $self = shift;
    if (@_) {
        my ($type) = @_;
        my $file_class = $FILE_FORMATS{$type};
        unless ($file_class) {
            croak "No class installed for type '$type'";
        }
        return $file_class;
    } else {
        return sort keys %FILE_FORMATS;
    }
}

=item $mol->mass

Return the molar mass. This is just the sum of the masses of the atoms.  See
L<Chemistry::Atom>::mass for details such as the handling of isotopes.

=cut

sub mass {
    my ($self) = @_;
    my $mass = 0;
    for my $atom ($self->atoms) {
        $mass += $atom->mass;
    }
    $mass;
}

=item $mol->charge

Return the charge of the molecule. By default it returns the sum of the formal
charges of the atoms. However, it is possible to set an arbitrary charge by
calling C<< $mol->charge($new_charge) >>

=cut

sub charge {
    my ($self) = shift;
    if (@_) {
        $self->{charge} = shift;
        $self;
    } else {
        return $self->{charge} if defined $self->{charge};
        my $charge = 0;
        $charge += $_->formal_charge || 0 for $self->atoms;
        $charge;
    }
}

=item $mol->formula_hash

Returns a hash reference describing the molecular formula. For methane it would
return { C => 1, H => 4 }.

=cut

sub formula_hash {
    my ($self) = @_;
    my $formula = {};
    for my $atom ($self->atoms) {
        $formula->{$atom->symbol}++;
        $formula->{H} += $atom->hydrogens if $atom->hydrogens;
    }
    $formula;
}

=item $mol->formula($format)

Returns a string with the formula. The format can be specified as a printf-like
string with the control sequences specified in the L<Chemistry::File::Formula>
documentation.

=cut

sub formula {
    my ($self, $format) = @_;
    require Chemistry::File::Formula;
    $self->print(format => "formula", formula_format => $format);
}

=item my $mol2 = $mol->clone;

Makes a copy of a molecule. Note that this is a B<deep> copy; if your molecule
has a pointer to the rest of the universe, the entire universe will be cloned!

=cut

sub clone {
    my ($self) = @_;
    my $clone = dclone $self;
    $clone->_weaken if Storable->VERSION < 2.14;
    $clone;
}

=item my $mol2 = $mol->safe_clone;

Like clone, it makes a deep copy of a molecule. The difference is that the copy
is not "exact" in that new molecule and its atoms and bonds get assigned new
IDs. This makes it safe to combine cloned molecules. For example, this is an
error:

    # XXX don't try this at home!
    my $mol2 = Chemistry::Mol->combine($mol1, $mol1);
    # the atoms in $mol1 will clash

But this is ok:

    # the "safe clone" of $mol1 will have new IDs
    my $mol2 = Chemistry::Mol->combine($mol1, $mol1->safe_clone);

=cut

sub safe_clone {
    my ($mol) = @_;
    my $clone = $mol->clone;
    for ($clone, $clone->atoms, $clone->bonds) {
        $_->id($_->nextID);
    }
    $clone;
} 

sub _weaken {
    my ($self) = @_;
    for ($self->atoms, $self->bonds) {
        $_->_weaken;
    }
    $self;
}

=item ($distance, $atom_here, $atom_there) = $mol->distance($obj)

Returns the minimum distance to $obj, which can be an atom, a molecule, or a
vector. In scalar context it returns only the distance; in list context it
also returns the atoms involved. The current implementation for calculating
the minimum distance between two molecules compares every possible pair of
atoms, so it's not efficient for large molecules.

=cut

sub distance {
    my ($self, $other) = @_;
    if ($other->isa("Chemistry::Mol")) {
        my @atoms = $self->atoms;
        my $atom = shift @atoms or return; # need at least one atom
        my $closest_here = $atom;
        my ($min_length, $closest_there) = $atom->distance($other);
        for $atom (@atoms) {
            my ($d, $o) = $atom->distance($other);
            if ($d < $min_length) {
                ($min_length, $closest_there, $closest_here) = ($d, $o, $atom);
            }
        }
        return wantarray ? 
            ($min_length, $closest_here, $closest_there) : $min_length;
    } elsif ($other->isa("Chemistry::Atom")) {
        return $other->distance($self);
    } elsif ($other->isa("Math::VectorReal")) {
        return Chemistry::Atom->new(coords => $other)->distance($self);
    }
}

=item my $bigmol = Chemistry::Mol->combine($mol1, $mol2, ...)

=item $mol1->combine($mol2, $mol3, ...)

Combines several molecules in one bigger molecule. If called as a class method,
as in the first example, it returns a new combined molecule without altering
any of the parameters. If called as an instance method, as in the second
example, all molecules are combined into $mol1 (but $mol2, $mol3, ...) are not
altered. B<Note>: Make sure you don't combine molecules which contain atoms
with duplicate IDs (for example, if they were cloned).

=cut

# joins several molecules into one
sub combine {
    my ($self, @others) = @_;
    my $mol;
    if (ref $self) {
        $mol = $self;
    } else {
        $mol = $self->new;
    }
    for my $other (@others) {
        my $mol2 = $other->clone;
        for my $atom ($mol2->atoms) {
            $mol->add_atom($atom);
        }
        for my $bond ($mol2->bonds) {
            $mol->add_bond($bond);
        }
    }
    $mol;
}

=item my @mols = $mol->separate

Separates a molecule into "connected fragments". The original object is not
modified; the fragments are clones of the original ones. Example: if you have
ethane (H3CCH3) and you delete the C-C bond, you have two CH3 radicals within
one molecule object ($mol). When you call $mol->separate you get two molecules,
each one with a CH3.

=cut

# splits a molecule into connected fragments
# returns a list of molecules. Does not touch the original copy.
sub separate {
    my ($self) = @_;
    $self = $self->clone;
    $self->{_paint_tab} = {};
    my $color = 0;
    for my $atom ($self->atoms) {
        next if defined $self->{_paint_tab}{$atom->id};
        $self->_paint($atom, $color++);
    }
    my @mols;
    push @mols, $self->new for (1 .. $color);
    for my $atom ($self->atoms) {
        $mols[$self->{_paint_tab}{$atom->id}]->add_atom($atom);
    }
    for my $bond ($self->bonds) {
        $mols[$self->{_paint_tab}{$bond->id}]->add_bond($bond);
    }
    @mols;
}

# this method fills the _paint_tab attribute for every atom connected
# to the given start atom $atom with $color. Used for separating
# connected fragments. Uses a depth-first search
sub _paint {
    my ($self, $atom, $color) = @_;
    return if defined $self->{_paint_tab}{$atom->id};
    $self->{_paint_tab}{$atom->id} = $color;
    $self->{_paint_tab}{$_->id}    = $color for ($atom->bonds);
    for my $neighbor ($atom->neighbors) {
        $self->_paint($neighbor, $color);
    }
}

=item $mol->sprout_hydrogens

Convert all the implicit hydrogen atoms in the molecule to explicit atoms.
It does B<not> generate coordinates for the atoms.

=cut

sub sprout_hydrogens {
    my ($self) = @_;
    $_->sprout_hydrogens for $self->atoms;
}

=item $mol->collapse_hydrogens

Convert all the explicit hydrogen atoms in the molecule to implicit hydrogens.
(Exception: hydrogen atoms that are adjacent to a hydrogen atom are not
collapsed.)

=cut

sub collapse_hydrogens {
    my ($self) = @_;
    for my $atom (grep { $_->symbol ne 'H' } $self->atoms) {
        $atom->collapse_hydrogens;
    }
}

=item $mol->add_implicit_hydrogens

Use heuristics to figure out how many implicit hydrogens should each atom in 
the molecule have to satisfy its normal "organic" valence.

=cut

sub add_implicit_hydrogens {
    my ($self) = @_;
    $_->add_implicit_hydrogens for $self->atoms;
}

1;

=back

=head1 VERSION

0.35

=head1 SEE ALSO

L<Chemistry::Atom>, L<Chemistry::Bond>, L<Chemistry::File>,
L<Chemistry::Tutorial>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

