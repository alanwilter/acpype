package Chemistry::Atom;

$VERSION = '0.35';
# $Id: Atom.pm,v 1.41 2005/05/20 19:01:04 itubert Exp $

=head1 NAME

Chemistry::Atom - Chemical atoms as objects in molecules

=head1 SYNOPSIS

    use Chemistry::Atom;

    my $atom = new Chemistry::Atom(
        id => 'a1',
        coords => [$x, $y, $z],
        symbol => 'Br'
    );

    print $atom->print;

=head1 DESCRIPTION

This module includes objects to describe chemical atoms. 
An atom is defined by its symbol and its coordinates, among other attributes.
Atomic coordinates are described by a Math::VectorReal
object, so that they can be easily used in vector operations.

=head2 Atom Attributes

In addition to common attributes such as id, name, and type, 
atoms have the following attributes, which are accessed or
modified through methods defined below: bonds, coords, internal_coords,
Z, symbol, etc.

In general, to get the value of a property, use $atom->method without
any parameters. To set the value, use $atom->method($new_value). When setting
an attribute, the accessor returns the atom object, so that accessors can be
chained:

    $atom->symbol("C")->name("CA")->coords(1,2,3);

=cut

# Considering to add the following attributes:
# mass_number (A)

use 5.006;
use strict;
use warnings;
use Scalar::Util 'weaken';
use Math::VectorReal qw(O vector);
use Math::Trig;
use Carp;
use base qw(Chemistry::Obj Exporter);
use List::Util qw(first);

our @EXPORT_OK = qw(distance angle dihedral angle_deg dihedral_deg);
our %EXPORT_TAGS = (
      all  => \@EXPORT_OK,
);



my $N = 0; # Number of atoms created so far, used to generate default IDs.

our @ELEMENTS = qw(
    n
    H                                                                   He
    Li  Be                                          B   C   N   O   F   Ne
    Na  Mg                                          Al  Si  P   S   Cl  Ar
    K   Ca  Sc  Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr
    Rb  Sr  Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe
    Cs  Ba
        La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb
            Lu  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn
    Fr  Ra
        Ac  Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No
            Lr  Rf  Db  Sg  Bh  Hs  Mt  Ds  Uuu Uub Uut Uuq Uup Uuh Uus Uuo
);

our %ELEMENTS;
for (my $i = 1; $i < @ELEMENTS; ++$i){
    $ELEMENTS{$ELEMENTS[$i]} = $i;
}
$ELEMENTS{D} = $ELEMENTS{T} = 1;

my %Atomic_masses = (
   "H" => 1.00794, "D" => 2.014101, "T" => 3.016049, "He" => 4.002602,
   "Li" => 6.941, "Be" => 9.012182, "B" => 10.811, "C" => 12.0107,
   "N" => 14.00674, "O" => 15.9994, "F" => 18.9984032, "Ne" => 20.1797,
   "Na" => 22.989770, "Mg" => 24.3050, "Al" => 26.981538, "Si" => 28.0855,
   "P" => 30.973761, "S" => 32.066, "Cl" => 35.4527, "Ar" => 39.948,
   "K" => 39.0983, "Ca" => 40.078, "Sc" => 44.955910, "Ti" => 47.867,
   "V" => 50.9415, "Cr" => 51.9961, "Mn" => 54.938049, "Fe" => 55.845,
   "Co" => 58.933200, "Ni" => 58.6934, "Cu" => 63.546, "Zn" => 65.39,
   "Ga" => 69.723, "Ge" => 72.61, "As" => 74.92160, "Se" => 78.96,
   "Br" => 79.904, "Kr" => 83.80, "Rb" => 85.4678, "Sr" => 87.62,
   "Y" => 88.90585, "Zr" => 91.224, "Nb" => 92.90638, "Mo" => 95.94,
   "Tc" => 98, "Ru" => 101.07, "Rh" => 102.90550, "Pd" => 106.42,
   "Ag" => 107.8682, "Cd" => 112.411, "In" => 114.818, "Sn" => 118.710,
   "Sb" => 121.760, "Te" => 127.60, "I" => 126.90447, "Xe" => 131.29,
   "Cs" => 132.90545, "Ba" => 137.327, "La" => 138.9055, "Ce" => 140.116,
   "Pr" => 140.90765, "Nd" => 144.24, "Pm" => 145, "Sm" => 150.36,
   "Eu" => 151.964, "Gd" => 157.25, "Tb" => 158.92534, "Dy" => 162.50,
   "Ho" => 164.93032, "Er" => 167.26, "Tm" => 168.93421, "Yb" => 173.04,
   "Lu" => 174.967, "Hf" => 178.49, "Ta" => 180.9479, "W" => 183.84,
   "Re" => 186.207, "Os" => 190.23, "Ir" => 192.217, "Pt" => 195.078,
   "Au" => 196.96655, "Hg" => 200.59, "Tl" => 204.3833, "Pb" => 207.2,
   "Bi" => 208.98038, "Po" => 209, "At" => 210, "Rn" => 222,
   "Fr" => 223, "Ra" => 226, "Ac" => 227, "Th" => 232.038,
   "Pa" => 231.03588, "U" => 238.0289, "Np" => 237, "Pu" => 244,
   "Am" => 243, "Cm" => 247, "Bk" => 247, "Cf" => 251,
   "Es" => 252, "Fm" => 257, "Md" => 258, "No" => 259,
   "Lr" => 262, "Rf" => 261, "Db" => 262, "Sg" => 266,
   "Bh" => 264, "Hs" => 269, "Mt" => 268, "Ds" => 271,
);

=head1 METHODS

=over 4

=item Chemistry::Atom->new(name => value, ...)

Create a new Atom object with the specified attributes.

=cut

sub new {
    my $class = shift;
    my %args = @_;

    my $self = bless {
        id => $class->nextID(),
        coords => vector(0, 0, 0),
        Z => 0,
        symbol => '',
        bonds => [],
    }, $class;

    $self->$_($args{$_}) for (keys %args);
    $self;
}

sub nextID {
    "a".++$N; 
}

sub reset_id {
    $N = 0; 
}


=item $atom->Z($new_Z)

Sets and returns the atomic number (Z). If the symbol of the atom doesn't
correspond to a known element, Z = undef.

=cut

sub Z {
    my $self = shift;

    if(@_) {
        $self->{symbol} = $ELEMENTS[$_[0]];
        $self->{Z} = $_[0];
        return $self;
    } else {
        return $self->{Z};
    }
}


=item $atom->symbol($new_symbol)

Sets and returns the atomic symbol.

=cut

sub symbol {
    my $self = shift;

    if(@_) {
        $_[0] =~ s/ //g;
        $self->{Z} = $ELEMENTS{$_[0]};
        $self->{symbol} = $_[0];
        return $self;
    } else {
        return $self->{symbol};
    }
}

=item $atom->mass($new_mass)

Sets and returns the atomic mass in atomic mass units. Any arbitrary mass may
be set explicitly by using this method. However, if no mass is set explicitly
and this method is called as an accessor, the return value is the following:

1) If the mass number is undefined (see the mass_number method below), the
relative atomic mass from the 1995 IUPAC recommendation is used. (Table stolen
from the Chemistry::MolecularMass module by Maksim A.  Khrapov).

2) If the mass number is defined and the L<Chemistry::Isotope> module is 
available and it knows the mass for the isotope, the exact mass of the isotope
is used; otherwise, the mass number is returned.

=cut

sub mass {
    my $self = shift;
    if (@_) {
        $self->{mass} = shift;
        return $self;
    } else {
        if (defined $self->{mass}) {
            return $self->{mass};
        } elsif (defined $self->{mass_number}) {
            if (eval { require Chemistry::Isotope } and 
                my $m = Chemistry::Isotope::isotope_mass(
                    $self->{mass_number}, $self->{Z})
            ) {
                return $m;
            } else {
                return $self->{mass_number};
            }
        } else {
            return $Atomic_masses{$self->symbol};
        }
    }
}

=item $atom->mass_number($new_mass_number)

Sets or gets the mass number. The mass number is undefined unless is 
set explicitly (this module does not try to guess a default mass number based
on the natural occuring isotope distribution).

=cut

Chemistry::Obj::accessor('mass_number');

=item $atom->coords

    my $vector = $atom->coords;  # get a Math::VectorReal object
    $atom->coords($vector);      # set a Math::VectorReal object 
    $atom->coords([$x, $y, $z]); # also accepts array refs 
    $atom->coords($x, $y, $z);   # also accepts lists

Sets or gets the atom's coordinates. It can take as a parameter a
Math::VectorReal object, a reference to an array, or the list of coordinates.

=cut

sub coords {
    my $self = shift;

    if(@_) {
        if (UNIVERSAL::isa($_[0], "Math::VectorReal")) {
            $self->{coords} = $_[0];
        } elsif (ref $_[0] eq "ARRAY") {
            $self->{coords} = vector(@{$_[0]});
        } else {
            $self->{coords} = vector(@_);
        }
    } else {
        return $self->{coords};
    }
    $self;
}

=item $atom->internal_coords

    # get a Chemistry::InternalCoords object
    my $ic = $atom->internal_coords;      

    # set a Chemistry::InternalCoords object 
    $atom->internal_coords($vic);         

    # also accepts array refs 
    $atom->internal_coords([8, 1.54, 7, 109.47, 6, 120.0]);   

    # also accepts lists
    $atom->internal_coords(8, 1.54, 7, 109.47, 6, 120.0);    

Sets or gets the atom's internal coordinates. It can take as a parameter a
Chemistry::InternalCoords object, a reference to an array, or the list of
coordinates. In the last two cases, the list elements are the following: atom 
number or reference for distance, distance, atom number or reference for angle,
angle in degrees, atom number or reference for dihedral, dihedral in degrees.

=cut

sub internal_coords {
    my $self = shift;

    if(@_) {
        if (UNIVERSAL::isa($_[0], "Chemistry::InternalCoords")) {
            $self->{internal_coords} = $_[0];
        } elsif (ref $_[0] eq "ARRAY") {
            require Chemistry::InternalCoords;
            $self->{internal_coords} = 
                Chemistry::InternalCoords->new($self, @{$_[0]});
        } else {
            require Chemistry::InternalCoords;
            $self->{internal_coords} = 
                Chemistry::InternalCoords->new($self, @_);
        }
    } else {
        return $self->{internal_coords};
    }
    $self;
}

=item $atom->x3, $atom->y3, $atom->z3

Get the x, y or z 3D coordinate of the atom. This methods are just accessors
that don't change the coordinates. $atom->x3 is short for 
($atom->coords->array)[0], and so on.

=cut

sub x3 { ($_[0]->coords->array)[0] }
sub y3 { ($_[0]->coords->array)[1] }
sub z3 { ($_[0]->coords->array)[2] }

=item $atom->formal_charge($charge)

Set or get the formal charge of the atom.

=cut

Chemistry::Obj::accessor('formal_charge');

=item $atom->formal_radical($radical)

Set or get the formal radical multiplicity for the atom.

=cut

Chemistry::Obj::accessor('formal_radical');

=item $atom->implicit_hydrogens($h_count)

Set or get the number of implicit hydrogen atoms bonded to the atom.

=cut

sub implicit_hydrogens { shift->hydrogens(@_) }

=item $atom->hydrogens($h_count)

Set or get the number of implicit hydrogen atoms bonded to the atom
(DEPRECATED: USE C<implicit_hydrogens> INSTEAD).

=cut

Chemistry::Obj::accessor('hydrogens');

=item $atom->total_hydrogens($h_count)

Get the total number of hydrogen atoms bonded to the atom (implicit +
explicit).

=cut

sub total_hydrogens {
    my ($self) = @_;
    no warnings 'uninitialized';
    $self->hydrogens + grep { $_->symbol eq 'H' } $self->neighbors;
}

=item $atom->sprout_hydrogens

Convert all the implicit hydrogens for this atom to explicit hydrogens. Note:
it does B<not> generate coordinates for the new atoms.

=cut

sub sprout_hydrogens {
    my ($self) = @_;
    for (1 .. $self->implicit_hydrogens || 0) {
        $self->parent->new_bond(
            atoms => [$self, $self->parent->new_atom(symbol => 'H')]);
    }
    $self->implicit_hydrogens(0);
}

=item $atom->collapse_hydrogens

Delete neighboring hydrogen atoms and add them as implicit hydrogens for this
atom.

=cut

sub collapse_hydrogens {
    my ($self) = @_;
    no warnings 'uninitialized';
    my $implicit = $self->implicit_hydrogens;
    for my $nei ($self->neighbors) {
        $nei->delete, $implicit++ if $nei->symbol eq 'H';
    }
    $self->implicit_hydrogens($implicit);
}

my %VALENCE_TABLE = (
    Br => 1, Cl => 1, B => 3, C => 4, N => 3, O => 2, P => 3, S => 2, 
    F => 1, I => 1,
);

# to make it easier to test
sub _calc_implicit_hydrogens {
    my ($self, $symbol, $valence, $charge, $radical) = @_;
    no warnings 'uninitialized';

    my $h_count = $VALENCE_TABLE{$symbol} - $valence;
    # should account for non-kekulized aromatic bonds

    # some common charge situations
    if (($symbol =~ /^[NOS]$/) && $charge == -1) {
        $h_count--;
    } elsif ($symbol =~ /^[NOSP]$/ && $charge == 1) {
        $h_count++;
    } elsif ($symbol eq 'C' && $charge) {
        $h_count--;
    } elsif ($symbol eq 'B' && $charge == -1) {
        $h_count++;
    }

    # some common radical situations
    if ($radical == 1 or $radical == 3) {
        $h_count -=2;
    } elsif ($radical == 2) {
        $h_count--;
    }

    $h_count = 0 if $h_count < 0;
    $h_count;
}

=item $atom->calc_implicit_hydrogens

Use heuristics to figure out how many implicit hydrogens should the atom have
to satisfy its normal "organic" valence. Returns the number of hydrogens but
does not affect the atom object.

=cut

sub calc_implicit_hydrogens {
    my ($self) = @_;
    $self->_calc_implicit_hydrogens(
        $self->symbol, $self->explicit_valence, 
        $self->formal_charge, $self->formal_radical,
    );
}

=item $atom->add_implicit_hydrogens

Similar to calc_implicit_hydrogens, but it also sets the number of implicit
hydrogens in the atom to the new calculated value. Equivalent to

    $atom->implicit_hydrogens($atom->calc_implicit_hydrogens);

It returns the atom object.

=cut

sub add_implicit_hydrogens {
    my ($self) = @_;
    my $h_count = $self->calc_implicit_hydrogens;
    $self->implicit_hydrogens($h_count);
}

=item $atom->aromatic($bool)

Set or get whether the atom is considered to be aromatic. This property may be
set arbitrarily, it doesn't imply any kind of "intelligent aromaticity
detection"! (For that, look at the L<Chemistry::Ring> module).

=cut

Chemistry::Obj::accessor('aromatic');

=item $atom->valence

Returns the sum of the bond orders of the bonds in which the atom participates,
including implicit hydrogens (which are assumed to have bond orders of one).

=cut

sub valence {
    my ($self) = @_;
    my $valence = 0;
    $valence += $_->order for $self->bonds;
    $valence += $self->hydrogens || 0;
    $valence;
}

=item $atom->explicit_valence

Like C<valence>, but excluding implicit hydrogen atoms. To get the raw number
of bonds, without counting bond orders, call $atom->bonds in scalar context.

=cut

sub explicit_valence {
    my ($self) = @_;
    my $valence = 0;
    $valence += $_->order for $self->bonds;
    $valence;
}

# this method is for internal use only; called by $mol->add_bond
sub add_bond {
    my $self = shift;
    my $bond = shift;
    my %seen;
    #return if first { $_ eq $bond } @{$self->{bonds}}; 

    for my $atom (@{$bond->{atoms}}){ #for each atom...
        if ($atom ne $self) {
            my $b = {to=>$atom, bond=>$bond};
            weaken($b->{to});
            weaken($b->{bond});
            push @{$self->{bonds}}, $b;
        }
    }
}

# make sure the atom doesn't cause circular references
sub _weaken {
    my $self = shift;
    for my $b (@{$self->{bonds}}) {
        weaken($b->{to});
        weaken($b->{bond});
    }
    weaken($self->{parent});
}

# This method is private. Bonds should be deleted from the 
# mol object. These methods should only be called by 
# $bond->delete_atoms, which is called by $mol->delete_bond
sub delete_bond {
    my ($self, $bond) = @_;
    $self->{bonds} = [ grep { $_->{bond} ne $bond } @{$self->{bonds}} ];
}

=item $atom->delete

Calls $mol->delete_atom($atom) on the atom's parent molecule.

=cut

sub delete {
    my ($self) = @_;
    $self->{parent}->_delete_atom($self);
}

=item $atom->parent

Returns the atom's containing object (the molecule to which the atom belongs).
An atom can only have one parent.

=cut

sub parent {
    my $self = shift;
    if (@_) {
        ($self->{parent}) = @_;
        weaken($self->{parent});
        $self;
    } else {
        $self->{parent};
    }
}

=item $atom->neighbors($from)

Return a list of neighbors. If an atom object $from is specified, it will be
excluded from the list (this is useful if an atom wants to know who its 
neighbor's neighbors are, without counting itself).

=cut

sub neighbors {
    my $self = shift;
    my $from = shift;
    my @ret = ();

    for my $b (@{$self->{bonds}}) {
        push @ret, $b->{to} unless $from && $b->{to} eq $from;
    }
    @ret;
}

=item $atom->bonds($from)

Return a list of bonds. If an atom object $from is specified, bonds to
that atom will be excluded from the list.

=cut

sub bonds {
    my $self = shift;
    my $from = shift;
    my @ret = ();

    for my $b (@{$self->{bonds}}) {
        push @ret, $b->{bond} unless $from && $b->{to} eq $from;
    }
    @ret;
}

=item $atom->bonds_neighbors($from)

Return a list of hash references, representing the bonds and neighbors from the
atom. If an atom object $from is specified, it will be excluded from the list.
The elements of the hash are 'to', and atom reference, and 'bond', a bond
reference. For example, 

    for my $bn ($atom->bonds_neighbors) {
        print "bond $bn->{bond} point to atom $bn->{to}\n";
    }

=cut

sub bonds_neighbors {
    my $self = shift;
    my $from = shift;
    my @ret = ();

    for my $b (@{$self->{bonds}}) {
        push @ret, {%$b} unless $from && $b->{to} eq $from;
    }
    @ret;
}

=item ($distance, $closest_atom) = $atom->distance($obj)

Returns the minimum distance to $obj, which can be an atom, a molecule, or a
vector. In scalar context it returns only the distance; in list context it
also returns the closest atom found. It can also be called as a function,
Chemistry::Atom::distance (which can be exported).

=cut

sub distance {
    my $self = shift;
    my $obj = shift;
    my $min_length;
    my $closest_atom = $obj;

    if ($obj->isa('Chemistry::Atom')) {
        my $v = $self->coords - $obj->coords;
        $min_length = $v->length;
    } elsif ($obj->isa('Math::VectorReal')) {
        my $v = $self->coords - $obj;
        $min_length = $v->length;
    } elsif ($obj->isa('Chemistry::Mol')) {
        my @atoms = $obj->atoms;
        my $a = shift @atoms or return undef; # ensure there's at least 1 atom
        $min_length = $self->distance($a);
        $closest_atom = $a;
        for $a (@atoms) {
            my $l = $self->distance($a);
            $min_length = $l, $closest_atom = $a if $l < $min_length;
        }
    } else {
        croak "atom->distance() undefined for objects of type '", ref $obj,"'";
    }
    wantarray ? ($min_length, $closest_atom) : $min_length;
}

=item $atom->angle($atom2, $atom3)

Returns the angle in radians between the atoms involved. $atom2 is the atom in
the middle. Can also be called as Chemistry::Atom::angle($atom1, $atom2,
$atom3). This function can be exported. Note: if you override this method,
you may also need to override angle_deg or strange things may happen.

=cut

# $a2 is the one in the center
sub angle {
    @_ == 3 or croak "Chemistry::Atom::angle requires three atoms!\n";
    my @c;
    for my $a (@_) { # extract coordinates
        ref $a or croak "Chemistry::Atom::angle: $a is not an object";
        push @c, $a->isa("Chemistry::Atom") ? $a->coords :
            $a->isa("Math::VectorReal") ? $a : 
                croak "angle: $a is neither an atom nor a vector!\n";
    }
    my $v1 = $c[0] - $c[1];
    my $v2 = $c[2] - $c[1];
    my $l = ($v1->length * $v2->length) or return 0;
    acos(($v1 . $v2) / $l);
}

=item $atom->angle_deg($atom2, $atom3)

Same as angle(), but returns the value in degrees. May be exported.

=cut

sub angle_deg {
    rad2deg(angle(@_));
}

=item $atom->dihedral($atom2, $atom3, $atom4)

Returns the dihedral angle in radians between the atoms involved.  Can also be
called as Chemistry::Atom::dihedral($atom1, $atom2, $atom3, $atom4). May be
exported. Note: if you override this method, you may also need to override 
dihedral_deg and angle or strange things may happen.

=cut

sub dihedral {
    @_ == 4 or croak "Chemistry::Atom::dihedral requires four atoms!\n";
    my @c;
    for my $a (@_) { # extract coordinates
        push @c, $a->isa("Chemistry::Atom") ? $a->coords :
            $a->isa("Math::VectorReal") ? $a : 
                croak "angle: $a is neither an atom nor a vector!\n";
    }
    my $v1 = $c[0] - $c[1];
    my $v2 = $c[2] - $c[1];
    my $v3 = $c[3] - $c[2];
    my $x1 = $v1 x $v2;
    my $x2 = $v3 x $v2;
    my $abs_dih = angle($x1, O(), $x2);
    $v1 . $x2 > 0 ? $abs_dih : -$abs_dih;
}

=item $atom->dihedral_deg($atom2, $atom3, $atom4)

Same as dihedral(), but returns the value in degrees. May be exported.

=cut

sub dihedral_deg {
    rad2deg(dihedral(@_));
}

=item $atom->print

Convert the atom to a string representation (used for debugging).

=cut

sub print {
    my $self = shift;
    my ($indent) = @_;

    no warnings; 
    $indent ||= 0;
    my $bonds = join " ", map {$_->id} $self->bonds;
    my $neighbors = join " ", map {$_->id} $self->neighbors;
    my $coords = $self->{coords}->stringify(
    'x3: %g
    y3: %g
    z3: %g'
    );

    my $ret = <<EOF;
$self->{id}:
    symbol: $self->{symbol}
    name  : $self->{name}
    $coords
    formal_charge: $self->{formal_charge}
    bonds: "$bonds"
    neighbors: "$neighbors"
EOF
    $ret .= "    attr:\n";
    $ret .= $self->print_attr($indent+2);
    $ret =~ s/^/"    "x$indent/gem;
    $ret;
}

=item my $info = $atom->sprintf($format)

Format interesting atomic information in a concise way, as specified by
a printf-like format.

    %s - symbol
    %Z - atomic number
    %n - name
    %q - formal charge
    %h - implicit hydrogen count
    %v - valence
    %i - id
    %8.3m - mass, formatted as %8.3f with core sprintf
    %8.3x - x coordinate, formatted as %8.3f with core sprintf
    %8.3y - y coordinate, formatted as %8.3f with core sprintf
    %8.3z - z coordinate, formatted as %8.3f with core sprintf
    %% - %

=cut

sub sprintf {
    my ($atom, $format) = @_;
    no warnings 'uninitialized'; # don't care if some properties are undefined
    $format ||= "%f";
    $format =~ s/%%/\\%/g;              # escape %% with a \
    $format =~ s/(?<!\\)%q/$atom->formal_charge || 0/eg;        # %q
    $format =~ s/(?<!\\)%s/$atom->symbol/eg;                    # %s
    $format =~ s/(?<!\\)%Z/$atom->Z/eg;                         # %Z
    $format =~ s/(?<!\\)%n/$atom->name/eg;                      # %n
    $format =~ s/(?<!\\)%h/$atom->hydrogens/eg;                 # %h
    $format =~ s/(?<!\\)%v/$atom->valence/eg;                   # %v
    $format =~ s/(?<!\\)%(\d*\.?\d*)m/
        $1 ? sprintf "%$1f", $atom->mass : $atom->mass/eg;      # %m
    $format =~ s/(?<!\\)%(\d*\.?\d*)x/
        $1 ? sprintf "%$1f", $atom->x3 : $atom->x3/eg;          # %x
    $format =~ s/(?<!\\)%(\d*\.?\d*)y/
        $1 ? sprintf "%$1f", $atom->y3 : $atom->y3/eg;          # %y
    $format =~ s/(?<!\\)%(\d*\.?\d*)z/
        $1 ? sprintf "%$1f", $atom->z3 : $atom->z3/eg;          # %z
    $format =~ s/(?<!\\)%i/$atom->id/eg;                        # %i
    $format =~ s/\\(.)/$1/g;                             # other \ escapes
    $format;
}

=item $atom->printf($format)

Same as $atom->sprintf, but prints to standard output automatically. Used
for quick and dirty atomic information dumping.

=cut

sub printf {
    my ($atom, $format) = @_;
    print $atom->sprintf($format);
}

1;

=back

=head1 VERSION

0.35

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Bond>, 
L<Math::VectorReal>, L<Chemistry::Tutorial>,
L<Chemistry::InternalCoords>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

