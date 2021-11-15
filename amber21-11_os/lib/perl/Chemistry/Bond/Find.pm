package Chemistry::Bond::Find;

$VERSION = '0.21';
our $DEBUG = 0;
# $Id: Find.pm,v 1.6 2004/06/16 19:36:29 itubert Exp $

=head1 NAME

Chemistry::Bond::Find - Detect bonds in a molecule from atomic 3D coordinates and assign formal bond orders

=head1 SYNOPSIS

    use Chemistry::Bond::Find ':all'; # export all available functions

    # $mol is a Chemistry::Mol object
    find_bonds($mol);
    assign_bond_orders($mol);

=head1 DESCRIPTION

This module provides functions for detecting the bonds in a molecule from its
3D coordinates by using simple cutoffs, and for guessing the formal bond
orders.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 FUNCTIONS

These functions may be exported, although nothing is exported by default.

=over

=cut

use strict;
use Chemistry::Mol;
use Exporter;
use Carp;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(find_bonds assign_bond_orders );
our %EXPORT_TAGS = (all => \@EXPORT_OK);

# table taken from
# http://environmentalchemistry.com/yogi/periodic/covalentradius.html
my %Covalent_Radius = (
    Ag => 1.34, Al => 1.18, Ar => 0.98, As => 1.20, At => 1.45, Au => 1.34,
    B  => 0.82, Ba => 1.98, Be => 0.90, Bi => 1.46, Br => 1.14, C  => 0.77,
    Ca => 1.74, Cd => 1.48, Ce => 1.65, Cl => 0.99, Co => 1.16, Cr => 1.18,
    Cs => 2.35, Cu => 1.17, Dy => 1.59, Er => 1.57, Eu => 1.85, F  => 0.72,
    Fe => 1.17, Ga => 1.26, Gd => 1.61, Ge => 1.22, H  => 0.32, He => 0.93,
    Hf => 1.44, Hg => 1.49, Ho => 1.58, I  => 1.33, In => 1.44, Ir => 1.27,
    K  => 2.03, Kr => 1.12, La => 1.69, Li => 1.23, Lu => 1.56, Mg => 1.36,
    Mn => 1.17, Mo => 1.30, N  => 0.75, Na => 1.54, Nb => 1.34, Nd => 1.64,
    Ne => 0.71, Ni => 1.15, O  => 0.73, Os => 1.26, P  => 1.06, Pb => 1.47,
    Pd => 1.28, Pm => 1.63, Po => 1.46, Pr => 1.65, Pt => 1.30, Rb => 2.16,
    Re => 1.28, Rh => 1.25, Ru => 1.25, S  => 1.02, Sb => 1.40, Sc => 1.44,
    Se => 1.16, Si => 1.11, Sm => 1.62, Sn => 1.41, Sr => 1.91, Ta => 1.34,
    Tb => 1.59, Tc => 1.27, Te => 1.36, Th => 1.65, Ti => 1.32, Tl => 1.48,
    Tm => 1.56, U  => 1.42, V  => 1.22, W  => 1.30, Xe => 1.31, Y  => 1.62,
    Yb => 1.74, Zn => 1.25, Zr => 1.45,
);

my $Default_Radius = 1.5; # radius for unknown elements

# I considered inlining this function, but the performance gain was minimal
# ( < 5 % ), so it's probably better to leave it here
# $opts->{cuttof} Hash Table of cutoffs. 
# Example: key = "H C", value = [0.76, 1.42]
sub are_bonded {
    my ($symbol_a, $symbol_b, $r, $opts) = @_;
    return $r < ($opts->{cuttoffs}{"$symbol_a $symbol_b"} ||= 
        (($Covalent_Radius{$symbol_a} || $Default_Radius) 
            + ($Covalent_Radius{$symbol_b} || $Default_Radius)) 
        * $opts->{tolerance});
}


=item C<find_bonds($mol, %options)>

Finds and adds the bonds in a molecule. Only use it in molecules that have no 
explicit bonds; for example, after reading a file with 3D coordinates but no
bond orders.

Available options:

=over

=item tolerance

Defaults to 1.1. Two atoms are considered to be bound if the distance between
them is less than the sum of their covalent radii multiplied by the tolerance.

=item margin

NOTE: in general setting this option is not recommended, unless you know what
you are doing. It is used by the space partitioning algorithm to determine the
"bucket size". It defaults to 2 * Rmax * tolerance, where Rmax is the largest
covalent radius among the elements found in the molecule. For example, if a
molecule has C, H, N, O, and I, Rmax = R(I) = 1.33, so the margin defaults to 2
* 1.33 * 1.1 = 2.926. This margin ensures that no bonds are missed by the
partitioning algorithm. 

Using a smaller value gives faster results, but at the risk of missing some
bonds. In this example, if you are certain that your molecule doesn't contain
I-I bonds (but it has C-I bonds), you can set margin to (0.77 + 1.33) * 1.1 =
2.31 and you still won't miss any bonds (0.77 is the radius of carbon).  This
only has a significant impact for molecules with a thousand atoms or more, but
it can reduce the execution time by 50% in some cases.

=item orders

If true, assign the bond orders after finding them, by calling
C<assign_bond_orders($mol, %opts)>.

=item bond_class

The class that will be used for creating the new bonds. The default is 
L<Chemistry::Bond>.

=back

=cut

# The new algorithm based on a suggestion by BrowserUK
sub find_bonds {
    my ($mol, %opts) = @_;
    %opts = (min_atoms => 20, tolerance => 1.1, %opts,
        cutoffs => {});  # set defaults
    my $margin = guess_margin($mol, \%opts);    
    %opts = (margin => $margin, %opts);
    my $grid = {};
    partition_space($mol, $grid, \%opts);
    find_bonds_grid($mol, $grid, \%opts);
    if ($opts{orders}) {
        assign_bond_orders($mol, %opts);
    }
}

use POSIX 'floor';
my $Y = 1000;
my $X = $Y * $Y;
my $Z = 1;
my $ORIGIN = int(($Y**3 + $Y**2 + $Y)/2);

# used by the new BrowserUK algorithm
sub partition_space {
    my ($mol, $grid, $opts) = @_;
    my $margin = $opts->{margin};
    for my $atom ($mol->atoms) {
        my (@vec) = $atom->coords->array;
        my (@norm_vec) = map { floor($_ / $margin) } @vec;
        my $n = $X * $norm_vec[0] + $Y * $norm_vec[1]
            + $norm_vec[2] + $ORIGIN;
        push @{$grid->{$n}}, $atom;
    }
}

# used by the new BrowserUK algorithm
sub find_bonds_grid {
    my ($mol, $grid, $opts) = @_;
    while (my ($n, $atoms) = each %$grid) {
        #print "Cell $n has ". @$atoms . " atoms\n";
        find_bonds_n2_one_set($mol, $atoms, $opts);
        for my $neigh_n (
            $n+$Z, 
            $n+$Z+$Y, $n+$Z-$Y, $n+$Z+$X, $n+$Z-$X, 
            $n+$Z+$Y+$X, $n+$Z+$Y-$X, $n+$Z-$Y+$X, $n+$Z-$Y-$X,
            $n+$Y, $n+$Y+$X, $n+$Y-$X,
            $n+$X,
            ) {
            if ($grid->{$neigh_n}) {
                find_bonds_n2_two_sets($mol, $atoms, $grid->{$neigh_n}, $opts);
            }
        }
    }
}

# used by both find_bonds variants to figure out the maximum cutoff
sub guess_margin {
    my ($mol, $opts) = @_;
    my $formula = $mol->formula_hash;
    my $max = 0;
    for my $elem (keys %$formula) {
        $max = $Covalent_Radius{$elem} if $Covalent_Radius{$elem} > $max;
    }
    $max *= 2 * $opts->{tolerance};
    #printf "MARGIN guessed at (%.2f)\n", $max;
    $max;
}

# brute-force N^2 algorithm
sub find_bonds_n2_one_set {
    my ($mol, $atoms, $opts) = @_;
    my $bond_class = $opts->{bond_class} || "Chemistry::Bond";
    for (my $i = 0; $i < @$atoms; ++$i) {
	for (my $j = $i + 1; $j < @$atoms; ++$j) {
	    my ($a1, $a2) = ($atoms->[$i], $atoms->[$j]);
	    if (are_bonded($a1->symbol, $a2->symbol, scalar $a1->distance($a2), $opts)) {
		$mol->add_bond($bond_class->new(atoms=>[$a1, $a2]));
	    }
	}
    }
}

# brute force N*M algorithm for finding bonds between to sets of atoms
sub find_bonds_n2_two_sets {
    my ($mol, $atoms1, $atoms2, $opts) = @_;
    my $bond_class = $opts->{bond_class} || "Chemistry::Bond";
    for my $a1 (@$atoms1) {
	for my $a2 (@$atoms2) {
	    if (are_bonded($a1->symbol, $a2->symbol, scalar $a1->distance($a2), $opts)) {
		$mol->add_bond($bond_class->new(atoms=>[$a1, $a2]));
	    }
	}
    }
}

# The recursive algorithm used in version 0.05, here for historical purposes
sub find_bonds_rec {
    my ($mol, %opts) = @_;
    %opts = (min_atoms => 20, tolerance => 1.1, %opts,
        cutoffs => {});  # set defaults
    my $margin = guess_margin($mol, \%opts);    
    %opts = (margin => $margin, %opts);
    _partition($mol, [$mol->atoms], 0, \%opts);
}

# partition space recursively, used by find_bonds_rec
sub _partition {
    my ($mol, $atoms, $dir, $opts) = @_;

    #printf "_partition(%s, $dir)\n", scalar(@$atoms);

    if (@$atoms < $opts->{min_atoms}) {
        #print "BOTTOM!\n";
        find_bonds_n2_one_set($mol, $atoms, $opts);
        return;
    }

    my $min = ($atoms->[0]->coords->array)[$dir];
    my $max = $min;
    my $center;
    for my $atom (@$atoms) {
        my $coord = ($atom->coords->array)[$dir]; 
        $center += $coord;
        $min = $coord if $coord < $min;
        $max = $coord if $coord > $min;
    }
    $center /= @$atoms;
    #printf "center($dir)=%.2f; range=(%.2f, %.2f)\n", $center, $min, $max;

    my @left = grep { ($_->coords->array)[$dir] < $center } @$atoms;
    my @right = grep { ($_->coords->array)[$dir] >= $center } @$atoms;
    _partition($mol, \@left, ($dir + 1) % 3, $opts);
    _partition($mol, \@right, ($dir + 1) % 3, $opts);

    # merge the interface between the two halves
    my $margin = $opts->{margin};
    my @left_margin = 
        grep { ($_->coords->array)[$dir] > $center - $margin } @left;
    my @right_margin = 
        grep { ($_->coords->array)[$dir] < $center + $margin } @right;
    find_bonds_n2_two_sets($mol, \@left_margin, \@right_margin, $opts);
}


=item C<assign_bond_orders($mol, %opts)>

Assign the formal bond orders in a molecule. The bonds must already be defined,
either by C<find_bonds> or because the molecule was read from a file that
included bonds but no bond orders. If the bond orders were already defined
(maybe the molecule came from a file that did include bond orders after all),
the original bond orders are erased and the process begins from scratch. Two
different algorithms are available, and may be selected by using the "method"
option:

    assign_bond_orders($mol, method => 'itub');
    assign_bond_orders($mol, method => 'baber');

=over

=item itub

This is the default if no method is specified. Developed from scratch by the
author of this module, this algorithm requires only the connection table
information, and it requires that all hydrogen atoms be explicit. It looks for
an atom with unsatisfied valence, increases a bond order, and then does the
same recursively for the neighbors. If everybody's not happy at the end, it
backtracks and tries another bond. The recursive part does not cover the whole
molecule, but only the contiguous region of "unhappy" atoms next to the
starting atom and their neighbors. This permits separating the molecule into
independent regions, so that if one is solved and there's a problem in another,
we don't have to backtrack to the first one.

The itub algorithm has the following additional options:

=over

=item use_coords

Although the algorithm does not I<require> 3D coordinates, it uses them by
default to improve the initial guesses of which bond orders should be
increased. To avoid using coordinates, add the C<use_coords> option with a
false value:

    assign_bond_orders($mol, use_coords => 0);

The results are the same most of the time, but using good coordinates improves
the results for complicated cases such as fused heteroaromatic systems.

=item scratch

If true, start the bond order assignment from scratch by assuming that all bond
orders are 1. If false, start from the current bond orders and try to fix the
unsatisfied valences. This option is true by default.

=back

=item baber

A bond order assignment algorithm based on Baber, J. C.; Hodgkin, E. E.
J. Chem. Inf. Comp. Sci. 1992, 32, 401-406 (with some interpretation).

This algorithm uses the 3D coordinates along with various cutoffs and
confidence formulas to guess the bond orders. It then tries to resolve
conflicts by looping through the atoms (but is not recursive or backtracking).
It does not require explicit hydrogens (although it's better when they are
available) because it was designed for use with real crystallographic data
which often doesn't have hydrogen atoms.

This method doesn't always give a good answer, especially for conjugated and
aromatic systems. The variation used in this module adds some random numbers to
resolve some ambiguities and break loops, so the results are not even entirely
deterministic (the 'itub' method is deterministic but the result may depend on
the input order of the atoms).

=back

=cut

sub assign_bond_orders  {
    my ($mol, %opts) = @_;
    if ($opts{method} and $opts{method} eq 'baber') {
        assign_bond_orders_baber($mol, %opts);
    } else {
        assign_bond_orders_itub($mol, %opts);
    }
}

####### Bond order assignment algorithm by Ivan Tubert-Brohman

# The "typical" valence that we expect an atom to have satisfied. If not
# given, a value of 1 is assumed.
my %MIN_VALENCES = ( O => 2, C => 4, S => 2, H => 1, N => 3, P => 3, Si => 2,
    F => 1, Cl => 1, Br => 1, I => 1 );

# $ALLOWED_INCREASES{$from}{$to} means that element $from is willing to 
# exceed its minimum valence by having a multiple bond to element $to.
my %ALLOWED_INCREASES = ( 
    Cl => { O => 7 },
    Br => { O => 7 },
    I  => { O => 7 },
    S  => { O => 6, C => 3, S => 4 },
    N  => { O => 5, C => 4, N => 4 },
    Si => { O => 4, C => 4 },
    O  => { C => 3, O => 3 },
    P  => { O => 5, C => 4 },
);

sub assign_bond_orders_itub {
    my ($mol, %opts) = @_;
    %opts = (use_coords => 1, scratch => 1, funny => {}, %opts);
    #my @funny_atoms;

    if ($opts{scratch}) {
        $_->order(1) for $mol->bonds;
    }
    for my $atom ($mol->atoms) {
        if (wants_more_bonds($atom)) {
            my $ret = make_happy($atom, \%opts, []);
            print "Atom $atom made happy? '$ret'\n" if $DEBUG;
            unless ($ret) { # atom is funny (has formal charge or unpaired e-)
                $opts{funny}{$atom} = 1; 
                #push @funny_atoms, $atom;
            }
        }
    }
}

sub wants_more_bonds {
    my ($atom) = @_;
    my $valence = $atom->valence;
    #$n_bonds += $_->order for $atom->bonds;

    my $ret = ($valence < ($MIN_VALENCES{$atom->symbol} || 1));
    print "    $atom wants more bonds? '$ret'\n" if $DEBUG;
    $ret;
}

# the heart of the algorithm; increase a bond order, and then do the same
# recursively for the neighbors. If everybody's not happy at the end, 
# backtrack and try another bond. Note that this search does not cover the
# whole molecule, but only the "unhappy" atoms and their neighbors. This
# permits separating the molecule into independent regions, so that if one
# is solved and there's a problem in another, we don't have to backtrack
# to the first one.
sub make_happy {
    my ($atom, $opts, $q) = @_;

    if (!$opts->{funny}{$atom} and wants_more_bonds($atom)) {
        push @$q, ($atom->neighbors); # the queue of atoms to be checked
        for my $bn (sorted_neighbors($atom, $opts)) {
            my ($nei, $bond) = @$bn{'to', 'bond'};
            # note: it would be better to find the bond that is most likely to
            # be increased by taking bond length into account
            if (accepts_more_bonds($nei, $atom)) {
                my $order = $bond->order;
                print "increasing $bond($atom-$nei)\n" if $DEBUG;
                $bond->order($order+1); # increase order, be happy
                # now make sure everybody's happy
                push @$q, $nei->neighbors($atom);   # our friend's neighbors
                                                    # need to be happy too
                if (make_happy($atom, $opts, $q)) { # are we happy now?
                    return 1;  # everybody's happy?
                } else {
                    # something's wrong, will have to backtrack
                    print "decreasing $bond($atom-$nei)\n" if $DEBUG;
                    $bond->order($order);
                    push @$q, $nei; # remember; neighbor's not happy either
                }
            }
        }
        # couldn't find happiness
    } else {
        my $next = shift @$q;
        unless ($next) { # no one left; everybody is happy!
            print "no more atoms left to check at $atom\n" if $DEBUG;
            return 1;
        }
        return (make_happy($next, $opts, $q)); # happy if next atom is happy
    }
    print "not happy at $atom\n" if $DEBUG;
    0; # not happy
}

sub sorted_neighbors {
    my ($atom, $opts) = @_;
    my $use_coords = $opts->{use_coords};

    print "sorting neighbors\n" if $DEBUG;
    return map {
            $_->{bn}
        } sort { 
            ($b->{wants} cmp $a->{wants}) # those who want go first
            || ($use_coords ? ($a->{len} <=> $b->{len}) : 0);
                # if they both want to the same degree, the shortest
                # bond goes first if we are using coords
        } map { 
            +{ 
                wants => wants_more_bonds($_->{to}), 
                bn => $_, 
                len => !$use_coords || $_->{bond}->length,
            }
        } $atom->bonds_neighbors;
}

sub accepts_more_bonds {
    my ($atom, $to) = @_;
    my ($symbol, $to_symbol) = ($atom->symbol, $to->symbol);

    #my $n_bonds = 0;
    #$n_bonds += $_->order for $atom->bonds;
    my $valence = $atom->valence;

    # not enough bonds even for the minimum valence?
    return 1 if ($valence < ($MIN_VALENCES{$atom->symbol} || 1));

    if ($valence < ($ALLOWED_INCREASES{$symbol}{$to_symbol} || 0)) {
        # make sure we are willing to make multiple bonds with this guy
        return 1;
    } else {
        return 0; # max. valence satisfied
    }
}


############
# Bond order assignment algorithm based on Baber, J. C.; Hodgkin, E. E.
# J. Chem. Inf. Comp. Sci. 1992, 32, 401-406

# this algorithm uses the 3D coordinates along with various cutoffs and
# confidence formulas to assign the bond orders. It does not require all
# explicit hydrogens (although it's better when they are available) because
# it's design for use with real crystallographic data which often doesn't 
# have hydrogen.

my %Valences = (
    C => [4], N => [3,4], O => [2],
    P => [4, 5], S => [2, 4, 6], As => [4, 5], Se => [2, 4, 6],
    Te => [2, 4, 6],
    F => [1], Cl => [1, 3], Br => [1, 3, 5], I => [1, 3, 5, 7],
    H => [1],
);

# note; this table should only have single-letter symbols
my %Bond_Orders = (
    "C C" => { is => 1.49, id => 1.31, it => 1.18, wsd => 1.38, wdt => 1.21 },
    "C N" => { is => 1.42, id => 1.32, it => 1.14, wsd => 1.34, wdt => 1.20 },
    "C O" => { is => 1.41, id => 1.22, wsd => 1.28 },
    "C S" => { is => 1.78, id => 1.68, wsd => 1.70 },
    "N N" => { is => 1.40, id => 1.22, wsd => 1.32 },
    "N O" => { is => 1.39, id => 1.22, wsd => 1.25 },
    "O S" => { is => 1.58, id => 1.45, wsd => 1.54 },
    #"O S" => { is => 1.58, id => 1.45, wsd => 1.51 }, # modified value
    "O P" => { is => 1.60, id => 1.48, wsd => 1.52 },
);

# add keys for reverse bond for ease of lookup
for (keys %Bond_Orders) {
    $Bond_Orders{ scalar reverse $_ } = $Bond_Orders{$_};
}

sub assign_bond_orders_baber {
    my ($mol, %opts) = @_;
    my $max_tries = $opts{max_tries} || 10;

    assign_initial_bonds($mol);
    assign_initial_coordinations($mol);
    my $tries = 0;
    while (resolve_conflicts($mol)) {
        last if $tries++ > $max_tries;
        print "try again\n" if $DEBUG;
    }
    $tries;
}

sub assign_initial_bonds {
    my ($mol) = @_;

    for my $bond ($mol->bonds) {
        #my %bond_has = map { ($_->symbol, 1 ) } $bond->atoms;
        my $symbols = join " ", map { $_->symbol } $bond->atoms;
        my ($order, $confidence);
        my $l_obs = $bond->length;
        if ($symbols =~ /\b(H|F|Cl|Br|I)\b/) {
            $order      = 1;
            $confidence = 10000;
        } elsif ($symbols =~ /\bSi\b/) {
            $order      = 1;
            $confidence = 20 * ($l_obs - 1.4); 
        } elsif ($symbols =~ /\bB\b/) {
            $order      = 1;
            $confidence = 20 * ($l_obs - 1.2); 
        } elsif (my $pars = $Bond_Orders{$symbols}) {
            if ($pars->{it}) { # may have triple bonds
                if ($l_obs > $pars->{wsd}) {
                    $order      = 1;
                    $confidence = ($l_obs - $pars->{wsd}) /
                        (2 * ($pars->{is} - $pars->{wsd}));
                } elsif ($pars->{id} < $l_obs and $l_obs < $pars->{wsd}) {
                    $order      = 2;
                    $confidence = ($pars->{wsd} - $l_obs) /
                        (2 * ($pars->{wsd} - $pars->{id}));
                } elsif ($pars->{id} > $l_obs and $l_obs > $pars->{wdt}) {
                    $order      = 2;
                    $confidence = ($l_obs - $pars->{wdt}) /
                        (2 * ($pars->{id} - $pars->{wdt}));
                } else {
                    $order      = 3;
                    $confidence = ($pars->{wdt} - $l_obs) /
                        (2 * ($pars->{wdt} - $pars->{it}));
                }
            } else { # only single and double
                if ($l_obs > $pars->{wsd}) {
                    $order      = 1;
                    $confidence = ($l_obs - $pars->{wsd}) /
                        (2 * ($pars->{is} - $pars->{wsd}));
                } else {
                    $order      = 2;
                    $confidence = ($pars->{wsd} - $l_obs) /
                        (2 * ($pars->{wsd} - $pars->{id}));
                }
            }
        } else { # unknown atom pair; let's assume it's always single
            $order      = 1;
            $confidence = 100_000;
        }

        $confidence *= 10; # normalization constant
        $bond->order($order);
        $bond->attr("bond-find/confidence", $confidence);
        print "bond $bond($symbols) len=$l_obs, order=$order, ",
            "conf=$confidence\n"
            if $DEBUG;
    }
}

sub assign_initial_coordinations {
    my ($mol) = @_;

    for my $atom ($mol->atoms) {
        my @neighbors = $atom->neighbors;
        my $a_obs = 0;
        my $n = 0;
        my ($max_conns, $confidence);
        for (my $i = 0; $i < @neighbors - 1; $i++) { 
            for (my $j = $i + 1; $j < @neighbors; $j++) { 
                my $angle = Chemistry::Atom::angle_deg($neighbors[$i],
                    $atom, $neighbors[$j]);
                # don't count linear angles for octahedral geometries
                $n++, $a_obs += $angle unless @neighbors > 2 and $angle > 150;
            }
        }
        $a_obs /= $n if $n;

        if ($n == 0) {
            $max_conns  = 1;        # linear
            $confidence = 10000;
            $n = 1;
        } elsif ($a_obs > 150) {
            $max_conns  = 2;        # linear
            $confidence = ($a_obs - 150) / 30;
        } elsif (120 > $a_obs and $a_obs >= 115) {
            $max_conns  = 3;        # trigonal
            $confidence = ($a_obs - 115) / 5;
        } elsif (150 > $a_obs and $a_obs >= 120) {
            $max_conns  = 3;        # trigonal
            $confidence = (150 - $a_obs) / 30;
        } elsif (109.5 > $a_obs and $a_obs >= 99) {
            $max_conns  = 4;        # tetrahedral
            $confidence = ($a_obs - 99) / 10.5;
        } elsif (115 > $a_obs and $a_obs >= 109.5) {
            $max_conns  = 4;        # tetrahedral
            $confidence = (115 - $a_obs) / 5.5;
        } elsif ($a_obs < 99) {
            $max_conns  = 6;        # octahedral
            $confidence = (9 - $a_obs) / 9;
        } else {
            confess("impossible coordination angle $a_obs!");
        }

        $confidence *= 100 * $n;
        if ($Valences{$atom->symbol}) {
            $atom->attr("bond-find/valence", $Valences{$atom->symbol}[0]);
        } else {
            $atom->attr("bond-find/valence", $max_conns);
        }
        $atom->attr("bond-find/confidence", $confidence);
        $atom->attr("bond-find/max_conns", $max_conns);
        print "Atom $atom has CN $max_conns with a conf. of $confidence\n"
            if $DEBUG;
    }
}

sub resolve_conflicts {
    my ($mol) = @_;
   
    my $changes = 0;
    for my $atom ($mol->atoms) {
        my $valence    = $atom->attr('bond-find/valence');
        my $max_conns  = $atom->attr('bond-find/max_conns');
        my $confidence = $atom->attr('bond-find/confidence');
        my $n_conns    = $atom->bonds;

        my $n_bonds    = 0;
        $n_bonds += $_->order for $atom->bonds;

        if ($n_conns > $valence) {
            my $next_valence = next_valence($atom->symbol, $valence);
            if ($next_valence) {
                print "increasing valence of $atom to $next_valence\n" 
                    if $DEBUG;
                $atom->attr('bond-find/valence', $next_valence);
                ++$changes, redo;
            } else {
                warn "too many conns $n_conns to $atom with valence $valence\n";
                next;
            }
        } elsif ($n_bonds > $valence) {
            my $next_valence = next_valence($atom->symbol, $valence);
            if ($next_valence) {
                print "increasing valence of $atom to $next_valence\n" 
                    if $DEBUG;
                $atom->attr('bond-find/valence', $next_valence);
                ++$changes, next;
            }
            my ($min_conf, $min_bond);
            for my $bond ($atom->bonds) {
                if ($bond->order > 1 and 
                    (not defined $min_conf 
                        or $bond->attr("bond-find/confidence") < $min_conf)) {
                    $min_conf = $bond->attr("bond-find/confidence");
                    $min_bond = $bond;
                }
            }
            my $new_order = $min_bond->order - 1;
            $min_bond->order($new_order);
            print "Decreasing order of $min_bond to $new_order\n" if $DEBUG;
            $min_bond->attr('bond-find/confidence', $min_conf*1.2+rand());
            ++$changes, next;
        } elsif ($n_bonds + $max_conns - $n_conns < $valence) {
            my ($min_conf, $min_bond);
            for my $bond ($atom->bonds) {
                if (($bond->order == 1 
                        or $atom->symbol =~ /^[CN]$/ and $bond->order == 2) and 
                    (not defined $min_conf 
                        or $bond->attr("bond-find/confidence") < $min_conf)) {
                    $min_conf = $bond->attr("bond-find/confidence");
                    $min_bond = $bond;
                }
            }

            if ($confidence > 95 or $min_conf < $confidence) {
                my $new_order = $min_bond->order + 1;
                $min_bond->order($new_order);
                print "Increasing order of $min_bond to $new_order\n" if $DEBUG;
                $min_bond->attr('bond-find/confidence', $min_conf*1.2+rand());
            } else {
                $max_conns++;
                $atom->attr("bond-find/max_conns", $max_conns);
                print "Increasing coord. num. of $min_bond to $max_conns\n" 
                    if $DEBUG;
            }
            ++$changes, next;
        }
    }
    $changes;
}

sub next_valence {
    my ($symbol, $current) = @_;

    if ($Valences{$symbol}) {
        for my $v (@{$Valences{$symbol}}) {
            return $v if $v > $current;
        }
    }
    return undef;
}

1;

=back

=head1 VERSION

0.21

=head1 TO DO

Some future version should let the user specify the desired cutoffs, and 
not always create a bond but call a user-supplied function instead. This way
these functions could be used for other purposes such as finding hydrogen bonds
or neighbor lists.

Add some tests.

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Bond>,
L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

The new C<find_bonds> algorithm was loosely based on a suggestion by BrowserUK
on perlmonks.org (L<http://perlmonks.org/index.pl?node_id=352838>).

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

