package Chemistry::Reaction;
$VERSION = '0.02';

=head1 NAME

Chemistry::Reaction - Explicit chemical reactions

=head1 SYNOPSIS

    use Chemistry::Reaction;
    use Chemistry::File::SMILES;

    my $s = Chemistry::Pattern->parse('C=CC=C.C=C', format=>'smiles');
    my $p = Chemistry::Pattern->parse('C1=CCCCC1', format=>'smiles');
    my %m;
    for (my $i = 1; $i <= $s->atoms; $i++) {
      $m{$s->atoms($i)} = $p->atoms($i);
    }
    my $r = Chemistry::Reaction->new($s, $p, \%m);

=head1 DESCRIPTION

This package, along with Chemistry::Pattern, provides an
implementation of explicit chemical reactions.

An explicit chemical reaction is a representation of the
transformation that takes place in a given chemical reaction. In an
explicit chemical reaction, a substrate molecule is transformed into a
product molecule by breaking existing bonds and creating new bonds
between atoms.

The representation of an explicit chemical reaction is a molecule in
which the order of a bond before the chemical reaction is
distinguished from the order of the bond after the chemical
reaction. Thus, the breaking of an existing bond is represented by one
of the following before/after pairs:

  3/2, 2/1, 1/0 (breaking of a single bond or reduce order by one)
       3/1, 2/0 (breaking of a double bond or reduce order by two)
            3/0 (breaking of a triple bond)

The creation of a new bond is represented by one of the following
before/after pairs:

  0/1, 1/2, 2/3 (creation of a single bond or increase order by one)
       0/2, 1/3 (creation of a double bond or increase order by two)
            0/3 (creation of a triple bond)

An explicit chemical reaction $react can be forward or reverse applied
once to a molecule $mol at the first subgraph of $mol found which is
isomorphic to the substrate or product of $react:

    my $subst = $react->substrate;
    if ($subst->match($mol)) {
      $react->forward($mol, $subst->atom_map);
    }

Also, an explicit chemical reaction $react can be forward or reverse
applied once to a molecule $mol at each subgraph of $mol which is
isomorphic to the substrate or product of $react:

    my $subst = $react->substrate;
    my @products;
    while ($subst->match($mol)) {
      my $new_mol = $mol->clone; # start from a fresh molecule
      my @map = $subst->atom_map;
      # translate atom map to the clone
      my @m = map { $new_mol->by_id($_->id) } @map;
      $react->forward($new_mol, @m);
      push @products, $new_mol;
    }

Furthermore, an explicit chemical reaction $react can be forward or
reverse applied as long as possible to a molecule $mol at the first
subgraph of $mol found which is isomorphic to the substrate or product
of $react:

    my $subst = $react->substrate;
    while ($subst->match($mol)) {
      $react->forward($mol, $subst->atom_map);
    }

=cut

use 5.006;
use strict;
use warnings;
use base qw(Chemistry::Pattern);

=head1 METHODS

=over 4

=item Chemistry::Reaction->new($subst, $prod, \%map)

Create a new Reaction object that describes the transformation of the
$subst substrate into the $prod product, according to the %map mapping
of substrate atoms to product atoms.

=cut

sub new {
  my ($class, $substrate, $product, $mapping, %args) = @_;

  die sprintf(
    "$class substrate and product must coincide on atoms (%s ne %s)\n", 
        $substrate->formula, $product->formula) 
    if $substrate->formula ne $product->formula;

  my $order1 = 0;
  foreach my $bond ($substrate->bonds) {
    $order1 += $bond->order;
  }
  my $order2 = 0;
  foreach my $bond ($product->bonds) {
    $order2 += $bond->order;
  }
  die "$class substrate and product must coincide on total bond order\n"
    if $order1 != $order2;

  foreach my $atom ($substrate->atoms) {
    die sprintf(
        "$class substrate and product must coincide on atom symbols "
        ."(%s ne %s)\n", $atom->symbol, $mapping->{$atom}->symbol)
      if $atom->symbol ne $mapping->{$atom}->symbol;
  }

  foreach my $atom ($product->atoms) {
    $atom->attr("reaction/mapped", 1);
  }
  foreach my $atom ($substrate->atoms) {
    $mapping->{$atom}->attr("reaction/mapped", 0);
  }
  foreach my $atom ($product->atoms) {
    die "$class atom mapping must be bijective\n"
      if $atom->attr("reaction/mapped");
  }

  # the substrate of the reaction is cloned in order to isolate all
  # changes to bond orders from the given $substrate
  
  my $self = $substrate->clone;
  bless $self, ref $class || $class;

  $self->$_($args{$_}) for (keys %args);
  # the $mapping array gives the product atom which each substrate
  # atom is mapped to, and the %unmapping hash gives back the
  # substrate atom which each product atom is mapped to
  
  my %unmapping;
  foreach my $a1 (keys %$mapping) {
    my $a2 = ${$mapping}{$a1};
    $unmapping{$a2} = $a1;
  }

  # the %bonds hash gives an array of substrate and product bond
  # orders for each pair of atoms which are bonded in substrate or in
  # the product of the reaction

  # first of all, set $bonds{$a1}{$a2} to an array containing only the
  # substrate bond order, for each pair of atoms $a1 and $a2 which are
  # bonded in the substrate

  my %bonds;
  foreach my $bond ($self->bonds) {
    my @atoms = $bond->atoms;
    @atoms[0,1] = @atoms[1,0] unless $atoms[0] lt $atoms[1];
    $bonds{$atoms[0]}{$atoms[1]} = [$bond->order];
  }

  # then, for each pair of atoms $a1 and $a2 which are bonded in the
  # product, append their product bond order to the array
  # $bonds{$a1}{$a2}, preceded by zero (as substrate bond order) if
  # these atoms are not bonded in the substrate
  
  foreach my $bond ($product->bonds) {
    my @atoms = $bond->atoms;
    my $a1 = $unmapping{$atoms[0]};
    my $a2 = $unmapping{$atoms[1]};
    @atoms[0,1] = @atoms[1,0] unless $atoms[0] lt $atoms[1];
    $bonds{$a1}{$a2} = [0] unless defined $bonds{$a1}{$a2};
    push @{$bonds{$a1}{$a2}}, $bond->order;
  }

  # finally, for each pair of atoms $a1 and $a2 which are bonded in
  # the substrate but not in the product, append zero (as product bond
  # order) to the array $bonds{$a1}{$a2}

  foreach my $a1 (keys %bonds) {
    foreach my $a2 (keys %{$bonds{$a1}}) {
      my $a = $bonds{$a1}{$a2};
      push @$a, 0 unless defined $a->[1];
    }
  }

  # now, for each pair of atoms $a1 and $a2 which are bonded in the
  # substrate or in the product, the array array $bonds{$a1}{$a2}
  # contains exactly two entries: the substrate bond order (zero if
  # not bonded) and the product bond order (zero if not bonded)

  # for each bond $bond in the substrate, the substrate bond order is
  # stored in $bond->attr("reaction/before"), and the product bond
  # order (if any) is stored in $bond->attr("reaction/after")
  
  foreach my $bond ($self->bonds) {
    my @atoms = $bond->atoms;
    my @a = @{$bonds{$atoms[0]}{$atoms[1]}};
    $bond->attr("reaction/before", $a[0]);
    $bond->attr("reaction/after", $a[1]);
  }
  
  # further, for each bond $bond in the product but not in the
  # substrate, $bond->attr("reaction/before") is set to zero and the
  # product bond order is stored in $bond->attr("reaction/after")
  
  foreach my $bond ($product->bonds) {
    my @atoms = $bond->atoms;
    my $a1 = $unmapping{$atoms[0]};
    my $a2 = $unmapping{$atoms[1]};
    @atoms[0,1] = @atoms[1,0] unless $atoms[0] lt $atoms[1];
    my @a = @{$bonds{$a1}{$a2}};
    if ($a[0] == 0) {
      my $b = $self->new_bond(atoms =>
			      [$self->by_id($a1), $self->by_id($a2)]);
      $b->attr("reaction/before", $a[0]);
      $b->attr("reaction/after", $a[1]);
    }
  }
  
  return $self;
}

=item $react->substrate

Return a Chemistry::Pattern object that represents the substrate
molecules of the explicit chemical reaction $react.

=cut

# the substrate molecule is obtained from a clone of the reaction by
# breaking all bonds with substrate order equal to zero and setting
# the order of each remaining $bond to $bond->attr("reaction/before")

sub substrate {
  my $react = shift;
  my $self = $react->clone;
  foreach my $bond ($self->bonds) {
    if ($bond->attr("reaction/before") == 0) {
      $bond->delete;
    } else {
      $bond->order($bond->attr("reaction/before"));
      $bond->del_attr("reaction/before");
      $bond->del_attr("reaction/after");
    }
  }
  return $self;
}

=item $react->product

Return a Chemistry::Pattern object that represents the product
molecules of the explicit chemical reaction $react.

=cut

# the product molecule is obtained from a clone of the reaction by
# breaking all bonds with product order equal to zero and setting the
# order of each remaining $bond to $bond->attr("reaction/after")

sub product {
  my $react = shift;
  my $self = $react->clone;
  foreach my $bond ($self->bonds) {
    if ($bond->attr("reaction/after") == 0) {
      $bond->delete;
    } else {
      $bond->order($bond->attr("reaction/after"));
      $bond->del_attr("reaction/before");
      $bond->del_attr("reaction/after");
    }
  }
  return $self;
}

# the map of substrate atoms to product atoms is just the identity
# mapping from $react->substrate->atoms to $react->product->atoms

=item $react->forward($mol, @map)

Forward application of the explicit chemical reaction $react to the
molecule $mol, according to the mapping @map of substrate atoms to
$mol atoms. The substrate of the explicit chemical reaction $react
must be a subgraph of the molecule $mol. Return the modified molecule
$mol.

=cut

sub forward {
  my ($react, $mol, @map) = @_;

  # the %occ hash gives the occurrence of each $react atom in $mol

  my %occ;
  for (my $i = 0; $i < $react->atoms; $i++) {
    $occ{$react->atoms($i+1)} = $map[$i];
  }

  # for each bond $bond in the reaction $react, either change the
  # corresponding bond $b in the molecule $mol to the bond resulting
  # from the forward application of the reaction, or break an existing
  # bond $b, or form a new bond $b between the corresponding atoms $a1
  # and $a2 in $mol

  foreach my $bond ($react->bonds) {
    my @atoms = $bond->atoms;
    my $a1 = $atoms[0];
    my $a2 = $atoms[1];
    my $b; # bond between $occ{$a1} and $occ{$a2} in $mol
    foreach my $bb ($occ{$a1}->bonds) {
      foreach my $aa ($bb->atoms) {	
	if ($aa eq $occ{$a2}) {
	  $b = $bb;
	  last;
	}
      }
      if ($b) { last; }
    }
    if ($b) {
      $b->order($b->order
		-$bond->attr('reaction/before')
		+$bond->attr('reaction/after'));
      if ($b->order == 0) {
	$mol->delete_bond($b);
      }
    } else {
      $mol->new_bond(atoms =>
		     [$mol->by_id($occ{$a1}), $mol->by_id($occ{$a2})],
		    order => $bond->attr('reaction/after'));
    }
  }
}

=item $react->reverse($mol, @map)

Reverse application of the explicit chemical reaction $react to the
molecule $mol, according to the mapping @map of product atoms to $mol
atoms. The product of the explicit chemical reaction $react must be a
subgraph of the molecule $mol. Return the modified molecule $mol.

=cut

sub reverse {
  my ($react, $mol, @map) = @_;

  # the %occ hash gives the occurrence of each $react atom in $mol

  my %occ;
  for (my $i = 0; $i < $react->atoms; $i++) {
    $occ{$react->atoms($i+1)} = $map[$i];
  }

  # for each bond $bond in the reaction $react, either change the
  # corresponding bond $b in the molecule $mol to the bond resulting
  # from the reverse application of the reaction, or break an existing
  # bond $b, or form a new bond $b between the corresponding atoms $a1
  # and $a2 in $mol

  foreach my $bond ($react->bonds) {
    my @atoms = $bond->atoms;
    my $a1 = $atoms[0];
    my $a2 = $atoms[1];
    my $b; # bond between $occ{$a1} and $occ{$a2} in $mol
    foreach my $bb ($occ{$a1}->bonds) {
      foreach my $aa ($bb->atoms) {	
	if ($aa eq $occ{$a2}) {
	  $b = $bb;
	  last;
	}
      }
      if ($b) { last; }
    }
    if ($b) {
      $b->order($b->order
		-$bond->attr('reaction/after')
		+$bond->attr('reaction/before'));
      if ($b->order == 0) {
	$mol->delete_bond($b);
      }
    } else {
      $mol->new_bond(atoms =>
		     [$mol->by_id($occ{$a1}), $mol->by_id($occ{$a2})],
		    order => $bond->attr('reaction/before'));
    }
  }
}

=back

=head1 VERSION

0.02

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Pattern>, L<Chemistry::Tutorial>

Rosselló, F. and G. Valiente, Analysis of metabolic pathways by graph
transformation, in: Proc. 2nd Int. Conf. Graph Transformation, Lecture
Notes in Computer Science 3256 (2004), pp. 73--85.

Rosselló, F. and G. Valiente, Chemical graphs, chemical reaction
graphs, and chemical graph transformation, in: Proc. 2nd Int. Workshop
on Graph-Based Tools, Electronic Notes in Theoretical Computer Science
(2004), in press.

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt> and Gabriel Valiente
E<lt>valiente@lsi.upc.esE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert-Brohman and Gabriel Valiente. All
rights reserved. This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

=cut
