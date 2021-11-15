package PerlMol;

$VERSION = '0.3500';
# $Id: PerlMol.pm,v 1.5 2005/05/13 20:40:33 itubert Exp $

1;

=head1 NAME

PerlMol - Perl modules for molecular chemistry

=head1 SYNOPSIS

    # This is a bundle containing all of the modules
    # of the PerlMol Project and their dependencies.
    # This is not a real module; it is the main index
    # to the documentation of the PerlMol modules.

=head1 DESCRIPTION

PerlMol is a collection of Perl modules for chemoinformatics and computational
chemistry with the philosophy that "simple things should be simple". It should
be possible to write one-liners that use this toolkit to do meaningful
"molecular munging". The PerlMol toolkit provides objects and methods for
representing molecules, atoms, and bonds in Perl; doing substructure matching;
and reading and writing files in various formats.

=head1 DOCUMENTATION

What follows is an index of the relevant documentation.

=head2 Tutorial

=over

=item L<Chemistry::Tutorial> - A good place to start reading the documentation.

=back

=head2 Object oriented modules

The following modules are indented according to the class hierarchy:

=over

=item L<Chemistry::Obj>

=over

=item L<Chemistry::Mol>

=over

=item L<Chemistry::Pattern>

=over

=item L<Chemistry::FormulaPattern>

=item L<Chemistry::MidasPattern>

=back

=item L<Chemistry::Domain>

=item L<Chemistry::MacroMol>

=item L<Chemistry::Ring>
        
=back

=item L<Chemistry::Atom>
        
=over

=item L<Chemistry::Pattern::Atom>
        
=back

=item L<Chemistry::Bond>
        
=over

=item L<Chemistry::Pattern::Bond>

=back
        
=back
    
=item L<Chemistry::File>
    
=over 

=item L<Chemistry::File::Formula>

=item L<Chemistry::File::FormulaPattern>

=item L<Chemistry::File::MDLMol>

=item L<Chemistry::File::MidasPattern>

=item L<Chemistry::File::Mopac>

=item L<Chemistry::File::PDB>

=item L<Chemistry::File::SDF>

=item L<Chemistry::File::SLN>

=item L<Chemistry::File::SMARTS>

=item L<Chemistry::File::SMILES>
    
=item L<Chemistry::File::VRML>

=item L<Chemistry::File::XYZ>
    
=back

=item L<Chemistry::InternalCoords>

=item L<Chemistry::Mok>

=back

=head2 Procedural modules

These are auxiliary modules for which object classes seemed overkill

=over

=item L<Chemistry::3DBuilder>

=item L<Chemistry::Bond::Find>

=item L<Chemistry::InternalCoords::Builder>

=item L<Chemistry::Ring::Find>

=item L<Chemistry::Canonicalize>

=back

=head2 Programs (scripts)

=over

=item L<mok> - an AWK for molecules

=back

=head2 PerlMol bundle description and contents

=over

=item L<PerlMol> - This document

=back

=head2 Publications

=over

=item *

Tubert-Brohman, I. Perl and Chemistry. The Perl Journal 2004-06 (subscription
required) 

=item *

Cozens, S. Molecular Biology With Perl. The Perl Journal 2004, 8[8], 15-19 
(L<http://www.tpj.com/documents/s=7618/tpj0408/>; requires subscription)

=item *

F. Rosselló, G. Valiente. Chemical Graphs, Chemical Reaction Graphs, and
Chemical Graph Transformation. 2nd Int. Workshop on Graph-Based Tools,
Electronic Notes in Computer Science 2005, 127, 157-166. (abstract:
L<http://www.lsi.upc.es/%7Evaliente/abs-grabats-2004.html>; preprint full text:
L<http://tfs.cs.tu-berlin.de/grabats/Final04/valiente.pdf>; published version:
L<http://dx.doi.org/10.1016/j.entcs.2004.12.033>).

=item *

Rosselló, F.; Valiente, G. Graph Transformation in Molecular Biology. (Full
text: L<http://bioinfo.uib.es/~cesc/recerca/he-paper.pdf>).

=back

=head1 EXAMPLES

The "examples" directory in the PerlMol distribution file has several sample
scripts with lots of comments and a few input and output files that show how
one can use PerlMol for common tasks. They can also be browsed online at
L<http://www.perlmol.org/examples/>. Some of the examples are:

    * combinatorial_enumeration
    * file_conversion
    * molgrep
    * pdb_viewer
    * peptide_builder
    * polar_surface_area 

=head1 VERSION INFORMATION

This is the PerlMol bundled release version 0.3500. It includes the following
distributions:

    Chemistry-3DBuilder             0.10
    Chemistry-Bond-Find             0.21
    Chemistry-Canonicalize          0.10
    Chemistry-File-MDLMol           0.20
    Chemistry-File-Mopac            0.15
    Chemistry-File-PDB              0.21
    Chemistry-File-SLN              0.10
    Chemistry-File-SMARTS           0.22
    Chemistry-File-SMILES           0.44
    Chemistry-File-VRML             0.10
    Chemistry-File-XYZ              0.11
    Chemistry-FormulaPattern        0.10
    Chemistry-InternalCoords        0.18
    Chemistry-Isotope               0.11
    Chemistry-MacroMol              0.06
    Chemistry-MidasPattern          0.11
    Chemistry-Mok                   0.25
    Chemistry-Mol                   0.35
    Chemistry-Pattern               0.26
    Chemistry-Reaction              0.02
    Chemistry-Ring                  0.18
    Math-VectorReal                 1.02
    Parse-Yapp                      1.05
    Statistics-Regression           0.15

The version number of a PerlMol bundle is always the same as the version number
of the included Chemistry-Mol distribution, plus two extra digits that
distinguish between different bundles based on the same Chemistry-Mol
distribution. 

=head1 SEE ALSO

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

