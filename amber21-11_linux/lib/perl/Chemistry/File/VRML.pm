package Chemistry::File::VRML;

$VERSION = '0.10';

use strict;
use warnings;
use base qw(Chemistry::File);
use Chemistry::Mol;
use POSIX qw(acos);

Chemistry::Mol->register_format(vrml => __PACKAGE__);

my %OPTS = (
    center => 'centerAtoms',
    color  => 'setColor',
    style  => 'setStyle',
    stick_radius => 'setStickRadius',
    ball_radius  => 'setBallRadius',
    compression  => 'setCompression',
);

sub write_mol {
    my ($self, $fh, $mol, %opts) = @_;
    my $vrml = Chemistry::File::VRML::PDB2VRML->new($fh);
    $vrml->add_mol($mol);
    while (my ($key, $val) = each %opts) {
        if (my $method = $OPTS{$key}) {
            $vrml->$method($val);
        }
    }
    $vrml->printVRML;
}

package Chemistry::File::VRML::PDB2VRML;

# Defaults variables
my $Compression = 0;
my $Style = 'Wireframe';    # default display style
my $Color = 'byAtom';       # default color
my $PI    = 3.14159265;

my $RadiusBNS   = 0.2;
my $RadiusStick = 0.15;

# Some global tables

# color indizes per atom type
my %AtomColors = qw(
  '' 5  H  3  HE 5  LI 5  BE 5  B  1  C  4  N  1  O  2  F  7  NE 5  NA 5
  MG 5  AL 5  SI 5  P  6  S  0  CL 7  AR 5  K  5  CA 5  SC 5  TI 5  V  5
  CR 5  MN 5  FE 5  CO 5  NI 5  CU 5  ZN 5  GA 5  GE 5  AS 5  SE 5  BR 7
  KR 5  RB 5  SR 5  Y  5  ZR 5  NB 5  MO 5  TC 5  RU 5  RH 5  PD 5  AG 5
  CD 5  IN 5  SN 5  SB 5  TE 5  I  7  XE 5  CS 5  BA 5  LA 5  CE 5  PR 5
  ND 5  PM 5  SM 5  EU 5  GD 5  TB 5  DY 5  HO 5  ER 5  TM 5  YB 5  LU 5
  HF 5  TA 5  W  5  RE 5  OS 5  IR 5  PT 5  AU 5  HG 5  TL 5  PB 5  BI 5
  PO 5  AT 7  RN 5  FR 5  RA 5  AC 5  TH 5  PA 5  U  5  NP 5  PU 5  AM 5
  CM 5  BK 5  CF 5  XX 5  FM 5  MD 5  NO 5  LW 5
);

# VDW radius per atom type
my %VDWRadius = qw(
  '' 1.00  H  1.08  HE 1.00  LI 1.00  BE 1.00  B  1.00  C  1.54  N  1.48
  O  1.36  F  1.30  NE 1.00  NA 2.30  MG 1.00  AL 2.86  SI 2.10  P  1.75
  S  1.70  CL 1.65  AR 1     K  1     CA 2.75  SC 1     TI 1     V  1
  CR 1     MN 1     FE 2.27  CO 1     NI 1     CU 1.4   ZN 1.4   GA 1
  GE 1     AS 1     SE 1     BR 1.8   KR 1     RB 1     SR 1     Y  1
  ZR 1     NB 1     MO 1     TC 1     RU 1     RH 1     PD 1     AG 1
  CD 1     IN 1     SN 1     SB 1     TE 1     I  1     XE 1     CS 1
  BA 1     LA 1     CE 1     PR 1     ND 1     PM 1     SM 1     EU 1
  GD 1     TB 1     DY 1     HO 1     ER 1     TM 1     YB 1     LU 1
  HF 1     TA 1     W  1     RE 1     OS 1     IR 1     PT 1     AU 1
  HG 1     TL 1     PB 1     BI 1     PO 1     AT 1     RN 1     FR 1
  RA 1     AC 1     TH 1     PA 1     U  1     NP 1     PU 1     AM 1
  CM 1     BK 1     CF 1     XX 1     FM 1     MD 1     NO 1     LW 1
);

#  atom radius per atom type
my %AtomRadius = qw(
  '' 0.00  H  0.37  HE 0.70  LI 1.23  BE 0.89  B  0.80  C  0.77  N  0.74
  O  0.74  F  0.72  NE 0.70  NA 1.57  MG 1.36  AL 1.25  SI 1.17  P  1.10
  S  1.04  CL 0.99  AR 0.70  K  2.03  CA 1.74  SC 1.44  TI 1.32  V  1.22
  CR 1.17  MN 1.16  FE 1.16  CO 1.15  NI 1.17  CU 1.25  ZN 1.25  GA 1.22
  GE 1.21  AS 1.17  SE 0.70  BR 1.24  KR 1.91  RB 1.62  SR 1.45  Y  1.34
  ZR 1.29  NB 1.29  MO 1.24  TC 1.25  RU 1.28  RH 1.34  PD 1.41  AG 1.50
  CD 1.40  IN 1.41  SN 1.37  SB 1.33  TE 0.70  I  1.33  XE 1.98  CS 1.69
  BA 1.69  LA 1.69  CE 1.69  PR 1.69  ND 1.69  PM 1.69  SM 1.69  EU 1.69
  GD 1.69  TB 1.69  DY 1.69  HO 1.69  ER 1.69  TM 1.69  YB 1.69  LU 1.69
  HF 1.44  TA 1.34  W  1.30  RE 1.28  OS 1.26  IR 1.29  PT 1.34  AU 1.44
  HG 1.55  TL 1.54  PB 1.52  BI 1.52  PO 1.40  AT 0.70  RN 2.40  FR 2.00
  RA 1.90  AC 1.90  TH 1.90  PA 1.90  U  1.90  NP 0.70  PU 0.26  AM 1.00
  CM 1.00  BK 1.00  CF 1.00  XX 1.00  FM 1.00  MD 1.00  NO 1.00  LW 1.00
);

#    S       N        O  H
my (@RGBColors) =
  ('1 1 0', '0 0 1', '1 0 0', '1 1 1', '.5 .5 .5', '1 0 1', '1 .5 0', '0 1 0');

#    C       rest     P  Hal

my (%ColorNames) = (
    'yellow' => 0,
    'blue'   => 1,
    'red'    => 2,
    'white'  => 3,
    'grey'   => 4,
    'purple' => 5,
    'brown'  => 6,
    'green'  => 7
);

sub new {
    my ($class, $fh) = @_;
    $class = ref($class) if (ref($class));

    my $this = bless {
        fh         => $fh,
        atoms      => [],
        bonds      => [],
        points     => [],
        lines      => [],
        style      => $Style,
        color      => $Color,
        lineSets   => [],
        lineColors => [],
        indent     => 0,
        DefUse     => {},       # stores shared instances
        BLT        => {},       # bond lookup table for CONECT list
        RadiusStick => $RadiusStick,
        Compression => $Compression,
        RadiusBNS   => $RadiusBNS,
    }, $class;

    return $this;
}

sub add_mol {
    my ($this, $mol) = @_;

    my %atoms;
    for my $atom ($mol->atoms) {
        my $label = uc $atom->symbol;
        my ($x, $y, $z) = $atom->coords->array;
        my $vrml_atom = {
            x      => $x,
            y      => $y,
            z      => $z,
            nr     => scalar(@{$this->{'atoms'}}),
            label  => $label,
            radius => $VDWRadius{$label},
        };
        $atoms{$atom} = $vrml_atom;
        push(@{$this->{'atoms'}}, $vrml_atom);
    }

    for my $bond ($mol->bonds) {
        my ($from, $to) = map { $atoms{$_} } $bond->atoms;
        push(@{$this->{'bonds'}}, {from => $from, to => $to});

    }
    1;
}
###########################################################################
#
# Read PDB file
# Old atoms and bonds won't be deleted, allowing to
# merge several PDB file.
# Syntax: $object->readPDB($fileName);
# Diag: returns undef on error
#
#         1     2       3         4   5     6       7
#1234567890123456789012345678901234567890123456789012345678901234567890123456789
#TOM      5  O5*   A A   1     -16.851  -5.543  74.981  1.00 55.62      3CRO 148
###########################################################################
sub readPDB {
    my $this     = shift;
    my $fileName = shift;

    open(FILE, "$fileName") or return undef;
    while (<FILE>) {
        if (/^ATOM  /) {    # only C,H,O,P,N,S allowed
            my $label = substr($_, 12, 4);
            $label =~ s/\s//g;
            $label = substr(uc $label, 0, 1);    # only the first letter
            my $x = substr($_, 30, 8);
            my $y = substr($_, 38, 8);
            my $z = substr($_, 46, 8);
            my (%atom) = (
                'x'      => $x,
                'y'      => $y,
                'z'      => $z,
                'nr'     => scalar(@{$this->{'atoms'}}),
                'label'  => $label,
                'radius' => $VDWRadius{$label}
            );
            push(@{$this->{'atoms'}}, \%atom);
        } elsif (/^HETATM/) {
            my $label = substr($_, 12, 4);
            $label =~ s/\s//g;
            $label =~ s/[^A-Za-z].*$//g;
            $label = uc $label;
            my $x = substr($_, 30, 8);
            my $y = substr($_, 38, 8);
            my $z = substr($_, 46, 8);
            my (%atom) = (
                'x'      => $x,
                'y'      => $y,
                'z'      => $z,
                'nr'     => scalar(@{$this->{'atoms'}}),
                'label'  => $label,
                'radius' => $VDWRadius{$label}
            );
            push(@{$this->{'atoms'}}, \%atom);
        } elsif (/^CONECT/) {
            my $tmp = substr($_, 0, 69);
            my ($null, $a, @b) = split(/\s+/, $tmp);
            foreach (@b) {
                my $n1 = $a - 1;
                my $n2 = $_ - 1;
                if ($n1 > $n2) { my $n3 = $n1; $n1 = $n2; $n2 = $n3; }
                my $label = $n1 . '_' . $n2;
                next if (exists($this->{'BLT'}->{$label}));
                my $from = $this->{'atoms'}->[$n1];
                my $to   = $this->{'atoms'}->[$n2];
                my (%bond) = ('from' => $from, 'to' => $to);
                push(@{$this->{'bonds'}}, \%bond);
                $this->{'BLT'}->{$label} = 1;
            }
        }
    }
    CORE::close FILE;

    1;
}

###########################################################################
sub setStyle { 
    my ($this, $style) = @_;
    $style = lc $style;
    $style =~ s/[_ ]//g;
    $this->{'style'} = $style; 
}
sub setColor       { my $this = shift; $this->{'color'} = shift; }
sub setStickRadius { my $this = shift; $this->{RadiusStick}  = shift; }
sub setBallRadius  { my $this = shift; $this->{RadiusBNS}    = shift; }
sub setCompression { my $this = shift; $this->{Compression}  = shift; }

###########################################################################
#
# Center all atoms
# Syntax: $object->centerAtoms();
#
###########################################################################
sub centerAtoms {
    my $this = shift;

    my ($cogX, $cogY, $cogZ) = (0, 0, 0);
    foreach (@{$this->{'atoms'}}) {
        $cogX += $_->{'x'};
        $cogY += $_->{'y'};
        $cogZ += $_->{'z'};
    }
    my $numAtoms = @{$this->{'atoms'}};
    $cogX /= $numAtoms;
    $cogY /= $numAtoms;
    $cogZ /= $numAtoms;
    foreach (($cogX, $cogY, $cogZ)) { $_ = sprintf("%.4f", $_); }
    foreach (@{$this->{'atoms'}}) {
        $_->{'x'} -= $cogX;
        $_->{'y'} -= $cogY;
        $_->{'z'} -= $cogZ;
    }
}

###########################################################################
#
# Generate list of lines (cylinders) and points.
# Syntax: $object->_genDisplay();
#
###########################################################################
sub _genDisplay {
    my $this  = shift;
    my $color = $this->{'color'};

    # determine atom colors
    if ($color eq 'byAtom') {
        foreach (@{$this->{'atoms'}}) {
            $_->{'color'} = $AtomColors{$_->{'label'}};
        }
    } else {
        my $c = $ColorNames{$color};
        foreach (@{$this->{'atoms'}}) { $_->{'color'} = $c; }
    }

    # create one point foreach atom
    @{$this->{'points'}} = ();
    foreach (@{$this->{'atoms'}}) {
        my (%point) = (%$_, 'lines' => []);
        push(@{$this->{'points'}}, \%point);
        $_->{'point'} = \%point;
    }

    # create one (or two) lines per bond
    @{$this->{'lines'}} = ();
    foreach (@{$this->{'bonds'}}) {
        my ($at1, $at2) = ($_->{'from'}, $_->{'to'});
        if ($at1->{'color'} == $at2->{'color'}) {
            my $from = $at1->{'point'};
            my (%line) = (
                'from'  => $from,
                'to'    => $at2->{'point'},
                'label' => $at1->{'nr'} . '_' . $at2->{'nr'}
            );
            push(@{$this->{'lines'}}, \%line);
            push(@{$from->{'lines'}}, \%line);
            next;
        }

        # split bond
        my $pt1 = $at1->{'point'};
        my $pt2 = $at2->{'point'};

        my $x = 0.5 * ($at1->{'x'} + $at2->{'x'});
        my $y = 0.5 * ($at1->{'y'} + $at2->{'y'});
        my $z = 0.5 * ($at1->{'z'} + $at2->{'z'});
        my (%point3) = (
            'x'     => $x,
            'y'     => $y,
            'z'     => $z,                            # no label,
            'color' => $at2->{'color'},               # radius or
            'nr'    => scalar(@{$this->{'points'}})
        );                                            # bonds needed
        my (%line1) = (
            'from'  => $pt1,
            'to'    => \%point3,
            'label' => $at1->{'nr'} . '_' . $at2->{'nr'}
        );
        my (%line2) = (
            'from'  => $pt2,
            'to'    => \%point3,
            'label' => $at2->{'nr'} . '_' . $at1->{'nr'}
        );
        push(@{$this->{'lines'}},  \%line1);
        push(@{$this->{'lines'}},  \%line2);
        push(@{$pt1->{'lines'}},   \%line1);
        push(@{$pt2->{'lines'}},   \%line2);
        push(@{$this->{'points'}}, \%point3);
    }
}

###########################################################################
#
# Optimize the wireframe representation.
# Create longer line strips instead of single lines.
#
###########################################################################
sub _genLineSets {
    my $this = shift;

    @{$this->{'lineSets'}}   = ();
    @{$this->{'lineColors'}} = ();
    foreach (@{$this->{'lines'}}) { $_->{'used'} = 0; }
    foreach (@{$this->{'lines'}}) {
        next if $_->{'used'};
        my ($from, $to) = ($_->{'from'}, $_->{'to'});
        push(@{$this->{'lineColors'}}, $from->{'color'});
        my $set = [$from->{'nr'}, $to->{'nr'}];
        push(@{$this->{'lineSets'}}, $set);
        $_->{'used'} = 1;
        my $next  = 1;
        my $bonds = $to->{'lines'};

        while ($next and $bonds) {
            $next = 0;
            my $b;
            foreach $b (@$bonds) {
                next if $b->{'used'};
                my $to = $b->{'to'};
                push(@$set, $to->{'nr'});
                $bonds       = $to->{'lines'};
                $b->{'used'} = 1;
                $next        = 1;
                last;
            }
        }
    }
}

###########################################################################
#
# print VRML SceneGraph
#
###########################################################################
sub printVRML {
    my $this = shift;

    %{$this->{'DefUse'}} = ();
    $this->{'indent'} = 0;
    $this->_printHead();

    $this->_genDisplay();

    $this->_genLineSets(), $this->_printWire()
      if ($this->{'style'} =~ /wire/);
    $this->_printAtoms()
      if ($this->{'style'} =~ /(ball|stick|cpk)/);

    $this->_printTail();
}

###########################################################################
#
# print VRML header
#
###########################################################################
sub _printHead {
    my $this = shift;
    my $fh = $this->{fh};

    print $fh <<EOT;
#VRML V2.0 utf8

Transform {
    children [
        NavigationInfo { type "EXAMINE" }
EOT

    $this->{'indent'} = 2;
}

###########################################################################
sub _printWire {
    my $this = shift;

    $this->_printWireShape();
}

###########################################################################
sub _printWireShape {
    my $this = shift;

    if ($this->{'DefUse'}->{'WireShape'}) {
        $this->_printLine('USE WIRESHAPE');
        return;
    }

    $this->_printLine('DEF WIRESHAPE Shape {');
    $this->{'indent'}++;
    $this->_printWireAppearance();
    $this->_printWireGeometry();
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{'WireShape'} = 1;
}

###########################################################################
sub _printWireAppearance {
    my $this = shift;

    if ($this->{'DefUse'}->{'WireApp'}) {
        $this->_printLine('USE WIREAPP');
        return;
    }

    $this->_printLine("appearance DEF WIREAPP Appearance {");
    $this->{'indent'}++;
    $this->_printLine("material Material { diffuseColor 1 1 1 }");    # dummy
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{'WireApp'} = 1;
}

###########################################################################
sub _printWireGeometry {
    my $this = shift;

    if ($this->{'DefUse'}->{'WireGeo'}) {
        $this->_printLine('geometry USE WIREGEO');
        return;
    }

    $this->_printLine('geometry DEF WIREGEO IndexedLineSet {');
    $this->{'indent'}++;
    $this->_printWireColor();
    $this->_printWireColorIndex();
    $this->_printLine('colorPerVertex FALSE');
    $this->_printWireCoordinate();
    $this->_printWireCoordIndex();
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{'WireGeo'} = 1;
}

###########################################################################
sub _printWireColor {
    my $this = shift;

    if ($this->{'DefUse'}->{'WireCol'}) {
        $this->_printLine('color USE WIRECOL');
        return;
    }

    $this->_printLine('color DEF WIRECOL Color {');
    $this->{'indent'}++;
    $this->_printLine('color [');
    $this->{'indent'}++;
    foreach (@RGBColors) { $this->_printLine("$_,"); }
    $this->{'indent'}--;
    $this->_printLine(']');
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{'WireCol'} = 1;
}

###########################################################################
sub _printWireColorIndex {
    my $this = shift;

    $this->_printLine('colorIndex [');
    $this->{'indent'}++;
    my $lc = $this->{'lineColors'};
    my $i;
    for ($i = 0 ; $i < (@$lc - 8) ; $i += 8) {
        $this->_printLine(join(', ', @$lc[$i .. ($i + 7)]) . ',');
    }
    $this->_printLine(join(', ', @$lc[$i .. $#$lc]) . ',')
      if ($i < @$lc);
    $this->{'indent'}--;
    $this->_printLine(']');
}

###########################################################################
sub _printWireCoordinate {
    my $this = shift;

    if ($this->{'DefUse'}->{'WireCoo'}) {
        $this->_printLine('coord USE WIRECOO');
        return;
    }

    $this->_printLine('coord DEF WIRECOO Coordinate {');
    $this->{'indent'}++;
    $this->_printLine('point [');
    $this->{'indent'}++;
    my ($x, $y, $z);
    foreach (@{$this->{'points'}}) {
        my $x = sprintf("%.4g", $_->{'x'});
        my $y = sprintf("%.4g", $_->{'y'});
        my $z = sprintf("%.4g", $_->{'z'});
        $this->_printLine("$x $y $z,");
    }
    $this->{'indent'}--;
    $this->_printLine(']');
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{'WireCoo'} = 1;
}

###########################################################################
sub _printWireCoordIndex {
    my $this = shift;

    $this->_printLine('coordIndex [');
    $this->{'indent'}++;
    my $ls = $this->{'lineSets'};
    foreach (@$ls) { $this->_printLine(join(', ', @$_) . ', -1,'); }
    $this->{'indent'}--;
    $this->_printLine(']');
}

###########################################################################
sub _printAtoms {
    my $this = shift;

    foreach (@{$this->{'atoms'}}) {
        $this->_printAtom($_) if ($_->{'label'});
    }
}

###########################################################################
sub _printAtom {
    my $this = shift;
    my $atom = shift;

    $this->_printLine("DEF ATOM_$atom->{'nr'} Transform {");
    $this->{'indent'}++;
    $this->_printLine("translation $atom->{'x'} $atom->{'y'} $atom->{'z'}");
    $this->_printLine('children [');
    $this->{'indent'}++;
    $this->_printAtomShape($atom);
    if ($this->{'style'} =~ /stick/) {
        my $point = $atom->{'point'};
        foreach (@{$point->{'lines'}}) { $this->_printBond($_); }
    }
    $this->{'indent'}--;
    $this->_printLine(']');
    $this->{'indent'}--;
    $this->_printLine('}');
}

###########################################################################
sub _printAtomShape {
    my $this = shift;
    my $atom = shift;

    my $l = $atom->{'label'};
    if ($this->{'DefUse'}->{"AtomShape$l"}) {
        $this->_printLine("USE ATOMSHAPE_$l");
        return;
    }

    $this->_printLine("DEF ATOMSHAPE_$l Shape {");
    $this->{'indent'}++;
    $this->_printAtomAppearance($atom);
    $this->_printAtomGeometry($atom);
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{"AtomShape$l"} = 1;
}

###########################################################################
sub _printAtomAppearance {
    my $this = shift;
    my $atom = shift;

    my $c = $atom->{'color'};
    if ($this->{'DefUse'}->{"AtomApp$c"}) {
        $this->_printLine("appearance USE ATOMAPP_$c");
        return;
    }

    $this->_printLine("appearance DEF ATOMAPP_$c Appearance {");
    $this->{'indent'}++;
    $this->_printLine('material Material {');
    $this->{'indent'}++;
    $this->_printLine("diffuseColor $RGBColors[$c]");
    $this->_printLine('specularColor 1 1 1');
    $this->_printLine('shininess 0.75');
    $this->{'indent'}--;
    $this->_printLine('}');
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{"AtomApp$c"} = 1;
}

###########################################################################
sub _printAtomGeometry {
    my $this = shift;
    my $atom = shift;

    my $l = $atom->{'label'};
    if ($this->{'DefUse'}->{"AtomGeo$l"}) {
        $this->_printLine("geometry USE ATOMGEO_$l");
        return;
    }

    my $r = $this->{RadiusStick};
    $r = $this->{RadiusBNS} * $atom->{'radius'} if ($this->{'style'} =~ /ball/);
    $r = $atom->{'radius'} if ($this->{'style'} =~ /cpk/);
    $this->_printLine("geometry DEF ATOMGEO_$l Sphere { radius $r }");

    $this->{'DefUse'}->{"AtomGeo$l"} = 1;
}

###########################################################################
sub _printBond {
    my $this = shift;
    my $bond = shift;

    my ($from, $to) = ($bond->{'from'}, $bond->{'to'});
    my ($tx, $ty, $tz, $s, $ax, $ay, $az, $angle) =
      $this->_calcBond($from, $to);

    foreach ($tx, $ty, $tz, $s, $ax, $ay, $az, $angle) {
        $_ = sprintf("%.5g", $_);
    }

    my $l = $bond->{'label'};
    $this->_printLine("DEF BOND_$l Transform {");
    $this->{'indent'}++;
    $this->_printLine("translation $tx $ty $tz");
    $this->_printLine("scale 1 $s 1");
    $this->_printLine("rotation $ax $ay $az $angle");
    $this->_printLine('children [');
    $this->{'indent'}++;
    $this->_printBondShape($from);
    $this->{'indent'}--;
    $this->_printLine(']');
    $this->{'indent'}--;
    $this->_printLine('}');
}

###########################################################################
sub _printBondShape {
    my $this = shift;
    my $from = shift;

    my $c = $from->{'color'};
    if ($this->{'DefUse'}->{"BondShape$c"}) {
        $this->_printLine("USE BONDSHAPE_$c");
        return;
    }

    $this->_printLine("DEF BONDSHAPE_$c Shape {");
    $this->{'indent'}++;
    $this->_printAtomAppearance($from);   # bond color is the same as atom color
    $this->_printBondGeometry();
    $this->{'indent'}--;
    $this->_printLine('}');

    $this->{'DefUse'}->{"BondShape$c"} = 1;
}

###########################################################################
sub _printBondGeometry {
    my $this = shift;

    if ($this->{'DefUse'}->{'BondGeo'}) {
        $this->_printLine('geometry USE BONDGEO');
        return;
    }

    $this->_printLine(
        "geometry DEF BONDGEO Cylinder { radius $this->{RadiusStick} top FALSE bottom FALSE}"
    );

    $this->{'DefUse'}->{'BondGeo'} = 1;
}

###########################################################################
#
# print VRML tail
#
###########################################################################
sub _printTail {
    my $this = shift;
    my $fh = $this->{fh};
    print $fh <<EOT;
    ]
}
EOT
}

###########################################################################
sub _printLine {
    my $this = shift;
    my $str  = shift;
    my $fh = $this->{fh};

    if ($this->{Compression}) { print $fh "$str\n"; }
    else {
        my $i = "\t" x int($this->{'indent'} >> 1);
        $i .= '    ' if ($this->{'indent'} & 0x1);
        print $fh "$i$str\n";
    }
}

###########################################################################
#
# Calculate bond transformation parameters
# Syntax: @geometry = _calcBond(\%atom1, \%atom2);
#
###########################################################################
sub _calcBond {
    my $this  = shift;
    my $atom1 = shift;
    my $atom2 = shift;

    my ($x1, $y1, $z1) = ($atom1->{'x'}, $atom1->{'y'}, $atom1->{'z'});
    my ($x2, $y2, $z2) = ($atom2->{'x'}, $atom2->{'y'}, $atom2->{'z'});

    my ($dx, $dy, $dz) = ($x2 - $x1, $y2 - $y1, $z2 - $z1);

    # length
    my $s = sqrt($dx * $dx + $dy * $dy + $dz * $dz);

    # translation
    my ($tx, $ty, $tz) = (0.5 * $dx, 0.5 * $dy, 0.5 * $dz);

    ($dx, $dy, $dz) = ($dx / $s, $dy / $s, $dz / $s);

    # rotation axis and angle
    my ($ax, $ay, $az, $angle);
    if    ($dy > 0.9999)  { ($ax, $ay, $az, $angle) = (1, 0, 0, 0); }
    elsif ($dy < -0.9999) { ($ax, $ay, $az, $angle) = (1, 0, 0, $PI); }
    else { ($ax, $ay, $az, $angle) = ($dz, 0, -$dx, acos($dy)); }

    return $tx, $ty, $tz, 0.5 * $s, $ax, $ay, $az, $angle;
}

###########################################################################
#
# Generate connectivities
#
###########################################################################
sub genBonds {
    no warnings 'uninitialized';
    my $this = shift;

    # find largest possible distance
    my $maxR;
    foreach (values %AtomRadius) { $maxR = $_ if ($_ > $maxR); }
    $maxR *= 2.4;

    # find the most negative coordinates to avoid negative indizes
    my ($minX, $minY, $minZ);
    foreach (@{$this->{'atoms'}}) {
        $minX = $_->{'x'} if ($_->{'x'} < $minX);
        $minY = $_->{'y'} if ($_->{'y'} < $minY);
        $minZ = $_->{'z'} if ($_->{'z'} < $minZ);
    }
    $minX -= 2.5 * $maxR;
    $minY -= 2.5 * $maxR;
    $minZ -= 2.5 * $maxR;

    # distribute atoms in a grid with $maxR cell distance
    my (@grid, $maxI, $maxJ, $maxK);
    foreach (@{$this->{'atoms'}}) {
        my $i = int(($_->{'x'} - $minX) / $maxR);
        my $j = int(($_->{'y'} - $minY) / $maxR);
        my $k = int(($_->{'z'} - $minZ) / $maxR);
        push(@{$grid[$i][$j][$k]}, $_);
        $maxI = $i if ($i > $maxI);
        $maxJ = $j if ($j > $maxJ);
        $maxK = $k if ($k > $maxK);
    }

    # loop of grid cells and find bonds
    my ($i, $j, $k, $a, $b, $c);
    for ($i = 1 ; $i <= $maxI ; $i++) {
        for ($j = 1 ; $j <= $maxJ ; $j++) {
            for ($k = 1 ; $k <= $maxK ; $k++) {
                foreach (@{$grid[$i][$j][$k]}) {
                    foreach $a (-1 .. 1) {
                        foreach $b (-1 .. 1) {
                            foreach $c (-1 .. 1) {
                                $this->_atomToGrid($_,
                                    \@{$grid[$i + $a][$j + $b][$k + $c]});
                            }
                        }
                    }
                }
            }
        }
    }
}

###########################################################################
sub _atomToGrid {
    my $this  = shift;
    my $atom1 = shift;
    my $grid  = shift;

    my $n1 = $atom1->{'nr'};
    my ($x1, $y1, $z1) = ($atom1->{'x'}, $atom1->{'y'}, $atom1->{'z'});
    my $ar1 = $AtomRadius{$atom1->{'label'}};

    my $atom2;
    foreach $atom2 (@$grid) {
        my $n2 = $atom2->{'nr'};
        next unless ($n1 < $n2);

        my $ar2 = $ar1 + $AtomRadius{$atom2->{'label'}};
        $ar2 *= 1.2;
        $ar2 *= $ar2;

        my ($x2, $y2, $z2) = ($atom2->{'x'}, $atom2->{'y'}, $atom2->{'z'});
        my ($dx, $dy, $dz) = ($x2 - $x1, $y2 - $y1, $z2 - $z1);
        my $dist = $dx * $dx + $dy * $dy + $dz * $dz;
        next if ($dist > $ar2);
        my $label = $atom1->{'nr'} . '_' . $atom2->{'nr'};
        next if (exists($this->{'BLT'}->{$label}));
        my (%bond) = ('from' => $atom1, 'to' => $atom2);
        push(@{$this->{'bonds'}}, \%bond);
    }
}

###########################################################################

1;

__END__

=head1 NAME

Chemistry::File::VRML - Generate VRML models for molecules

=head1 SYNOPSIS

    use Chemistry::File::PDB;
    use Chemistry::File::VRML;
    use Chemistry::Bond::Find 'find_bonds';

    my $mol = Chemistry::Mol->read('test.pdb');
    find_bonds($mol, orders => 1);
    $mol->write('test.wrl', format => 'vrml', 
        center => 1,
        style  => 'ballAndWire',
        color  => 'byAtom',
    );

=head1 DESCRIPTION

This module generates a VRML (Virtual Reality Modeling Language) representation
of a molecule, which can then be visualized with any VRML viewer. This is a
PerlMol file I/O plugin, and registers the 'vrml' format with
L<Chemistry::Mol>. Note however that this file plugin is write-only; there's no
way of reading a VRML file back into a molecule.

This module is a modification of PDB2VRML by Horst Vollhardt, adapted to the
L<Chemistry::File> interface.

=head2 OPTIONS

The following options may be passed to $mol->write.

=over

=item B<center>

If true, shift the molecules center of geometry into the origin of the
coordinate system. Note: this only affects the output; it does not affect
the coordinates of the atoms in the original Chemistry::Mol object.

=item B<style>

Sets the style for the VRML representation of the molecular structure.
Default is 'Wireframe'. Currently supported styles are:

    Wireframe, BallAndWire,
    Stick, BallAndStick,
    CPK

=item B<color>

Set the overall color of the molecular structure. If the color is
set to 'byAtom', the color the for atoms and bonds is defined by
the atom type. Default is 'byAtom'. Currently supported colors are:

    byAtom,
    yellow, blue, red,
    green, white, brown,
    grey, purple

=item B<stick_radius>

Defines the radius in Angstrom for the cylinders in the 'Stick'
and 'BallAndStick' style. Default is 0.15 .

=item B<ball_radius>

Defines the factor which is multiplied with the VDW radius for
the spheres in the 'BallAndWire' and 'BallAndStick' style. Default
is 0.2 .

=item B<compression>

Turns on/off compression of the output. If turned on, all leading
whitespaces are removed. This produces a less readable but approx.
20% smaller output, the speed is increased by 10% as well.

=back

=head1 AUTHOR

PDB2VRML originally by Horst Vollhardt, horstv@yahoo.com, 1998.
Modified and adapted as Chemistry::File::VRML by Ivan Tubert-Brohman,
itub@cpan.org, 2005.

=head1 COPYRIGHT

PDB2VRML Copyright (c) 1998 by Horst Vollhardt. All rights reserved.
Chemistry::File::VRML modifications Copyright (c) 2005 by Ivan Tubert-Brohman.
All rights reserved.  This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.


=head1 SEE ALSO

PDB2VRML found at 
http://www.realitydiluted.com/mirrors/reality.sgi.com/horstv_basel/pdb2vrml/

PerlMol project at http://www.perlmol.org/

=cut

