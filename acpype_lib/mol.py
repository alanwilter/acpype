from acpype_lib.params import Pi


class Atom:

    """
    Charges in prmtop file has to be divide by 18.2223 to convert to charge
    in units of the electron charge.
    To convert ACOEF and BCOEF to r0 (Ang.) and epsilon (kcal/mol), as seen
    in gaff.dat for example; same atom type (i = j):
        r0 = 1/2 * (2 * ACOEF/BCOEF)^(1/6)
        epsilon = 1/(4 * A) * BCOEF^2
    To convert r0 and epsilon to ACOEF and BCOEF
        ACOEF = sqrt(ep_i * ep_j) * (r0_i + r0_j)^12
        BCOEF = 2 * sqrt(ep_i * ep_j) * (r0_i + r0_j)^6
              = 2 * ACOEF/(r0_i + r0_j)^6
    where index i and j for atom types.
    Coord is given in Ang. and mass in Atomic Mass Unit.
    """

    def __init__(self, atomName, atomType, id_, resid, mass, charge, coord):
        self.atomName = atomName
        self.atomType = atomType
        self.id = id_
        self.cgnr = id_
        self.resid = resid
        self.mass = mass
        self.charge = charge  # / qConv
        self.coords = coord

    def __str__(self):
        return "<Atom id=%s, name=%s, %s>" % (self.id, self.atomName, self.atomType)

    def __repr__(self):
        return "<Atom id=%s, name=%s, %s>" % (self.id, self.atomName, self.atomType)


class AtomType:

    """
    AtomType per atom in gaff or amber.
    """

    def __init__(self, atomTypeName, mass, ACOEF, BCOEF):
        self.atomTypeName = atomTypeName
        self.mass = mass
        self.ACOEF = ACOEF
        self.BCOEF = BCOEF

    def __str__(self):
        return "<AtomType=%s>" % self.atomTypeName

    def __repr__(self):
        return "<AtomType=%s>" % self.atomTypeName


class Bond:

    """
    attributes: pair of Atoms, spring constant (kcal/mol), dist. eq. (Ang)
    """

    def __init__(self, atoms, kBond, rEq):
        self.atoms = atoms
        self.kBond = kBond
        self.rEq = rEq

    def __str__(self):
        return "<%s, r=%s>" % (self.atoms, self.rEq)

    def __repr__(self):
        return "<%s, r=%s>" % (self.atoms, self.rEq)


class Angle:

    """
    attributes: 3 Atoms, spring constant (kcal/mol/rad^2), angle eq. (rad)
    """

    def __init__(self, atoms, kTheta, thetaEq):
        self.atoms = atoms
        self.kTheta = kTheta
        self.thetaEq = thetaEq  # rad, to convert to degree: thetaEq * 180/Pi

    def __str__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.thetaEq * 180 / Pi)

    def __repr__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.thetaEq * 180 / Pi)


class Dihedral:

    """
    attributes: 4 Atoms, spring constant (kcal/mol), periodicity,
    phase (rad)
    """

    def __init__(self, atoms, kPhi, period, phase):
        self.atoms = atoms
        self.kPhi = kPhi
        self.period = period
        self.phase = phase  # rad, to convert to degree: kPhi * 180/Pi

    def __str__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.phase * 180 / Pi)

    def __repr__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.phase * 180 / Pi)
