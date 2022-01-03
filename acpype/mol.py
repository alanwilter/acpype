"""
    Constructors to define and store the system's topology

    It will create instances for Atoms, AtomTypes, Bonds, Angles and Dihedrals
    where the topology (the relationships between atoms) is defined and
    paramenters are stored.

    Example:

        >>> atom = acpype.mol.Atom(...) # to be improved

    Attributes:
        acpype.mol.Atom     : define Atom
        acpype.mol.AtomType : define AtomType
"""

from typing import List

from acpype.params import Pi


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


class Atom:
    r"""
    Atom Object Definition

    Charges in *prmtop* file are divided by ``18.2223`` to be converted
    in units of the electron charge.

    To convert ``ACOEF`` and ``BCOEF`` to ``r0`` (Å) and ``epsilon`` (ε: kcal/mol), as seen
    in ``gaff.dat`` for example, for a same atom type (``i = j``):

    .. math::
        r_0 &= 1/2 * (2 * A_{coef}/B_{coef})^{1/6} \\
        \epsilon &= 1/(4 * A_{coef}) * B_{coef}^2

    To convert ``r0`` and ``epsilon`` to ``ACOEF`` and ``BCOEF``:

    .. math::
        A_{coef} &= \sqrt{\epsilon_i * \epsilon_j} * (r_{0i} + r_{0j})^{12} \\
        B_{coef} &= 2 * \sqrt{\epsilon_i * \epsilon_j} * (r_{0i} + r_{0j})^6 \\
                          &= 2 * A_{coef}/(r_{0i} + r_{0j})^6

    where index ``i`` and ``j`` for atom types.
    Coordinates are given in Å and masses in Atomic Mass Unit.

    Returns:
        acpype.mol.Atom: atom object
    """

    def __init__(
        self, atomName: str, atomType: AtomType, id_: int, resid: int, mass: float, charge: float, coord: List[float]
    ):
        """
        Args:
            atomName (str): atom name
            atomType (AtomType): atomType object
            id_ (int): atom number index
            resid (int): residues number index
            mass (float): atom mass
            charge (float): atom charge
            coord (List[float]): atom (x,y,z) coordinates
        """
        self.atomName = atomName
        self.atomType = atomType
        self.id = id_
        self.cgnr = id_
        self.resid = resid
        self.mass = mass
        self.charge = charge  # / qConv
        self.coords = coord

    def __str__(self):
        return f"<Atom id={self.id}, name={self.atomName}, {self.atomType}>"

    def __repr__(self):
        return f"<Atom id={self.id}, name={self.atomName}, {self.atomType}>"


class Bond:

    """
    attributes: pair of Atoms, spring constant (kcal/mol), dist. eq. (Ang)
    """

    def __init__(self, atoms, kBond, rEq):
        self.atoms = atoms
        self.kBond = kBond
        self.rEq = rEq

    def __str__(self):
        return f"<{self.atoms}, r={self.rEq}>"

    def __repr__(self):
        return f"<{self.atoms}, r={self.rEq}>"


class Angle:

    """
    attributes: 3 Atoms, spring constant (kcal/mol/rad^2), angle eq. (rad)
    """

    def __init__(self, atoms, kTheta, thetaEq):
        self.atoms = atoms
        self.kTheta = kTheta
        self.thetaEq = thetaEq  # rad, to convert to degree: thetaEq * 180/Pi

    def __str__(self):
        return f"<{self.atoms}, ang={self.thetaEq * 180 / Pi:.2f}>"

    def __repr__(self):
        return f"<{self.atoms}, ang={self.thetaEq * 180 / Pi:.2f}>"


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
        return f"<{self.atoms}, ang={self.phase * 180 / Pi:.2f}>"

    def __repr__(self):
        return f"<{self.atoms}, ang={self.phase * 180 / Pi:.2f}>"
