from acpype.mol import Angle, Atom, Bond, Dihedral

r = "<[<Atom id=3, name=C1, c3>, <Atom id=4, name=H, hc>, <Atom id=5, name=H, hc>, <Atom id=6, name=H, hc>], ang=0.00>"


def test_atom():
    atom = Atom("C1", "c3", 3, 5, 12.0, 0.5, (1.0, 0.9, 0.5))
    assert str(atom) == "<Atom id=3, name=C1, c3>"


def test_bond():
    a1 = Atom("C1", "c3", 3, 5, 12.0, 0.5, (1.0, 0.9, 0.5))
    a2 = Atom("H", "hc", 4, 5, 1.0, -0.5, (2.0, 1.9, 1.5))
    bond = Bond([a1, a2], 330.0, 1.0969)
    assert str(bond) == repr(bond) == "<[<Atom id=3, name=C1, c3>, <Atom id=4, name=H, hc>], r=1.0969>"


def test_angle():
    a1 = Atom("C1", "c3", 3, 5, 12.0, 0.5, (1.0, 0.9, 0.5))
    a2 = Atom("H", "hc", 4, 5, 1.0, -0.5, (2.0, 1.9, 1.5))
    a3 = Atom("H", "hc", 5, 5, 1.0, 0.1, (0.0, 0.1, -0.5))
    angle = Angle([a1, a2, a3], 46.3, 1.91637234)
    assert (
        str(angle)
        == repr(angle)
        == "<[<Atom id=3, name=C1, c3>, <Atom id=4, name=H, hc>, <Atom id=5, name=H, hc>], ang=109.80>"
    )


def test_dihedral():
    a1 = Atom("C1", "c3", 3, 5, 12.0, 0.5, (1.0, 0.9, 0.5))
    a2 = Atom("H", "hc", 4, 5, 1.0, -0.5, (2.0, 1.9, 1.5))
    a3 = Atom("H", "hc", 5, 5, 1.0, 0.1, (0.0, 0.1, -0.5))
    a4 = Atom("H", "hc", 6, 5, 1.0, -0.1, (0.0, 0.1, 2.0))
    dih = Dihedral([a1, a2, a3, a4], 0.16, 3, 0.0)
    assert str(dih) == repr(dih) == r
