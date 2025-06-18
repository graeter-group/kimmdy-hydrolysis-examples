from math import sqrt
from pathlib import Path
from src.units import nm

import logging

logger = logging.getLogger(__name__)


class GroAtom:
    """An atom in a gromacs gro file

    The fixed format for each atom is:

        - residue number (5 positions, integer)
        - residue name (5 characters)
        - atom name (5 characters)
        - atom number (5 positions, integer)
        - position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
        - velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)

    e.g. C format: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

    NOTE: never trust the atom numbers of a gro file! They wrap around at 99999.

    Attributes
    ----------
    residue_number
        The residue number
    residue_name
        The residue name
    atom_name
        The atom name
    atom_number
        The atom number
    position
        The position in nm
    velocity
        The velocity in nm/ps
    ix
        The index of the atom in the gro file, 0-based, fed during parsing.
    """

    def __init__(
        self,
        residue_number: int,
        residue_name: str,
        atom_name: str,
        atom_number: int,
        position: tuple[float, float, float],
        velocity: tuple[float, float, float] | None = None,
        ix: int = 0,
    ):
        """Create a GroAtom"""
        self.residue_number = residue_number
        self.residue_name = residue_name
        self.atom_name = atom_name
        self.atom_number = atom_number
        self.position = position
        self.velocity = velocity
        self.ix = ix

    def distance(self, other: "GroAtom") -> nm:
        """Calculate the distance between two atoms

        Parameters
        ----------
        other
            The other atom

        Returns
        -------
        nm
            The distance between the two atoms in nm
        """
        if type(other) != GroAtom:
            m = f"Expected GroAtom, got {type(other)}"
            logger.error(m)
            raise TypeError(m)

        s = self.position
        o = other.position

        return nm(sqrt(sum((x - y) ** 2 for x, y in zip(s, o))))

    def angle(self, left: "GroAtom", right: "GroAtom") -> float:
        """Calculate the angle between three atoms

        With self, left and right as positions s, l, r, the angle is calculated of
        l -- s -- r
        i.e. `self` is at the apex of the angle.
        If the atoms are colliding, a warning is logged and 0 is returned.

        Parameters
        ----------
        left
            The first other atom

        right
            The second other atom

        Returns
        -------
        float
            The angle in degrees
        """
        from math import acos, degrees

        s = self.position
        l = left.position
        r = right.position

        if s == l or s == r or l == r:
            m = "Atoms are colliding. Returning 0 degrees."
            logger.warning(m)
            return 0

        ls = [s - l for s, l in zip(s, l)]
        rs = [s - r for s, r in zip(s, r)]

        dot = sum(x * y for x, y in zip(ls, rs))
        norm_ab = sum(x**2 for x in ls) ** 0.5
        norm_cb = sum(x**2 for x in rs) ** 0.5

        return degrees(acos(dot / (norm_ab * norm_cb)))

    def __str__(self):
        # Format each field back to the original format string

        # would overflow the fixed-width fields
        if self.residue_number >= 1e5:
            raise ValueError("residue_number too large")

        if self.atom_number >= 1e5:
            raise ValueError("atom_number too large")

        s = (
            f"{self.residue_number:5d}"
            f"{self.residue_name:<5s}"
            f"{self.atom_name:>5s}"
            f"{self.atom_number:5d}"
            f"{self.position[0]:8.3f}"
            f"{self.position[1]:8.3f}"
            f"{self.position[2]:8.3f}"
        )
        if self.velocity is not None:
            s += (
                f"{self.velocity[0]:8.4f}"
                f"{self.velocity[1]:8.4f}"
                f"{self.velocity[2]:8.4f}"
            )
        return s

    def __repr__(self):
        return f"GroAtom({self.residue_number}, {self.residue_name}, {self.atom_name}, {self.atom_number}, {self.position}, {self.velocity} at {self.ix})"

    @classmethod
    def from_str(cls, line: str, ix: int = 0):
        # lines can also not have all the fields, like
        #   828ACE   HH31    1   3.623   4.154  19.949
        # vs
        # 1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
        l = len(line)
        f1 = int(line[0:5].strip())  # %5d
        f2 = line[5:10].strip()  # %-5s
        f3 = line[10:15].strip()  # %5s
        f4 = int(line[15:20].strip())  # %5d
        f5 = float(line[20:28].strip())  # %8.3f
        f6 = float(line[28:36].strip())  # %8.3f
        f7 = float(line[36:44].strip())  # %8.3f
        residue_number = f1
        residue_name = f2
        atom_name = f3
        atom_number = f4
        position = (f5, f6, f7)
        if l < 52:
            velocity = None
        else:
            f8 = float(line[44:52].strip())  # %8.4f
            f9 = float(line[52:60].strip())  # %8.4f
            f10 = float(line[60:68].strip())  # %8.4f
            velocity = (f8, f9, f10)
        return cls(
            residue_number, residue_name, atom_name, atom_number, position, velocity, ix
        )


class Gro:
    """A gromacs gro file

    Lines contain the following information (top to bottom):
        - title string (free format string, optional time in ps after ‘t=’)
        - number of atoms (free format integer)
        - one line for each atom (fixed format, see GroAtom)
        - box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
          the last 6 values may be omitted (they will be set to zero).
          GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

    For Example:

    ```
    MD of 2 waters, t= 0.0
        6
        1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
        1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
        1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
        2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
        2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
        2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
       1.82060   1.82060   1.82060
    ```
    """

    def __init__(
        self, title: str, atoms: list[GroAtom], box: tuple[float, float, float]
    ):
        self.title = title
        self.n_atoms = len(atoms)
        self.atoms = atoms
        self.box = box

    def __repr__(self):
        return f"Gro({self.title}, with {self.n_atoms} atoms and box {self.box})"

    def __str__(self):
        s = f"{self.title}\n"
        s += f"{self.n_atoms}\n"
        for atom in self.atoms:
            s += f"{atom}\n"
        s += " ".join(map(str, self.box)) + "\n"
        return s


def read_gro(path: str | Path) -> Gro:
    with open(path) as f:
        title = f.readline().strip()
        n_atoms = int(f.readline().strip())
        atoms = [GroAtom.from_str(f.readline(), i) for i in range(n_atoms)]
        b = list(map(float, f.readline().split()))
        box = (b[0], b[1], b[2])

    return Gro(title, atoms, box)


def write_gro(gro: Gro, path: str | Path):
    with open(path, "w") as f:
        f.write(str(gro))
