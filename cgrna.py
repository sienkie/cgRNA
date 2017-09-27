#!/usr/bin/env python

from cabsDock.vector3d import Vector3d
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select


class Atom:
    def __init__(self, v, kind, res):
        self.coord = v
        self.kind = kind
        self.res = res

    def __str__(self):
        return str(self.kind) + '\t' + str(self.coord) + '\t' + str(self.res)


class RNAChain:
    def __init__(self, kind):
        self.atoms = []
        self.kind = kind

    def add_atom(self, atom):
        self.atoms.append(atom)

    def __str__(self):
        s = 'RNA chain ' + str(self.kind) + '\n'
        for at in self.atoms:
            s += str(at) + '\n'
        return s


class RNA:
    def __init__(self, pdbfile):
        self.chains = []
        self.pdb_atoms = []

        mainchain = ["C4'", "P"]

        p = PDBParser()
        if len(pdbfile.split('/')) > 1:
            self.name = pdbfile.split('/')[-1]
        self.name = self.name.split('.')[0]
        print(pdbfile)

        chains = []
        self.structure = p.get_structure(self.name, pdbfile)
        for chain in self.structure[0]:  # model1
            ch = RNAChain(chain.id)
            chains.append(ch)
            for residue in chain.get_residues():
                for atom in residue:
                    if atom.id in mainchain:
                        coords = list(atom.get_vector())
                        ch.add_atom(Atom(Vector3d(coords[0], coords[1], coords[2]), atom.id, residue.get_resname()))
                        self.pdb_atoms.append(atom)

        for chain in chains:
            if len(chain.atoms) >= 1:
                self.chains.append(chain)

    def __str__(self):
        s = 'RNA ' + self.name + '\n'
        for chain in self.chains:
            s += str(chain) + '\n'
        return s


class Lattice:
    """
    This class represent a CABS-like lattice. It is initialized with:
    grid_spacing: distance between grid nodes, default 0.67 A
    r12: tuple with min and max allowed values for P-C4' pseudo-bond length
    r13: tuple with min and max allowed values for P-C4'-P or C4'-P-C4' end distance
    """

    def __init__(self, grid_spacing=0.67, r12=(3.37, 4.34), r13=(4.97, 7.29)):
        self.grid = grid_spacing
        r12min = round((r12[0] / self.grid) ** 2)
        r12max = round((r12[1] / self.grid) ** 2)
        r13min = round((r13[0] / self.grid) ** 2)
        r13max = round((r13[1] / self.grid) ** 2)
        dim = int(r12max ** 0.5)

        self.vectors = []
        for i in range(-dim, dim + 1):
            for j in range(-dim, dim + 1):
                for k in range(-dim, dim + 1):
                    l = i * i + j * j + k * k
                    if r12min <= float(l) <= r12max:
                        self.vectors.append(Vector3d(i, j, k))

        n = len(self.vectors)
        self.good = np.zeros((n, n))
        for i in range(n):
            vi = self.vectors[i]
            for j in range(n):
                vj = self.vectors[j]
                if r13min < (vi + vj).mod2() < r13max and vi.cross(vj).mod2():
                    self.good[i, j] = 1

    def cast(self, ch):
        """
        Function that casts a single protein chain onto the lattice.
        Returns a list of tuples with (x, y, z) coordinates of CA atoms.
        """

        if len(ch.atoms) < 3:
            raise Exception('Protein chain too short!')

        prev = None
        coord = [Vector3d(
            round(ch.atoms[0].coord.x / self.grid),
            round(ch.atoms[0].coord.y / self.grid),
            round(ch.atoms[0].coord.z / self.grid)
        )]

        for atom in ch.atoms[1:]:
            #  iterate over atoms
            min_dr = 1e12
            min_i = -1

            for i, v in enumerate(self.vectors):
                #  iterate over all possible vectors

                # if len(coord) > 2 and self.good[prev, i] == 0:  #??
                #     continue

                if len(coord) < 2 or self.good[prev, i] == 1:

                    new = coord[-1] + v
                    dr = (self.grid * new - atom.coord).mod2()
                    if dr < min_dr:
                        min_dr = dr
                        min_i = i

            if min_i < 0:
                raise Exception('Unsolvable geometric problem!')
            else:
                coord.append(coord[-1] + self.vectors[min_i])
                prev = min_i

        # ??
        coord.insert(0, coord[0] + coord[1] - coord[2])
        coord.append(coord[-1] + coord[-2] - coord[-3])

        return coord


class SelectCGAtoms(Select):
    def __init__(self, atoms):
        self.atoms = atoms

    def accept_atom(self, atom):
        if atom in self.atoms:
            return True
        else:
            return False


def change_coords(RNA, cl, coords):
    for i, vector in enumerate(coords[1:-1]):
        cor = np.array([float(vector.x)*cl.grid, float(vector.y)*cl.grid, float(vector.z)*cl.grid])
        RNA.pdb_atoms[i].set_coord(cor)


def create_main_chain_files(RNA, cl, coords=None):
    io = PDBIO()
    io.set_structure(RNA.structure)
    io.save('structures/' + RNA.name + '_RNA.pdb', SelectCGAtoms(RNA.pdb_atoms))

    if coords:
        change_coords(RNA, cl, coords)
        io = PDBIO()
        io.set_structure(RNA.structure)
        io.save('structures/' + RNA.name + '_cgRNA.pdb', SelectCGAtoms(RNA.pdb_atoms))


RNAstruct = RNA('Ding-data-set/1a51.pdb')
chain = RNAstruct.chains[0]
#lattice = Lattice()
#coords = lattice.cast(chain)
print chain
#print len(coords), len(chain.atoms)

#create_main_chain_files(RNAstruct, lattice, coords)
