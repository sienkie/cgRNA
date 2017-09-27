#!/usr/bin/env python

from cabsDock.vector3d import Vector3d
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO, Select, StructureBuilder
from Bio.PDB.Atom import Atom
import glob
import warnings
from cgrna import Lattice

warnings.filterwarnings('error', message='.*discontinuous at.*')

global nucleid_acid_atoms
nucleid_acid_atoms = ['N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'C2', 'C4', 'C5', 'C6', 'C8', 'O2', 'O4', 'O6', 'H1',
                      'H2', 'H3', 'H5', 'H6', 'H8', 'H41', 'H42', 'H21', 'H22', 'H61', 'H62']

global cg_atoms_pur
global cg_atoms_pir

cg_atoms_pur = ["C4'", "P", "N9"]  # A | G
cg_atoms_pir = ["C4'", "P", "N1"]  # C | U


class Atom:
    def __init__(self, v, kind, mass):
        self.coord = v
        self.kind = kind
        self.mass = mass

    def __str__(self, res=None):
        return str(self.kind) + '\t' + str(self.coord) + '\t' + str(res)

    def calculate_dist(self, other):
        self_vec = np.array([float(self.coord.x), float(self.coord.y), float(self.coord.z)])
        other_vec = np.array([float(other.coord.x), float(other.coord.y), float(other.coord.z)])
        return plt.mlab.dist(self_vec, other_vec)


class Residue:
    def __init__(self, kind):
        self.kind = kind
        self.atoms = []
        self.mass_center = None
        self.coords = []

    def add_atom(self, atom):
        self.atoms.append(atom)

    def calculate_mass_center(self):
        mass_sum = 0
        factored_vectors_sum = Vector3d(0, 0, 0)
        n_atom = None
        for atom in self.atoms:
            if atom.kind in nucleid_acid_atoms:
                factored_vectors_sum += float(atom.mass) * atom.coord
                mass_sum += atom.mass
                if atom.kind in ['N1', 'N9']:
                    n_atom = atom
        self.mass_center = factored_vectors_sum / mass_sum
        # print(self.mass_center)
        # TODO creating new atom for mass center now i overwrite N1/N9
        for atom in self.atoms:
            if self.kind in ['  A', '  G'] and atom.kind == 'N9':
                n_atom.coord = self.mass_center
            if self.kind in ['  C', '  U'] and atom.kind == 'N1':
                n_atom.coord = self.mass_center

    def __str__(self):
        s = ''
        for at in self.atoms:
            s += at.__str__(self.kind) + '\n'
        return s

    def set_coords(self):
        for atom in self.atoms:
            if self.kind in ['  A', '  G'] and atom.kind in cg_atoms_pur:
                self.coords.append(atom.coord)
            if self.kind in ['  C', '  U'] and atom.kind in cg_atoms_pir:
                self.coords.append(atom.coord)

    def calculate_p_n_dist(self):
        p = 'P'
        n = 'N9' if self.kind in ['  A', '  G'] else 'N1'
        p_vec = None
        n_vec = None
        for atom in self.atoms:
            if atom.kind == p:
                vector = atom.coord
                p_vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
            elif atom.kind == n:
                vector = atom.coord
                n_vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
        return [(self.kind, plt.mlab.dist(p_vec, n_vec))]

    def get_atom(self, at_name):
        for atom in self.atoms:
            if atom.kind == at_name:
                return atom

    def calculate_c4_n_dist(self):
        p = "C4'"
        n = 'N9' if self.kind in ['  A', '  G'] else 'N1'
        p_vec = None
        n_vec = None
        for atom in self.atoms:
            if atom.kind == p:
                vector = atom.coord
                p_vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
            elif atom.kind == n:
                vector = atom.coord
                n_vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
        return [(self.kind, plt.mlab.dist(p_vec, n_vec))]


class RNAChain:
    def __init__(self, kind):
        self.residues = []
        self.kind = kind
        self.coords = []

    def add_residue(self, res):
        self.residues.append(res)

    def __str__(self):
        s = 'RNA chain ' + str(self.kind) + '\n'
        for res in self.residues:
            s += str(res) + '\n'
        return s

    def calculate_mass_center(self):
        for residue in self.residues:
            # print(residue.kind)
            residue.calculate_mass_center()

    def set_coords(self):
        for res in self.residues:
            res.set_coords()
            for cor in res.coords:
                self.coords.append(cor)

    def calculate_p_n_dist(self):
        inter_dists = []
        intra_dists = []
        for residue in self.residues:
            inter_dists += residue.calculate_p_n_dist()
        for i in range(len(self.residues))[1:]:
            at1 = self.residues[i - 1].get_atom('N9' if self.residues[i - 1].kind in ['  A', '  G'] else 'N1')
            at2 = self.residues[i].get_atom('P')
            intra_dists += [(self.residues[i - 1].kind, at1.calculate_dist(at2))]
        return inter_dists + intra_dists

    def calculate_c4_n_dist(self):
        inter_dists = []
        for residue in self.residues:
            inter_dists += residue.calculate_c4_n_dist()
        return inter_dists


class RNA:
    def __init__(self, pdbfile):
        self.chains = []
        self.cg_atoms = []  # C4' | P | N1 / N9 ==> changed in mass_center
        self.pdb_residues = []
        self.coords = []
        self.changed_coords = False

        p = PDBParser()
        if len(pdbfile.split('/')) > 1:
            self.name = pdbfile.split('/')[-1]
        self.name = self.name.split('.')[0]
        print(pdbfile)

        chains = []
        self.structure = p.get_structure(self.name, pdbfile)
        for chain in self.structure[0]:  # model_1
            ch = RNAChain(chain.id)
            chains.append(ch)
            for residue in chain.get_residues():
                res = Residue(residue.get_resname())
                cg_atoms_list = []
                for atom in residue:
                    coords = list(atom.get_vector())
                    res.add_atom(Atom(Vector3d(coords[0], coords[1], coords[2]), atom.id, atom.mass))
                    if residue.get_resname() in ['  A', '  G'] and atom.id in cg_atoms_pur:
                        cg_atoms_list.append(atom)
                    if residue.get_resname() in ['  C', '  U'] and atom.id in cg_atoms_pir:
                        cg_atoms_list.append(atom)
                if res.atoms != [] and len(cg_atoms_list) == 3:
                    ch.add_residue(res)
                    self.pdb_residues.append(residue)
                    for at in cg_atoms_list:
                        self.cg_atoms.append(at)

        for chain in chains:
            if len(chain.residues) >= 1:
                self.chains.append(chain)

    def __str__(self):
        s = 'RNA ' + self.name + '\n'
        for chain in self.chains:
            s += str(chain) + '\n'
        return s

    def calculate_mass_centers(self):
        for chain in self.chains:
            chain.calculate_mass_center()

    def set_coords(self):
        for chain in self.chains:
            chain.set_coords()
            for cor in chain.coords:
                self.coords.append(cor)

    def create_selected_atoms_file(self, prefix):
        if self.changed_coords:
            cor_list = self.coords
            for i, vector in enumerate(cor_list):
                cor = np.array([float(vector.x), float(vector.y), float(vector.z)])
                self.cg_atoms[i].set_coord(cor)
        self.changed_coords = False
        io = PDBIO()
        io.set_structure(self.structure)
        io.save('structures/' + self.name + '_' + prefix + '_RNA.pdb', SelectCGAtoms(self.cg_atoms))

    def move_nitrogens_to_mass_center(self):
        self.calculate_mass_centers()
        self.set_coords()
        self.changed_coords = True

    def calculate_p_n_dist(self):
        dists = []
        for chain in self.chains:
            dists += chain.calculate_p_n_dist()
        pir = []
        pur = []
        for (res, dist) in dists:
            pur.append(dist) if res in ['  A', '  G'] else pir.append(dist)
        return pur, pir

    def calculate_c4_n_dist(self):
        dists = []
        for chain in self.chains:
            dists += chain.calculate_c4_n_dist()
        pir = []
        pur = []
        for (res, dist) in dists:
            pur.append(dist) if res in ['  A', '  G'] else pir.append(dist)
        return pur, pir


class SelectCGAtoms(Select):
    def __init__(self, atoms):
        self.atoms = atoms

    def accept_atom(self, atom):
        if atom in self.atoms:
            return True
        else:
            return False


def write_csv(filename, colnames, cols):
    assert len(colnames) == len(cols)
    with open(filename, 'w') as f:
        for name in colnames:
            f.write(name + ',')
        f.write('\n')
        max_len = max([len(col) for col in cols])
        for i in range(max_len):
            for col in cols:
                a = ''
                try:
                    a += str(col[i])
                except:
                    pass
                a += ','
                f.write(a)
            f.write('\n')


def get_stats(dir):
    print 'GETTING STATS'
    files = glob.glob(dir + '/*.pdb')
    pur = {'p_mc': [], 'p_n': [], 'c4_n': [], 'c4_mc': []}
    pir = {'p_mc': [], 'p_n': [], 'c4_n': [], 'c4_mc': []}

    for file in files:
        try:
            RNAstruct = RNA(file)
            p_n_pur, p_n_pir = RNAstruct.calculate_p_n_dist()
            pur['p_n'] += p_n_pur
            pir['p_n'] += p_n_pir
            c4_n_pur, c4_n_pir = RNAstruct.calculate_c4_n_dist()
            pur['c4_n'] += c4_n_pur
            pir['c4_n'] += c4_n_pir

            RNAstruct.move_nitrogens_to_mass_center()
            p_mc_pur, p_mc_pir = RNAstruct.calculate_p_n_dist()
            pur['p_mc'] += p_mc_pur
            pir['p_mc'] += p_mc_pir
            c4_mc_pur, c4_mc_pir = RNAstruct.calculate_c4_n_dist()
            pur['c4_mc'] += c4_mc_pur
            pir['c4_mc'] += c4_mc_pir
        except:
            print('Warning was raised as an exception! ' + str(file) + ' excluded from analysis')
    write_csv('Kyu_dist_two_ways.csv',
              ['p_mc_pur_AG', 'p_mc_pir_CU', 'p_n_pur_AG', 'p_n_pir_CU',
               'c4_mc_pur_AG', 'c4_mc_pir_CU', 'c4_n_pur_AG', 'c4_n_pir_CU'],
              [pur['p_mc'], pir['p_mc'], pur['p_n'], pir['p_n'],
               pur['c4_mc'], pir['c4_mc'], pur['c4_n'], pir['c4_n']])


# RNAstruct = RNA('Ding-data-set/1a51.pdb')
# chain = RNAstruct.chains[0]
get_stats('Kyu-data-set')

# cl = Lattice(grid_spacing=0.67, r12=(4.53, 6.45), r13=(4.97, 7.29))
# grid spacing from P-C4' statistics ==> 704 vectors 704:64=11
# r12: P-N1/N9
# r13: P-(C4')-P
#       ==> vectors: 2482
# grid: 0.7 pc_vectors: 560 pn_vectors: 2146
# grid: 0.57 pc_vectors: 1016 pn_vectors: 3940

# cl = Lattice(grid_spacing=0.67, r12=(4.46, 6.4), r13=(4.97, 7.29))
# grid spacing from P-C4' statistics ==> 704 vectors 704:64=11
# r12: P-MC for PURINES A/G
# r13: P-(C4')-P
#       ==> vectors: 2530
# grid: 0.7 pc_vectors: 560 pn_vectors: 2194
# grid: 0.57 pc_vectors: 1016 pn_vectors: 4096

# cl = Lattice(grid_spacing=0.67, r12=(5.08, 7.34), r13=(4.97, 7.29))
# grid spacing from P-C4' statistics ==> 704 vectors 704:64=11
# r12: P-MC for PYRIMIDINES C/U
# r13: P-(C4')-P
#       ==> vectors: 3706
# grid: 0.7 pc_vectors: 560 pn_vectors: 3370
# grid: 0.57 pc_vectors: 1016 pn_vectors: 6100

# cl = Lattice(grid_spacing=0.57, r12=(5.08, 7.34), r13=(4.97, 7.29))
# print len(cl.vectors)
