#!/usr/bin/env python3

from CABS.cabsDock.vector3d import Vector3d
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO, Select, StructureBuilder
from Bio.PDB.Atom import Atom
import glob, rmsd
from random import sample
from statistics import median
import warnings
from cgrna import Lattice

warnings.filterwarnings('error', message='.*discontinuous at.*')

global nucleid_acid_atoms
nucleid_acid_atoms = ['N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'C2', 'C4', 'C5', 'C6', 'C8', 'O2', 'O4', 'O6', 'H1',
                      'H2', 'H3', 'H5', 'H6', 'H8', 'H41', 'H42', 'H21', 'H22', 'H61', 'H62']

global cg_atoms_pur
global cg_atoms_pir
global cg_tiads

cg_atoms_pur = ["C4'", "P", "N9"]  # A | G
cg_atoms_pir = ["C4'", "P", "N1"]  # C | U

cg_tiads = {
    '  A': ["C4'", "P", "N6", "C8", "C2"],
    '  G': ["C4'", "P", "N2", "C8", "O6"],
    '  C': ["C4'", "P", "N4", "O6", "C6"],
    '  U': ["C4'", "P", "O4", "O2", "C6"]
}


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
        # TODO creating diff_stats atom for mass center now i overwrite N1/N9
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

    def calculate_c4_cu_cg_dist(self):
        # C4-CU lub C4-CG ==> A,G:CG=C8 | C,U:CU=C6
        at1 = "C4'"
        at2 = "C8" if self.kind in ['  A', '  G'] else 'C6'
        at1vec = None
        at2vec = None
        for atom in self.atoms:
            if atom.kind == at1:
                vector = atom.coord
                at1vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
            elif atom.kind == at2:
                vector = atom.coord
                at2vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
        return [(self.kind, plt.mlab.dist(at1vec, at2vec))]

    def calculate_p_cu_cg_dist(self):
        # P-CU lub P-CG ==> A,G:CG=C8 | C,U:CU=C6
        at1 = "P"
        at2 = "C8" if self.kind in ['  A', '  G'] else 'C6'
        at1vec = None
        at2vec = None
        for atom in self.atoms:
            if atom.kind == at1:
                vector = atom.coord
                at1vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
            elif atom.kind == at2:
                vector = atom.coord
                at2vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
        return [(self.kind, plt.mlab.dist(at1vec, at2vec))]

    def get_norm_vector(self, at1, at2):
        at1vec = None
        at2vec = None
        for atom in self.atoms:
            if atom.kind == at1:
                vector = atom.coord
                at1vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
            elif atom.kind == at2:
                vector = atom.coord
                at2vec = np.array([float(vector.x), float(vector.y), float(vector.z)])
        return (at1vec - at2vec)/float(plt.mlab.dist(at1vec, at2vec))



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

    def calculate_c4_c4_dist(self):
        intra_dists = []
        for i in range(len(self.residues))[1:]:
            at1 = self.residues[i - 1].get_atom("C4'")
            at2 = self.residues[i].get_atom("C4'")
            intra_dists += [at1.calculate_dist(at2)]
        return intra_dists

    def calculate_c4_n_dist(self):
        inter_dists = []
        for residue in self.residues:
            inter_dists += residue.calculate_c4_n_dist()
        return inter_dists

    def calculate_c4_cu_cg_dist(self):
        # C4-CU lub C4-CG ==> A,G:CG=C8 | C,U:CU=C6
        inter_dists = []
        for residue in self.residues:
            inter_dists += residue.calculate_c4_cu_cg_dist()
        return inter_dists

    def calculate_p_cu_cg_dist(self):
        # P-CU lub P-CG ==> A,G:CG=C8 | C,U:CU=C6
        inter_dists = []
        for residue in self.residues:
            inter_dists += residue.calculate_p_cu_cg_dist()
        return inter_dists

    def get_triangles(self, base=None):
        # iterate through residue:
        # - get P and N from "n-1" residue and P from "n" residue
        # - next P and N from "n" residue and P from "n+1" residue
        # and so on ...
        trios = []
        for i in range(len(self.residues))[1:]:
            if base == None or self.residues[i - 1].kind == base:
                atp1 = self.residues[i - 1].get_atom('P')
                atn = self.residues[i - 1].get_atom('N9' if self.residues[i - 1].kind in ['  A', '  G'] else 'N1')
                atp2 = self.residues[i].get_atom('P')
                trios.append([atp1, atn, atp2])
        return trios


class RNA:
    def __init__(self, pdbfile):
        self.chains = []
        self.cg_atoms = []  # C4' | P | N1 / N9 ==> changed in mass_center
        self.atoms_triads = []  # according to cg_triads dictionary
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
                atoms_pdb = []
                for atom in residue:
                    coords = list(atom.get_vector())
                    res.add_atom(Atom(Vector3d(coords[0], coords[1], coords[2]), atom.id, atom.mass))
                    if residue.get_resname() in ['  A', '  G'] and atom.id in cg_atoms_pur:
                        cg_atoms_list.append(atom)
                    if residue.get_resname() in ['  C', '  U'] and atom.id in cg_atoms_pir:
                        cg_atoms_list.append(atom)
                    if residue.get_resname() in cg_tiads.keys() and atom.id in cg_tiads[residue.get_resname()]:
                        atoms_pdb.append(atom)
                if res.atoms != [] and len(cg_atoms_list) == 3 and len(atoms_pdb) == 5:
                    ch.add_residue(res)
                    for at in atoms_pdb:
                        self.atoms_triads.append(at)
                    self.pdb_residues.append(residue)
                    for at in cg_atoms_list:
                        self.cg_atoms.append(at)

        for chain in chains:
            if len(chain.residues) >= 1:
                self.chains.append(chain)

        self.atoms = [] #MAIN CHAIN
        for atom in self.get_atoms():
            if atom.kind in ["C4'","P"]:
                self.atoms.append(atom)
        self.cg_atoms_all = self.get_atoms()

    def __str__(self):
        s = 'RNA ' + self.name + '\n'
        for chain in self.chains:
            s += str(chain) + '\n'
        return s

    def get_atoms(self):
        return [
            Atom(Vector3d(list(at.get_vector())[0], list(at.get_vector())[1], list(at.get_vector())[2]), at.id, at.mass)
            for at in self.cg_atoms]

    def calculate_mass_centers(self):
        for chain in self.chains:
            chain.calculate_mass_center()

    def set_coords(self):
        for chain in self.chains:
            chain.set_coords()
            for cor in chain.coords:
                self.coords.append(cor)

    def create_selected_atoms_file(self, prefix, atoms=None):
        if not atoms:
            atoms = self.cg_atoms
        if self.changed_coords:
            cor_list = self.coords
            for i, vector in enumerate(cor_list):
                cor = np.array([float(vector.x), float(vector.y), float(vector.z)])
                self.cg_atoms[i].set_coord(cor)
        self.changed_coords = False
        io = PDBIO()
        io.set_structure(self.structure)
        io.save('structures/' + self.name + '_' + prefix + '_RNA.pdb', SelectCGAtoms(atoms))

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

    def calculate_c4_cu_cg_dist(self):
        # C4-CU lub C4-CG ==> A,G:CG=C8 | C,U:CU=C6
        dists = []
        for chain in self.chains:
            dists += chain.calculate_c4_cu_cg_dist()
        pir = []
        pur = []
        for (res, dist) in dists:
            pur.append(dist) if res in ['  A', '  G'] else pir.append(dist)
        return pur, pir

    def calculate_p_cu_cg_dist(self):
        # P-CU lub P-CG ==> A,G:CG=C8 | C,U:CU=C6
        dists = []
        for chain in self.chains:
            dists += chain.calculate_p_cu_cg_dist()
        pir = []
        pur = []
        for (res, dist) in dists:
            pur.append(dist) if res in ['  A', '  G'] else pir.append(dist)
        return pur, pir

    def calculate_c4_c4_dist(self):
        dists = []
        for chain in self.chains:
            dists += chain.calculate_c4_c4_dist()
        return dists

    def get_triangles(self, base=None):
        trios = []
        for chain in self.chains:
            trios += chain.get_triangles(base)
        return trios


class SelectCGAtoms(Select):
    def __init__(self, atoms):
        self.atoms = atoms

    def accept_atom(self, atom):
        if atom in self.atoms:
            return True
        else:
            return False
