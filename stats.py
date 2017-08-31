#!/usr/bin/env python3

import glob, math
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Vector import calc_dihedral, Vector
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import warnings

warnings.filterwarnings('error', message='.*discontinuous at.*')


def get_pdb(ids, dir):
    # Get structures ids
    ids_list = []
    with open(ids, 'r') as file:
        for line in file.readlines():
            ids_list.append(line.split('\n')[0])

    # Selecting structures from PDB
    pdbl = PDBList()

    for i in ids_list:
        print(i)
        pdbl.retrieve_pdb_file(i, pdir=dir, file_format='pdb')


def calculate_angle(atoms):
    a = np.array(atoms[0])
    b = np.array(atoms[1])
    c = np.array(atoms[2])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def bond_length(atoms):
    a = np.array(atoms[0])
    b = np.array(atoms[1])
    return plt.mlab.dist(a, b)


def degrees(rad_angle):
    if rad_angle is None:
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180:
        angle = angle - 360
    while angle < -180:
        angle = angle + 360
    return angle


def get_stats(dir, atoms_list, kind, plot):
    if kind == 'flat':
        return get_flat_angles(dir, atoms_list, plot)
    elif kind == 'torsion':
        return get_torsion_angles(dir, atoms_list, plot)
    elif kind == 'bond':
        points = get_bonds(dir, atoms_list)
        points_reversed = get_bonds(dir, atoms_list[::-1])
        if plot:
            plt.hist(points + points_reversed, color='b', alpha=0.5)
            plt.xlabel(
                'Length of bonds: ' + atoms_list[0] + '-' + atoms_list[1] + ' ' + atoms_list[1] + '-' + atoms_list[0])
            plt.legend()
            plt.show()
        return points + points_reversed


def get_torsion_angles(dir, atoms_list, plot):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    points = []

    for struct in pdbs:
        print(struct)
        structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)
        for chain in structure[0]:  # model1
            i = 0
            coords = []
            for residue in chain:
                for atom in residue:
                    if atom.id == atoms_list[i]:
                        coords.append(atom.get_vector())
                        i = (i + 1) % 4
            coords = coords[:int(len(coords) / 4) * 4]
            coords = [coords[x:x + 100] for x in range(0, len(coords), 4)]
            for four in coords:
                points.append(degrees(calc_dihedral(four[0], four[1], four[2], four[3])))
    if plot:
        plt.hist(points, alpha=0.5)
        plt.xlabel('Torsion angles ' + atoms_list[0] + '-' + atoms_list[1] + '-' + atoms_list[2] + '-' + atoms_list[3])
        plt.show()
    return points


def get_flat_angles(dir, atoms_list, plot):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    points = []

    for struct in pdbs:
        print(struct)
        structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)
        for chain in structure[0]:  # model1
            i = 0
            coords = []
            for residue in chain:
                for atom in residue:
                    if atom.id == atoms_list[i]:
                        coords.append(list(atom.get_vector()))
                        i = (i + 1) % 3
            coords = coords[:int(len(coords) / 3) * 3]
            coords = [coords[x:x + 100] for x in range(0, len(coords), 3)]
            for trio in coords:
                points.append(calculate_angle(trio))
    if plot:
        plt.hist(points)
        plt.xlabel('Angle values for ' + atoms_list[0] + '-' + atoms_list[1] + '-' + atoms_list[2])
        plt.show()
    return points


def get_bonds(dir, atoms_list):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    points = []

    for struct in pdbs:
        print(struct)
        structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)
        for chain in structure[0]:  # model1
            i = 0
            coords = []
            for residue in chain:
                for atom in residue:
                    if atom.get_id() == atoms_list[i]:
                        coords.append(list(atom.get_vector()))
                        i = (i + 1) % 2
            coords = coords[:int(len(coords) / 2) * 2]
            coords = [coords[x:x + 100] for x in range(0, len(coords), 2)]
            for duo in coords:
                points.append(bond_length(duo))
    return points


def dif_plot(datasets, labels, title):
    for i in range(len(datasets)):
        plt.hist(datasets[i], alpha=0.5, label=labels[i])
    plt.title(title)
    plt.legend()
    plt.show()


class NewRes:
    def __init__(self):
        self.p_1 = None
        self.c4 = None
        self.p_2 = None
        self.na = None
        self.norm_sum = None
        self.diff = None
        self.cross = None
        self.A = None

    def set_p_1(self, atom):
        self.p_1 = atom

    def set_p_2(self, atom):
        self.p_2 = atom

    def set_c4(self, atom):
        self.c4 = atom

    def set_na(self, atom):
        self.na = atom

    def local_coords_system(self):
        v1 = self.p_1.get_vector() - self.c4.get_vector()
        v2 = self.p_2.get_vector() - self.c4.get_vector()
        self.norm_sum = v1 + v2
        self.norm_sum = np.array(list(self.norm_sum)) / np.linalg.norm(self.norm_sum)
        self.diff = np.array(list(v1 - v2))
        self.cross = np.cross(self.norm_sum, self.diff)
        self.A = np.matrix([list(self.norm_sum), list(self.diff), list(self.cross)])

    def in_local(self):
        v = np.array(list(self.na.get_vector()))
        if self.A is not None:
            return v * inv(self.A)
        else:
            return None

    def p_dist(self):
        v1 = self.p_1.get_vector()
        v2 = self.p_2.get_vector()
        return plt.mlab.dist(v1, v2)

    def c_p_dist(self):
        v1 = self.p_1.get_vector()
        v2 = self.p_2.get_vector()
        vc = self.c4.get_vector()
        return plt.mlab.dist(v1, vc), plt.mlab.dist(vc, v2)

    def calc_angle(self):
        # P1-C4-P2
        a = calculate_angle([list(self.p_1.get_vector()), list(self.c4.get_vector()), list(self.p_2.get_vector())])
        # P2-C4-N1/N9
        b = calculate_angle([list(self.p_2.get_vector()), list(self.c4.get_vector()), list(self.na.get_vector())])
        # P1-C4-N1/N9
        c = calculate_angle([list(self.p_1.get_vector()), list(self.c4.get_vector()), list(self.na.get_vector())])
        return a, b, c


def get_local_stats(dir):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    atoms_list = ["P", "C4'"]
    nucleotides = ['  A', '  C', '  G', '  U']

    structs_for_stats = []

    for struct in pdbs:  # for every pdb structure
        print(struct)
        selected = []
        structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)

        for res in structure.get_residues():
            if str(res.get_resname()) in nucleotides:
                selected.append(res)

        selected = selected[1:]

        for i in range(len(selected)):
            if i < len(selected) - 1 and selected[i].get_parent() == selected[i + 1].get_parent():
                res1 = selected[i]
                res2 = selected[i + 1]

                new_res = NewRes()

                for at in res1.get_atoms():
                    if at.id == 'P':
                        new_res.set_p_1(at)
                    elif at.id == "C4'":
                        new_res.set_c4(at)
                    elif at.id == 'N9' and res1.get_resname() in ['  A', '  G']:
                        new_res.set_na(at)
                    elif at.id == 'N1' and res1.get_resname() in ['  C', '  U']:
                        new_res.set_na(at)
                for at in res2.get_atoms():
                    if at.id == 'P':
                        new_res.set_p_2(at)

                if new_res.p_1 != None and new_res.p_2 != None and new_res.c4 != None and new_res.na != None:
                    structs_for_stats.append(new_res)


    # print(len(structs_for_stats))
    # for i in range(len(structs_for_stats)):
    #     s = structs_for_stats[i]
    #     print(i, s.p_1.id, s.c4.id, s.na.id, s.p_2.id)

    return structs_for_stats

class NewRes2(NewRes):
    def __init__(self):
        NewRes.__init__(self)
        self.c4_before = None
        self.c4_after = None

    def set_c4_before(self, atom):
        self.c4_before = atom

    def set_c4_after(self, atom):
        self.c4_after = atom

    def c_dist(self):
        v1 = self.c4.get_vector()
        v2 = self.c4_after.get_vector()
        v0 = self.c4_before.get_vector()
        return plt.mlab.dist(v1, v2) + plt.mlab.dist(v0, v1)

    def cc_dist(self):
        v1 = self.c4.get_vector()
        v0 = self.c4_before.get_vector()
        return plt.mlab.dist(v0, v1)

    def calc_ccc_angle(self):
        return calculate_angle([list(self.c4_before.get_vector()), list(self.c4.get_vector()), list(self.c4_after.get_vector())])

def check_assumption(dir):
    o = open(dir.split('-')[0] + '_stats.csv', 'w')
    # (x,y,z) -> N1 for A & G / N9 for C & U
    o.write('x,y,z,distance_P-P,angle_P-C4-P,angle_P-C4-N1/9_down,angle_P-C4-N1/9_up' + '\n')
    structs_for_stats = get_local_stats(dir)
    for s in structs_for_stats:
        s.local_coords_system()
        cor = s.in_local()
        dist = s.p_dist()
        a, b, c = s.calc_angle()
        o.write(
            str(cor[0, 0]) + ',' + str(cor[0, 1]) + ',' + str(cor[0, 2]) + ',' + str(dist) + ',' + str(a) + ',' + str(
                b) + ',' + str(c) + '\n')
    o.close()

def get_local_stats_c_dist(dir):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    atoms_list = ["P", "C4'"]
    nucleotides = ['  A', '  C', '  G', '  U']

    structs_for_stats = []

    for struct in pdbs:  # for every pdb structure
        try:
            print(struct)
            selected = []
            structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)

            for res in structure.get_residues():
                if str(res.get_resname()) in nucleotides:
                    selected.append(res)

            selected = selected[1:]

            for i in range(len(selected))[1:]:
                if i < len(selected) - 1 and selected[i].get_parent() == selected[i + 1].get_parent() and selected[i].get_parent() == selected[i - 1].get_parent():
                    res0 = selected[i-1]
                    res1 = selected[i]
                    res2 = selected[i + 1]

                    new_res = NewRes2()

                    for at in res1.get_atoms():
                        if at.id == 'P':
                            new_res.set_p_1(at)
                        elif at.id == "C4'":
                            new_res.set_c4(at)
                        elif at.id == 'N9' and res1.get_resname() in ['  A', '  G']:
                            new_res.set_na(at)
                        elif at.id == 'N1' and res1.get_resname() in ['  C', '  U']:
                            new_res.set_na(at)
                    for at in res2.get_atoms():
                        if at.id == 'P':
                            new_res.set_p_2(at)
                        elif at.id == "C4'":
                            new_res.set_c4_after(at)
                    for at in res0.get_atoms():
                        if at.id == "C4'":
                            new_res.set_c4_before(at)

                    if new_res.p_1 != None and new_res.p_2 != None and new_res.c4 != None and new_res.na != None and new_res.c4_after != None and new_res.c4_before != None:
                        structs_for_stats.append(new_res)
        except:
            print('Warning was raised as an exception!')

    # print(len(structs_for_stats))
    # for i in range(len(structs_for_stats)):
    #     s = structs_for_stats[i]
    #     print(i, s.p_1.id, s.c4.id, s.na.id, s.p_2.id)

    return structs_for_stats

def get_local_stats_c_dist(dir):
    p = PDBParser()
    pdbs = glob.glob(dir + '/*', recursive=True)
    atoms_list = ["P", "C4'"]
    nucleotides = ['  A', '  C', '  G', '  U']

    structs_for_stats = []

    for struct in pdbs:  # for every pdb structure
        try:
            print(struct)
            selected = []
            structure = p.get_structure(struct.split('/')[1].split('.')[0], struct)

            for res in structure.get_residues():
                if str(res.get_resname()) in nucleotides:
                    selected.append(res)

            selected = selected[1:]

            for i in range(len(selected))[1:]:
                if i < len(selected) - 1 and selected[i].get_parent() == selected[i + 1].get_parent() and selected[i].get_parent() == selected[i - 1].get_parent():
                    res0 = selected[i-1]
                    res1 = selected[i]
                    res2 = selected[i + 1]

                    new_res = NewRes2()

                    for at in res1.get_atoms():
                        if at.id == 'P':
                            new_res.set_p_1(at)
                        elif at.id == "C4'":
                            new_res.set_c4(at)
                        elif at.id == 'N9' and res1.get_resname() in ['  A', '  G']:
                            new_res.set_na(at)
                        elif at.id == 'N1' and res1.get_resname() in ['  C', '  U']:
                            new_res.set_na(at)
                    for at in res2.get_atoms():
                        if at.id == 'P':
                            new_res.set_p_2(at)
                        elif at.id == "C4'":
                            new_res.set_c4_after(at)
                    for at in res0.get_atoms():
                        if at.id == "C4'":
                            new_res.set_c4_before(at)

                    if new_res.p_1 != None and new_res.p_2 != None and new_res.c4 != None and new_res.na != None and new_res.c4_after != None and new_res.c4_before != None:
                        structs_for_stats.append(new_res)
        except:
            print('Warning was raised as an exception!')

    # print(len(structs_for_stats))
    # for i in range(len(structs_for_stats)):
    #     s = structs_for_stats[i]
    #     print(i, s.p_1.id, s.c4.id, s.na.id, s.p_2.id)

    return structs_for_stats

def check_assumption_c_dist(dir):
    o = open(dir.split('-')[0] + '_ccc_stats.csv', 'w')
    # (x,y,z) -> N1 for A & G / N9 for C & U
    o.write('x,y,z,distance_C4-C4-C4, C4-C4-C4_angle' + '\n')
    structs_for_stats = get_local_stats(dir)
    for s in structs_for_stats:
        s.local_coords_system()
        cor = s.in_local()
        dist = s.c_dist()
        angle = s.calc_ccc_angle()
        o.write(
            str(cor[0, 0]) + ',' + str(cor[0, 1]) + ',' + str(cor[0, 2]) + ',' + str(dist) + ',' + str(angle) + '\n')
    o.close()

def dist_stats(dir):
    o = open(dir.split('-')[0] + '_dists.csv', 'w')
    o.write('distPP, distCP, distPC, distCC' + '\n')
    structs_for_stats = get_local_stats_c_dist(dir)
    for s in structs_for_stats:
        cp, pc = s.c_p_dist()
        pp = s.p_dist()
        cc = s.cc_dist()
        o.write(str(pp) + ',' + str(pc) + ',' + str(cp) + ',' + str(cc) + '\n')
    o.close()

# check_assumption('Kyu-data-set')
# check_assumption_c_dist('Kyu-data-set')

dist_stats('Ding-data-set')

# get_pdb('Ding-data-set.txt', 'Ding-data-set')
# get_stats('Ding-data-set', ["C4'", "P", "C4'"], 'flat')
# get_stats('Ding-data-set', ["P", "C4'", "P", "C4'"], 'torsion')  # theta
# get_stats('Ding-data-set', ["C4'", "P", "C4'", "P"], 'torsion')  # eta

# kyu = get_stats('Kyu-data-set', ["C4'", "P", "C4'", "P"], 'torsion', plot=False)
# ding = get_stats('Ding-data-set', ["C4'", "P", "C4'", "P"], 'torsion', plot=False)
# dif_plot([ding, kyu], ['Ding data set', 'Kyungsook data set'], "RNA torsion eta angles (C4'-P-C4'-P)")
