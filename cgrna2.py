#!/usr/bin/env python3

from modules import *

warnings.filterwarnings('error', message='.*discontinuous at.*')


def count_rmsd(t1, t2):
    P = np.array([[float(at.coord.x), float(at.coord.y), float(at.coord.z)] for at in t1])
    Q = np.array([[float(at.coord.x), float(at.coord.y), float(at.coord.z)] for at in t2])

    # print("RMSD before translation: ", rmsd.kabsch_rmsd(P, Q))
    P -= rmsd.centroid(P)
    Q -= rmsd.centroid(Q)
    # print("RMSD after translation: ", rmsd.kabsch_rmsd(P, Q))
    return rmsd.kabsch_rmsd(P, Q)


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
    print('GETTING STATS')
    files = glob.glob(dir + '/*.pdb')
    pur = {'p_mc': [], 'p_n': [], 'c4_n': [], 'c4_mc': [], 'c4_cg': [], 'p_cg': []}
    pir = {'p_mc': [], 'p_n': [], 'c4_n': [], 'c4_mc': [], 'c4_cu': [], 'p_cu': []}
    c4_c4 = []

    for file in files:
        try:
            RNAstruct = RNA(file)
            p_n_pur, p_n_pir = RNAstruct.calculate_p_n_dist()
            pur['p_n'] += p_n_pur
            pir['p_n'] += p_n_pir
            c4_n_pur, c4_n_pir = RNAstruct.calculate_c4_n_dist()
            pur['c4_n'] += c4_n_pur
            pir['c4_n'] += c4_n_pir
            c4_cg, c4_cu = RNAstruct.calculate_c4_cu_cg_dist()
            pur['c4_cg'] += c4_cg
            pir['c4_cu'] += c4_cu
            p_cg, p_cu = RNAstruct.calculate_p_cu_cg_dist()
            pur['p_cg'] += p_cg
            pir['p_cu'] += p_cu
            c4_c4 += RNAstruct.calculate_c4_c4_dist()

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
              ['p_mc_pur_AG', 'p_mc_pir_CU', 'p_n_pur_AG', 'p_n_pir_CU', 'c4_mc_pur_AG', 'c4_mc_pir_CU',
               'c4_n_pur_AG', 'c4_n_pir_CU', 'c4_cg_pur', 'c4_cu_pir', 'p_cg_pur', 'p_cu_pir', 'c4_c4'],
              [pur['p_mc'], pir['p_mc'], pur['p_n'], pir['p_n'], pur['c4_mc'], pir['c4_mc'],
               pur['c4_n'], pir['c4_n'], pur['c4_cg'], pir['c4_cu'], pur['p_cg'], pir['p_cu'], c4_c4])


def triangle_stats(dirs, random_from_struct=1):
    print('GETTING STATS')
    files = []
    for dir in dirs:
        files += glob.glob(dir + '/*.pdb')
    mean_rmsd_inter = []
    rmsd_intra = []
    triangles_whole_set = []
    excluded_files = 0
    for file in files:
        triangles = []
        rmsd_inter = []
        try:
            RNAstruct = RNA(file)
            triangles += RNAstruct.get_triangles()
            for n in range(len(triangles))[1:]:
                rmsd_inter.append(count_rmsd(triangles[n - 1], triangles[n]))
            mean_rmsd_inter.append(sum(rmsd_inter) / len(rmsd_inter))
            triangles_whole_set += sample(triangles, random_from_struct)
        except:
            excluded_files += 1
            print('Warning was raised as an exception! ' + str(file) + ' excluded from analysis')
    print('Excluded files: ', excluded_files)
    for n in range(len(triangles_whole_set)):
        for m in range(len(triangles_whole_set)):
            if n != m:
                rmsd_intra.append(count_rmsd(triangles_whole_set[n], triangles_whole_set[m]))

    print(len(triangles_whole_set), len(mean_rmsd_inter), len(rmsd_intra))

    print('Mean RMSD from random triangles: ', sum(rmsd_intra) / len(rmsd_intra), 'Median: ', median(rmsd_intra))

    write_csv('RMSDs_all.csv', ['mean_inter_structure_rmsd', 'intra_structure_rmsd'], [mean_rmsd_inter, rmsd_intra])


def get_triangles(dirs, output, base=None):
    files = []
    excluded_files = 0
    with open(output, 'w') as f:
        for dir in dirs:
            files += glob.glob(dir + '/*.pdb')
            for file in files:
                try:
                    RNAstruct = RNA(file)
                    triangles = RNAstruct.get_triangles(base)
                    for tri in triangles:
                        triangle = ''
                        for at in tri:
                            triangle += str(at.coord.x) + ',' + str(at.coord.y) + ',' + str(at.coord.z) + ','
                        f.write(triangle + '\n')
                except:
                    excluded_files += 1
                    print('Warning was raised as an exception! ' + str(file) + ' excluded from analysis')
            print('Excluded files: ', excluded_files)


def lattice_and_scalars(dirs):
    files = []
    for dir in dirs:
        files += glob.glob(dir + '/*.pdb')

    cross_products = []
    for file in files:
        try:
            RNAstruct = RNA(file)

            cl = Lattice(grid_spacing=0.67)
            cn_vectors = []
            # every possible vector from every C4' on lattice
            for chain in RNAstruct.chains:
                for residue in chain.residues:
                    cn_vectors.append(residue.get_norm_vector("C4'", 'N9' if residue.kind in ['  A', '  G'] else 'N1'))

            min_val = 10000
            for cn in cn_vectors:
                for vector in cl.vectors:
                    v = np.array([float(vector.x), float(vector.y), float(vector.z)]) / float(
                        plt.mlab.dist(np.array([float(vector.x), float(vector.y), float(vector.z)]),
                                      np.array([float(0), float(0), float(0)])))
                    dist = plt.mlab.dist(np.cross(cn, v), np.array([float(0), float(0), float(0)]))
                    if dist < min_val:
                        min_val = dist
                cross_products.append(min_val)
        except:
            print('Warning was raised as an exception! ' + str(file) + ' excluded from analysis')
    write_csv('cross_products_all.csv', ['cross_products'], [cross_products])


# triangle_stats(['Ding-data-set','Kyu-data-set'], 1)

# get_triangles(['Ding-data-set', 'Kyu-data-set'], "similarity/triangles.csv")

# lattice_and_scalars(['Ding-data-set', 'Kyu-data-set'])

# get_triangles(['Ding-data-set', 'Kyu-data-set'], "similarity/triangles_C.csv", "C")

# get_triangles(['Ding-data-set', 'Kyu-data-set'], "similarity/triangles_G.csv", "  G")
#
# get_triangles(['Ding-data-set', 'Kyu-data-set'], "similarity/triangles_U.csv", "  U")

def basic(dirs, base):
    files = []
    excluded_files = 0
    for dir in dirs:
        files += glob.glob(dir + '/*.pdb')
        for file in files:
            try:
                RNAstruct = RNA(file)
                for chain in RNAstruct.chains:
                    for residue in chain.residues:
                        if residue.kind == base:
                            print(base)
            except:
                excluded_files += 1
                print('Warning was raised as an exception! ' + str(file) + ' excluded from analysis')
    print('Excluded files: ', excluded_files)


# basic(['Ding-data-set', 'Kyu-data-set'], '  C')
