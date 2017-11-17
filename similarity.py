from CABS.cabsDock.vector3d import Vector3d
import numpy as np
import rmsd
import itertools
from scipy.cluster.hierarchy import fcluster, linkage


class Structure:
    def __init__(self, points):
        self.points = points  # [] each point is Vector3d

    def addPoint(self, v):
        self.points.append(v)

    def __str__(self):
        s = "Structure\n"
        for p in self.points:
            s += "\t"
            s += p.__str__()
        return s


### SIMILARITY METHODS

def count_rmsd(t1, t2):
    P = np.array([[float(v.x), float(v.y), float(v.z)] for v in t1.points])
    Q = np.array([[float(v.x), float(v.y), float(v.z)] for v in t2.points])

    # print("RMSD before translation: ", round(rmsd.kabsch_rmsd(P, Q), 6))
    P -= rmsd.centroid(P)
    Q -= rmsd.centroid(Q)
    # print("RMSD after translation: ", round(rmsd.kabsch_rmsd(P, Q), 6))
    return round(rmsd.kabsch_rmsd(P, Q), 4)  # round 6 in triangles_matrix


### CLUSTERING METHODS

def hierarchical_clustering(sim_array, threshold=None):
    tri = np.triu(sim_array, k=1)  # get upper triangle from array without diagonal
    if not threshold:
        flatten = np.array(list(itertools.chain.from_iterable(tri.tolist())))  # flatten list from tri
        var = np.var(flatten)
        print(var) #--> variance 1.75
        print(np.mean(flatten), np.std(flatten)) # --> std 1.32, mean 0.74
        threshold = np.mean(flatten)
        print(np.max(flatten))  # --> 14.82
    return
    linkage_matrix = linkage(tri)
    clusters = fcluster(linkage_matrix, threshold)
    print(clusters)
    print(len(set(clusters)))  # 2567
    return clusters


def k_medoids(sim_array): #TODO k = sqrt(n)/3 ?
    pass


### WORKFLOW
##  - creating similarity matrix from .csv file
##  - clustering on matrix


def create_similarity_matrix(coords_file, similarity_method):
    structs = []
    with open(coords_file, "r") as f:
        for line in f.readlines():
            l = [line.split(",")[i:i + 3] for i in range(0, len(line.split(",")), 3)][:3]

            structs.append(Structure([Vector3d(float(trio[0]), float(trio[1]), float(trio[2])) for trio in l]))
    sim_array = np.zeros((len(structs), len(structs)))
    for i in range(len(structs)):
        for j in range(len(structs)):
            # print(structs[i])
            # print(structs[j])
            sim_array[i][j] = similarity_method(structs[i], structs[j])
    np.savetxt("similarity/triangles_matrix_U.csv", sim_array, delimiter=",")


def cluster(similarity_matrix_file, clustering_method, treshold=None):
    sim_array = np.genfromtxt(similarity_matrix_file, delimiter=',')
    clustering_method(sim_array, treshold)


# create_similarity_matrix("similarity/triangles_U.csv", count_rmsd)
cluster("similarity/triangles_matrix.csv", hierarchical_clustering)


