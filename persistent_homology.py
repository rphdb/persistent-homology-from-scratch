"""Persistent homology from scratch: boundary matrix reduction and barcodes (Z/2Z coefficients).

Usage:
    python persistent_homology.py filtrations/torus.txt
    python persistent_homology.py filtrations/torus.txt -o outputs/torus_barcode.txt --time
    python persistent_homology.py --generate
"""

import argparse
import os
import time
from itertools import combinations


# 1. Reading the filtration file

class Simplex:
    """Represents one simplex with 3 fields: filtration value (val), dimension (dim) and the sorted list of vertices (vert)."""
    def __init__(self, val, dim, vert):
        self.val = val
        self.dim = dim
        self.vert = vert


def read_data(filename):
    """read_data(filename) -> list(Simplex). Reads the filtration file and returns the list of the corresponding simplices, sorted by filtration value (and dimension in case of equality)."""
    P = []
    with open(filename, "r") as data:
        for line in data:
            num = line.split()
            if not num:
                continue
            P.append(Simplex(float(num[0]), int(num[1]), sorted([int(v) for v in num[2:]])))
    return sorted(P, key=lambda i: (i.val, i.dim))


# 2. Computing the boundary matrix

# 2.1. Dense and naive version

def compute_boundary_naive(v):
    """compute_boundary_naive(v) -> list(list(int)). Returns the dense boundary matrix of the list of simplices v."""
    n = len(v)
    M = []
    for x in v:
        Mx = [0] * n
        for i in range(n):
            if (v[i].dim == x.dim - 1) and (set(v[i].vert) <= set(x.vert)):
                Mx[i] = 1
        M.append(Mx)
    return M


# 2.2. Sparse and optimized version

def compute_boundary(v):
    """compute_boundary(v) -> (list(set(int)), list(int), list(float)). Returns the sparse boundary matrix, the list of dimensions of the simplices v and the list of filtration values of the simplices v."""
    n = len(v)
    c = [set() for _ in range(n)]
    d = [x.dim for x in v]
    f = [x.val for x in v]
    index = {tuple(x.vert): i for i, x in enumerate(v)}
    for j, x in enumerate(v):
        if x.dim > 0:
            c[j] = {index[face] for face in combinations(x.vert, len(x.vert) - 1)}
    return (c, d, f)


# 3. Reduction algorithm and computation of the barcode

# 3.1. Dense and naive version

def last_one(x):
    """last_one(x) -> int. Returns the index of the lowest 1 of the column x, or -1 if the whole column is null."""
    d = len(x)
    for i in range(d):
        if x[d - 1 - i] == 1:
            return d - 1 - i
    return -1


def reduction_algorithm_naive(M):
    """reduction_algorithm_naive(M) -> list(list(int)). Returns the reduction of the dense boundary matrix M, computed with naive Gaussian elimination."""
    d = len(M)
    for jy in range(d):
        while True:
            k = last_one(M[jy])
            if k == -1:
                break
            u = 0
            while u < jy:
                if last_one(M[u]) == k:
                    M[jy] = [(M[jy][l] + M[u][l]) % 2 for l in range(d)]
                    break
                u += 1
            if u == jy:
                break
    return M


# 3.2. Sparse and optimized version

def reduction_algorithm(c):
    """reduction_algorithm(c) -> list(set(int)). Returns the reduction of the sparse boundary matrix c, computed with Gaussian elimination.

    Columns are sets of row indices, so adding two columns over Z/2Z is a symmetric difference (^=),
    and the pivot array gives in O(1) the column whose lowest 1 is in a given row."""
    n = len(c)
    pivot = [-1 for _ in range(n)]
    for j in range(n):
        while c[j] and pivot[max(c[j])] != -1:
            c[j] ^= c[pivot[max(c[j])]]
        if c[j]:
            pivot[max(c[j])] = j
    return c


def barcode(c, d, f, output_path):
    """barcode(c, d, f, output_path) -> None. Reads the persistence pairs from the reduced matrix c and writes the barcode to output_path, one interval per line: dimension, birth, death (inf if the feature never dies)."""
    n = len(c)
    killed = [False for _ in range(n)]
    with open(output_path, "w") as out:
        for i in range(n):
            j = n - 1 - i
            if not killed[j]:
                if c[j]:
                    out.write(f"{d[j] - 1} {f[max(c[j])]} {f[j]}\n")
                    killed[max(c[j])] = True
                else:
                    out.write(f"{d[j]} {f[j]} inf\n")


# 4. Experiments

def list_to_string(l):
    """list_to_string(l) -> str. Returns a string of the elements of l, in a format compatible with our filtration format."""
    s = ""
    for i in range(len(l)):
        s += f"{l[i]} "
    s += "\n"
    return s


def filtrations_spheres(folder="filtrations"):
    """filtrations_spheres() -> None. Writes filtrations for the d-spheres for d = 0 to d = 10 in the given folder."""
    for i in range(11):
        with open(os.path.join(folder, f"filtration_sphere_{i}.txt"), "w") as fi:
            for u in range(i + 2):
                fi.write(f"0 0 {u}\n")
            l = range(i + 2)
            for k in range(1, i + 1):
                for t in combinations(l, k + 1):
                    fi.write(f"{k} {k} ")
                    fi.write(list_to_string(t))


def filtrations_balls(folder="filtrations"):
    """filtrations_balls() -> None. Writes filtrations for the (d+1)-balls for d = 0 to d = 10 in the given folder (requires the sphere filtrations)."""
    for i in range(11):
        with open(os.path.join(folder, f"filtration_ball_{i}.txt"), "w") as fi:
            with open(os.path.join(folder, f"filtration_sphere_{i}.txt"), "r") as f_sphere_i:
                fi.write(f_sphere_i.read())
                fi.write(f"{i + 1} {i + 1} ")
                fi.write(list_to_string(range(i + 2)))


def compute_barcode(input_path, output_path):
    """compute_barcode(input_path, output_path) -> None. Full pipeline: read the filtration, build the sparse boundary matrix, reduce it and write the barcode."""
    c, d, f = compute_boundary(read_data(input_path))
    barcode(reduction_algorithm(c), d, f, output_path)


def main():
    parser = argparse.ArgumentParser(description="Compute the persistence barcode of a simplicial filtration (Z/2Z coefficients).")
    parser.add_argument("filtration", nargs="?", help="filtration file, one simplex per line: value dimension vertices...")
    parser.add_argument("-o", "--output", help="output barcode file (default: outputs/barcode_<name>.txt)")
    parser.add_argument("--time", action="store_true", help="print the execution time")
    parser.add_argument("--generate", action="store_true", help="generate the sphere and ball filtrations in filtrations/")
    args = parser.parse_args()

    if args.generate:
        os.makedirs("filtrations", exist_ok=True)
        filtrations_spheres()
        filtrations_balls()
        print("Generated filtrations/filtration_sphere_{0..10}.txt and filtration_ball_{0..10}.txt")
        if not args.filtration:
            return
    if not args.filtration:
        parser.error("give a filtration file, or use --generate")

    output = args.output
    if output is None:
        name = os.path.splitext(os.path.basename(args.filtration))[0]
        output = os.path.join("outputs", f"barcode_{name}.txt")
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)

    start = time.perf_counter()
    compute_barcode(args.filtration, output)
    elapsed = time.perf_counter() - start
    print(f"Barcode written to {output}")
    if args.time:
        print(f"Execution time: {elapsed:.4f} s")


if __name__ == "__main__":
    main()
