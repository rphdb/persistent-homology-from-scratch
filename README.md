# Persistent Homology from Scratch

A from-scratch Python implementation of **persistent homology** with coefficients in ℤ/2ℤ: it reads a simplicial filtration, builds its boundary matrix, reduces it with Gaussian elimination and outputs the **persistence barcode**. Pure Python, no dependencies.

📄 **Full report:** [`reports/report.pdf`](reports/report.pdf)

## How it works

1. **Read the filtration.** Each simplex comes with a filtration value; simplices are sorted by value, then by dimension.
2. **Build the boundary matrix.** A dense N × N matrix is out of reach: for a filtration of ~430k simplices it would take about 1.5 TB. The matrix is stored sparsely instead, one set of facet indices per column, built in O(N·k) with a dictionary from vertex tuples to simplex indices (k = maximal simplex dimension).
3. **Reduce it.** Standard column reduction over ℤ/2ℤ, where adding two columns is a set symmetric difference (`^=`, implemented in C) and a pivot array finds the column to add in O(1). Switching from a hand-written merge of sorted lists to Python sets made the reduction about 10× faster.
4. **Read the barcode.** Each pivot pairs a birth with a death; unpaired columns give infinite bars.

Naive dense versions of steps 2 and 3 (`compute_boundary_naive`, `reduction_algorithm_naive`) are kept for comparison.

## Results

**Homology of classic spaces.** The infinite bars of each barcode give the Betti numbers over ℤ/2ℤ, which match the theory:

| Space | Betti numbers | Filtration |
|---|---|---|
| d-sphere, d = 1…10 | β₀ = 1, β_d = 1 | generated |
| 0-sphere (two points) | β₀ = 2 | generated |
| (d+1)-ball, d = 0…10 | β₀ = 1 | generated |
| Torus | β₀ = 1, β₁ = 2, β₂ = 1 | hand-built |
| Möbius band | β₀ = 1, β₁ = 1 | hand-built |
| Projective plane | β₀ = 1, β₁ = 1, β₂ = 1 | hand-built |
| Klein bottle | β₀ = 1, β₁ = 2, β₂ = 1 | hand-built |

All other Betti numbers are 0.

**Performance** on the four larger filtrations provided by the course (not included in this repository), from reading the file to writing the barcode:

| Filtration | Simplices | Time |
|---|---:|---:|
| B | 108,161 | 6.4 s |
| C | 180,347 | 14.1 s |
| A | 428,643 | 4.3 min |
| D | 2,716,431 | 16 min |

The report also interprets the barcodes of these four filtrations as shapes (for instance, a curve winding around a torus that is progressively filled in), and gives the time and memory complexity of each function.

## Usage

Requires Python 3.8+.

```bash
# Compute a barcode (written to outputs/barcode_torus.txt by default)
python persistent_homology.py filtrations/torus.txt

# Choose the output file and print the execution time
python persistent_homology.py filtrations/torus.txt -o torus_barcode.txt --time

# Regenerate the sphere and ball filtrations (d = 0 to 10)
python persistent_homology.py --generate
```

**Input format:** one simplex per line, `value dimension vertex_1 ... vertex_k`, e.g. `2.0 1 0 3` for the edge {0, 3} entering at time 2.

**Output format:** one bar per line, `dimension birth death`, with `inf` for features that never die. For the torus:

```
2 3.0 inf
1 2.0 inf
1 2.0 inf
0 1.0 inf
...
```

## Repository structure

| Path | Content |
|---|---|
| [`persistent_homology.py`](persistent_homology.py) | The whole implementation and the command-line interface |
| [`filtrations/`](filtrations) | Generated filtrations (spheres and balls, d = 0 to 10) and hand-built ones (torus, Möbius band, projective plane, Klein bottle, small example) |
| [`outputs/`](outputs) | The barcodes of these filtrations |
| [`reports/`](reports/report.pdf) | Report: implementation, complexity analysis, answers and barcode interpretations |

## Acknowledgements

Developed by Raphaël des Boscs and Karel Nyobé as a lab assignment for the *Topological Data Analysis* course at École Polytechnique (2026). The problem statement and the test filtrations A–D were provided by the course and are not included in this repository; the implementation, the generated and hand-built filtrations, and the report are our own work.
