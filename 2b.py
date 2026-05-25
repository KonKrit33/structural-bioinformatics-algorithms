import math


# ============================================================
# Basic vector utilities
# ============================================================

def vec_sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(v):
    return math.sqrt(dot(v, v))

def distance(a, b):
    return norm(vec_sub(a, b))


# ============================================================
# PDB parsing: keep only CA atoms for residues 50-99
# ============================================================

def parse_pdb_ca_atoms(pdb_path, residue_start=50, residue_end=99, chain_id=None):
    """
    Extract CA atoms for residues residue_start..residue_end inclusive.
    Returns:
        selected_chain_id, coords_dict
    where coords_dict[resseq] = [x, y, z]
    """
    chains = {}

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            alt_loc = line[16].strip()
            if alt_loc not in ("", "A"):
                continue

            current_chain = line[21].strip()
            if current_chain == "":
                current_chain = "_"

            try:
                resseq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue

            if current_chain not in chains:
                chains[current_chain] = {}

            if resseq not in chains[current_chain]:
                chains[current_chain][resseq] = [x, y, z]

    needed = set(range(residue_start, residue_end + 1))

    if chain_id is not None:
        if chain_id not in chains:
            raise ValueError(f"Chain '{chain_id}' not found in {pdb_path}")
        missing = sorted(list(needed - set(chains[chain_id].keys())))
        if missing:
            raise ValueError(
                f"Chain '{chain_id}' in {pdb_path} is missing residues: {missing}"
            )
        return chain_id, chains[chain_id]

    # automatic chain selection: first chain with full coverage
    for ch in sorted(chains.keys()):
        if needed.issubset(chains[ch].keys()):
            return ch, chains[ch]

    coverage_info = {
        ch: len(needed.intersection(set(res_dict.keys())))
        for ch, res_dict in chains.items()
    }
    raise ValueError(
        f"No chain in {pdb_path} contains all residues {residue_start}-{residue_end}. "
        f"Coverage by chain: {coverage_info}"
    )


def get_ordered_segment(coords_dict, residue_start=50, residue_end=99):
    return [coords_dict[r] for r in range(residue_start, residue_end + 1)]


# ============================================================
# d-RMSD using ALL pairwise distances
# ============================================================

def all_pairwise_distances(points, residue_start=50):
    """
    Returns list of tuples:
        (res_i, res_j, dist)
    for all i < j
    """
    pairs = []
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(points[i], points[j])
            res_i = residue_start + i
            res_j = residue_start + j
            pairs.append((res_i, res_j, d))
    return pairs


def d_rmsd_all_pairs(points1, points2, residue_start=50):
    """
    Uses all pairwise distances between the 50 CA atoms.
    Returns:
        k, d_rmsd, pairs1, pairs2
    """
    pairs1 = all_pairwise_distances(points1, residue_start)
    pairs2 = all_pairwise_distances(points2, residue_start)

    if len(pairs1) != len(pairs2):
        raise ValueError("The two segments do not have the same number of pairwise distances.")

    k = len(pairs1)
    s = 0.0

    for (r1a, r1b, d1), (r2a, r2b, d2) in zip(pairs1, pairs2):
        if (r1a, r1b) != (r2a, r2b):
            raise ValueError("Residue-pair mismatch between the two structures.")
        s += (d1 - d2) ** 2

    drmsd = math.sqrt(s / k)
    return k, drmsd, pairs1, pairs2


# ============================================================
# Main
# ============================================================

def main():
    pdb1_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\6LU7.pdb"
    pdb2_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\2AMD.pdb"

    residue_start = 50
    residue_end = 99

    # Since from part (a) I found chain A in both:
    chain1 = "A"
    chain2 = "A"

    # --------------------------------------------------------
    # Parse structures
    # --------------------------------------------------------
    selected_chain1, coords1 = parse_pdb_ca_atoms(
        pdb1_path, residue_start, residue_end, chain1
    )
    selected_chain2, coords2 = parse_pdb_ca_atoms(
        pdb2_path, residue_start, residue_end, chain2
    )

    X = get_ordered_segment(coords1, residue_start, residue_end)
    Y = get_ordered_segment(coords2, residue_start, residue_end)

    if len(X) != 50 or len(Y) != 50:
        raise ValueError("Expected 50 CA atoms from each structure for residues 50-99.")

    k, drmsd, pairs1, pairs2 = d_rmsd_all_pairs(X, Y, residue_start)

    print("====================================================")
    print("Problem 2(b): d-RMSD using all pairwise CA distances")
    print("====================================================")
    print()
    print(f"Structure 1: {pdb1_path}")
    print(f"Selected chain in structure 1: {selected_chain1}")
    print(f"Structure 2: {pdb2_path}")
    print(f"Selected chain in structure 2: {selected_chain2}")
    print()
    print(f"Residue range used: {residue_start}-{residue_end}")
    print(f"Number of CA atoms in each segment: {len(X)}")
    print()

    print("Number of all pairwise CA distances:")
    print(f"k = C(50,2) = 50*49/2 = {k}")
    print()

    print("d-RMSD formula:")
    print("d-RMSD = sqrt( (1/k) * sum_{i<j} (d_ij - d'_ij)^2 )")
    print()

    print(f"Final d-RMSD = {drmsd:.10f} Angstrom")
    print()

if __name__ == "__main__":
    main()