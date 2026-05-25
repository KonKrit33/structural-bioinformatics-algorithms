import math
import random


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
# PDB parsing
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
# Pairwise distances and d-RMSD
# ============================================================

def build_pair_records(points1, points2, residue_start=50):
    """
    Builds all corresponding pairwise distances for the two segments.

    Returns a list of dicts with:
      res_i, res_j, d1, d2, sq_err, avg_dist
    """
    n = len(points1)
    records = []

    for i in range(n):
        for j in range(i + 1, n):
            d1 = distance(points1[i], points1[j])
            d2 = distance(points2[i], points2[j])

            records.append({
                "res_i": residue_start + i,
                "res_j": residue_start + j,
                "d1": d1,
                "d2": d2,
                "sq_err": (d1 - d2) ** 2,
                "avg_dist": 0.5 * (d1 + d2)   # symmetric choice for "shortest distances"
            })

    return records


def d_rmsd_from_records(records):
    if len(records) == 0:
        raise ValueError("Empty subset of distances.")
    s = sum(r["sq_err"] for r in records)
    return math.sqrt(s / len(records))


def rel_error_percent(approx, exact):
    if abs(exact) < 1e-15:
        return 0.0
    return 100.0 * abs(approx - exact) / abs(exact)


def accuracy_comment(rel_err_percent):
    if rel_err_percent < 5:
        return "excellent approximation"
    elif rel_err_percent < 15:
        return "good approximation"
    elif rel_err_percent < 30:
        return "moderate approximation"
    else:
        return "poor approximation"


# ============================================================
# Approximation strategies
# ============================================================

def shortest_distance_subsets(records, subset_sizes):
    """
    Select subsets using the shortest corresponding distances.
    To be symmetric across the two structures, pairs are sorted
    by avg_dist = (d1 + d2)/2.
    """
    sorted_records = sorted(records, key=lambda r: r["avg_dist"])
    results = []

    for m in subset_sizes:
        subset = sorted_records[:m]
        approx = d_rmsd_from_records(subset)
        results.append((m, approx))

    return results


def random_distance_subsets_single_sample(records, subset_sizes, seed=20260405):
    """
    One reproducible random subset for each size.
    """
    results = []

    for m in subset_sizes:
        rng = random.Random(seed + m)
        subset = rng.sample(records, m)
        approx = d_rmsd_from_records(subset)
        results.append((m, approx, seed + m))

    return results


def random_distance_subsets_trials(records, subset_sizes, trials=200, seed=20260405):
    """
    Many random trials per subset size, for a more robust accuracy comment.
    """
    results = []

    for m in subset_sizes:
        values = []
        for t in range(trials):
            rng = random.Random(seed + 10000*m + t)
            subset = rng.sample(records, m)
            values.append(d_rmsd_from_records(subset))

        mean_val = sum(values) / len(values)
        var_val = sum((x - mean_val) ** 2 for x in values) / len(values)
        std_val = math.sqrt(var_val)
        min_val = min(values)
        max_val = max(values)

        results.append((m, mean_val, std_val, min_val, max_val))

    return results


# ============================================================
# Main
# ============================================================

def main():
    pdb1_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\6LU7.pdb"
    pdb2_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\2AMD.pdb"

    residue_start = 50
    residue_end = 99

    # From parts (a) and (b), both selected chains were A
    chain1 = "A"
    chain2 = "A"

    subset_sizes = [500, 200, 100]
    random_trials = 200

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

    # --------------------------------------------------------
    # Exact d-RMSD using all 1225 pairwise distances
    # --------------------------------------------------------
    records = build_pair_records(X, Y, residue_start)
    k = len(records)  
    exact_drmsd = d_rmsd_from_records(records)

    # --------------------------------------------------------
    # Approximation 1: shortest distances
    # --------------------------------------------------------
    shortest_results = shortest_distance_subsets(records, subset_sizes)

    # --------------------------------------------------------
    # Approximation 2: one reproducible random subset for each size
    # --------------------------------------------------------
    random_single_results = random_distance_subsets_single_sample(records, subset_sizes)

    # --------------------------------------------------------
    # Extra: many random trials per size, to comment on stability
    # --------------------------------------------------------
    random_trial_results = random_distance_subsets_trials(
        records, subset_sizes, trials=random_trials
    )

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------
    print("====================================================")
    print("Problem 2(c): Approximate d-RMSD using subsets of CA distances")
    print("====================================================")
    print()
    print(f"Structure 1: {pdb1_path}")
    print(f"Selected chain in structure 1: {selected_chain1}")
    print(f"Structure 2: {pdb2_path}")
    print(f"Selected chain in structure 2: {selected_chain2}")
    print()
    print(f"Residue range used: {residue_start}-{residue_end}")
    print(f"Number of CA atoms in each segment: {len(X)}")
    print(f"Total number of pairwise distances: k = C(50,2) = {k}")
    print()
    print(f"Exact d-RMSD using all {k} pairwise distances = {exact_drmsd:.10f} Angstrom")
    print()

    print("----------------------------------------------------")
    print("Approximation 1: shortest distances")
    print("(pairs sorted by average corresponding distance (d_ij + d'_ij)/2)")
    print("----------------------------------------------------")
    for m, approx in shortest_results:
        abs_err = abs(approx - exact_drmsd)
        rel_err = rel_error_percent(approx, exact_drmsd)
        comment = accuracy_comment(rel_err)
        print(
            f"Subset size {m:>3d}: approx d-RMSD = {approx:.10f} A   |   "
            f"abs error = {abs_err:.10f}   |   rel error = {rel_err:.4f}%   |   {comment}"
        )
    print()

    print("----------------------------------------------------")
    print("Approximation 2: randomly selected distances")
    print("(one reproducible random subset for each size)")
    print("----------------------------------------------------")
    for m, approx, used_seed in random_single_results:
        abs_err = abs(approx - exact_drmsd)
        rel_err = rel_error_percent(approx, exact_drmsd)
        comment = accuracy_comment(rel_err)
        print(
            f"Subset size {m:>3d}: approx d-RMSD = {approx:.10f} A   |   "
            f"abs error = {abs_err:.10f}   |   rel error = {rel_err:.4f}%   |   "
            f"seed = {used_seed}   |   {comment}"
        )
    print()

    print("----------------------------------------------------")
    print(f"Random-subset stability over {random_trials} trials per subset size")
    print("----------------------------------------------------")
    for m, mean_val, std_val, min_val, max_val in random_trial_results:
        mean_abs_err = abs(mean_val - exact_drmsd)
        mean_rel_err = rel_error_percent(mean_val, exact_drmsd)
        comment = accuracy_comment(mean_rel_err)
        print(
            f"Subset size {m:>3d}: mean = {mean_val:.10f} A, std = {std_val:.10f}, "
            f"min = {min_val:.10f}, max = {max_val:.10f}   |   "
            f"mean rel error = {mean_rel_err:.4f}%   |   {comment}"
        )
    print()

    print("General expectation:")
    print("- larger subsets should usually give more accurate approximations")
    print("- shortest-distance subsets often behave differently from random subsets")
    print("- random subsets may vary depending on the sample, especially for 100 distances")
    print()

if __name__ == "__main__":
    main()