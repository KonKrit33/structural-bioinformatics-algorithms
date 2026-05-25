import math

# ============================================================
# Basic linear algebra utilities (pure Python, no external libs)
# ============================================================

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(v):
    return math.sqrt(dot(v, v))

def vec_add(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]

def vec_sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

def vec_scale(s, v):
    return [s*v[0], s*v[1], s*v[2]]

def mat_vec_mul(M, v):
    return [
        M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
        M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
        M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2],
    ]

def mat_mul(A, B):
    C = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            C[i][j] = sum(A[i][k] * B[k][j] for k in range(3))
    return C

def transpose(A):
    return [
        [A[0][0], A[1][0], A[2][0]],
        [A[0][1], A[1][1], A[2][1]],
        [A[0][2], A[1][2], A[2][2]],
    ]

def det3(A):
    return (
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
        - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
        + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0])
    )

def identity3():
    return [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]

def column(A, j):
    return [A[0][j], A[1][j], A[2][j]]

def set_column(A, j, v):
    A[0][j] = v[0]
    A[1][j] = v[1]
    A[2][j] = v[2]

def outer(u, v):
    return [
        [u[0]*v[0], u[0]*v[1], u[0]*v[2]],
        [u[1]*v[0], u[1]*v[1], u[1]*v[2]],
        [u[2]*v[0], u[2]*v[1], u[2]*v[2]],
    ]

def mat_sub(A, B):
    C = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            C[i][j] = A[i][j] - B[i][j]
    return C

def mat_scale(s, A):
    C = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            C[i][j] = s * A[i][j]
    return C


# ============================================================
# Jacobi eigenvalue algorithm for symmetric 3x3 matrices
# ============================================================

def jacobi_eigen_symmetric_3x3(A, max_iter=100, tol=1e-12):
    """
    Returns eigenvalues and eigenvectors of a symmetric 3x3 matrix A.
    Eigenvectors returned as columns of V.
    """
    D = [row[:] for row in A]
    V = identity3()

    for _ in range(max_iter):
        # Find largest off-diagonal absolute value
        p, q = 0, 1
        max_off = abs(D[0][1])
        for i in range(3):
            for j in range(i+1, 3):
                if abs(D[i][j]) > max_off:
                    max_off = abs(D[i][j])
                    p, q = i, j

        if max_off < tol:
            break

        if abs(D[p][p] - D[q][q]) < tol:
            angle = math.pi / 4.0
        else:
            tau = (D[q][q] - D[p][p]) / (2.0 * D[p][q])
            t = 1.0 / (abs(tau) + math.sqrt(1.0 + tau*tau))
            if tau < 0:
                t = -t
            angle = math.atan(t)

        c = math.cos(angle)
        s = math.sin(angle)

        # Apply Jacobi rotation to D: D' = J^T D J
        J = identity3()
        J[p][p] = c
        J[q][q] = c
        J[p][q] = s
        J[q][p] = -s

        JT = transpose(J)
        D = mat_mul(mat_mul(JT, D), J)
        V = mat_mul(V, J)

    eigenvalues = [D[0][0], D[1][1], D[2][2]]

    # Sort descending by eigenvalue
    idx = sorted(range(3), key=lambda i: eigenvalues[i], reverse=True)
    eigenvalues_sorted = [eigenvalues[i] for i in idx]

    V_sorted = [[0.0]*3 for _ in range(3)]
    for new_j, old_j in enumerate(idx):
        col = column(V, old_j)
        # normalize
        nrm = norm(col)
        if nrm > 1e-15:
            col = [x / nrm for x in col]
        set_column(V_sorted, new_j, col)

    return eigenvalues_sorted, V_sorted


# ============================================================
# 3x3 SVD via eigen-decomposition of H^T H
# ============================================================

def svd_3x3(H):
    """
    Compute SVD of a 3x3 matrix H:
        H = U * S * V^T
    using eigen-decomposition of H^T H.
    Returns U, singular_values, V.
    """
    HT = transpose(H)
    HTH = mat_mul(HT, H)  # symmetric positive semidefinite

    eigenvalues, V = jacobi_eigen_symmetric_3x3(HTH)

    singular_values = [math.sqrt(max(ev, 0.0)) for ev in eigenvalues]

    # Compute U = H * V * S^{-1}
    U = [[0.0]*3 for _ in range(3)]
    for j in range(3):
        vj = column(V, j)
        Hv = mat_vec_mul(H, vj)
        sj = singular_values[j]
        if sj > 1e-12:
            uj = [x / sj for x in Hv]
        else:
            # fallback: fill later if singular
            uj = [0.0, 0.0, 0.0]
        # normalize
        nrm = norm(uj)
        if nrm > 1e-15:
            uj = [x / nrm for x in uj]
        set_column(U, j, uj)

    # If rank deficient, complete orthonormal basis manually
    u0 = column(U, 0)
    u1 = column(U, 1)
    u2 = column(U, 2)

    # Fill missing vectors if needed
    if norm(u0) < 1e-10:
        u0 = [1.0, 0.0, 0.0]
    if norm(u1) < 1e-10:
        # choose vector orthogonal to u0
        candidate = [0.0, 1.0, 0.0]
        proj = vec_scale(dot(candidate, u0), u0)
        u1 = vec_sub(candidate, proj)
        if norm(u1) < 1e-10:
            candidate = [0.0, 0.0, 1.0]
            proj = vec_scale(dot(candidate, u0), u0)
            u1 = vec_sub(candidate, proj)
    # Gram-Schmidt
    u0 = vec_scale(1.0 / norm(u0), u0)
    u1 = vec_sub(u1, vec_scale(dot(u1, u0), u0))
    if norm(u1) < 1e-10:
        if abs(u0[0]) < 0.9:
            u1 = [1.0, 0.0, 0.0]
        else:
            u1 = [0.0, 1.0, 0.0]
        u1 = vec_sub(u1, vec_scale(dot(u1, u0), u0))
    u1 = vec_scale(1.0 / norm(u1), u1)

    # u2 = u0 x u1
    u2 = [
        u0[1]*u1[2] - u0[2]*u1[1],
        u0[2]*u1[0] - u0[0]*u1[2],
        u0[0]*u1[1] - u0[1]*u1[0]
    ]
    if norm(u2) > 1e-15:
        u2 = vec_scale(1.0 / norm(u2), u2)

    set_column(U, 0, u0)
    set_column(U, 1, u1)
    set_column(U, 2, u2)

    return U, singular_values, V


# ============================================================
# PDB parsing
# ============================================================

def parse_pdb_ca_atoms(pdb_path, residue_start=50, residue_end=99, chain_id=None):
    """
    Extract CA atoms for residues residue_start..residue_end inclusive.
    If chain_id is None, uses the first chain that contains the full set.
    Returns:
        selected_chain_id, coords_dict
    where coords_dict[resseq] = [x, y, z]
    """

    # First pass: collect CA atoms by chain and residue
    chains = {}

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            alt_loc = line[16].strip()
            # Keep blank or A altloc
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

            # Keep first encountered CA per residue
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

    # Auto-select first chain with complete coverage
    for ch in sorted(chains.keys()):
        if needed.issubset(chains[ch].keys()):
            return ch, chains[ch]

    # If no chain has full coverage, report what exists
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
# Alignment functions
# ============================================================

def centroid(points):
    n = len(points)
    c = [0.0, 0.0, 0.0]
    for p in points:
        c[0] += p[0]
        c[1] += p[1]
        c[2] += p[2]
    return [c[0]/n, c[1]/n, c[2]/n]

def center_points(points, c):
    return [vec_sub(p, c) for p in points]

def covariance_matrix(Xc, Yc):
    """
    H = X^T Y
    Xc, Yc are centered point lists, each n x 3.
    """
    H = [[0.0]*3 for _ in range(3)]
    for x, y in zip(Xc, Yc):
        for i in range(3):
            for j in range(3):
                H[i][j] += x[i] * y[j]
    return H

def kabsch_rotation(H):
    """
    Returns optimal proper rotation Q using SVD:
        H = U S V^T
        Q = U V^T
    If det(Q) < 0, flip 3rd column of U.
    """
    U, singular_values, V = svd_3x3(H)
    VT = transpose(V)
    Q = mat_mul(U, VT)

    if det3(Q) < 0:
        # flip third column of U
        U_flipped = [row[:] for row in U]
        for i in range(3):
            U_flipped[i][2] *= -1.0
        Q = mat_mul(U_flipped, VT)

    return Q, singular_values

def apply_transform(points, Q, t):
    transformed = []
    for p in points:
        qp = mat_vec_mul(Q, p)
        transformed.append(vec_add(qp, t))
    return transformed

def rmsd(A, B):
    n = len(A)
    s = 0.0
    for a, b in zip(A, B):
        d = vec_sub(a, b)
        s += dot(d, d)
    return math.sqrt(s / n)

def optimal_alignment(X, Y):
    """
    Given corresponding point sets X and Y (same length),
    compute:
      - centroids
      - optimal rotation Q
      - optimal translation t
      - c-RMSD
    mapping X onto Y:
      x -> Qx + t
    """
    xc = centroid(X)
    yc = centroid(Y)

    Xc = center_points(X, xc)
    Yc = center_points(Y, yc)

    H = covariance_matrix(Xc, Yc)
    Q, singular_values = kabsch_rotation(H)

    Qxc = mat_vec_mul(Q, xc)
    t = vec_sub(yc, Qxc)

    X_aligned = apply_transform(X, Q, t)
    crmsd = rmsd(X_aligned, Y)

    return {
        "xc": xc,
        "yc": yc,
        "H": H,
        "Q": Q,
        "t": t,
        "singular_values": singular_values,
        "c_rmsd": crmsd,
        "X_aligned": X_aligned,
    }


# ============================================================
# Pretty printing
# ============================================================

def print_matrix(name, M):
    print(name)
    for row in M:
        print("  [{: .10f}, {: .10f}, {: .10f}]".format(row[0], row[1], row[2]))
    print()

def print_vector(name, v):
    print("{} = [{: .10f}, {: .10f}, {: .10f}]".format(name, v[0], v[1], v[2]))

def main():

    pdb1_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\6LU7.pdb"
    pdb2_path = r"C:\Users\Konstantinos Krit\Desktop\MSc\Β' Εξάμηνο\Algorithms in Structural Bioinformatics\Project_A\Problem_2\2AMD.pdb"

    residue_start = 50
    residue_end = 99

    chain1 = None
    chain2 = None

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

    result = optimal_alignment(X, Y)

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------
    print("====================================================")
    print("Problem 2(a): Optimal alignment using CA atoms of residues 50-99")
    print("====================================================")
    print()
    print(f"Structure 1: {pdb1_path}")
    print(f"Selected chain in structure 1: {selected_chain1}")
    print(f"Structure 2: {pdb2_path}")
    print(f"Selected chain in structure 2: {selected_chain2}")
    print()
    print(f"Number of corresponding CA atoms used: {len(X)}")
    print()

    print_vector("Centroid x_c", result["xc"])
    print_vector("Centroid y_c", result["yc"])
    print()

    print_matrix("Covariance matrix H = X^T Y:", result["H"])
    print("Singular values:", ["{:.10f}".format(s) for s in result["singular_values"]])
    print()

    print_matrix("Optimal rotation matrix Q:", result["Q"])
    print("det(Q) = {:.10f}".format(det3(result["Q"])))
    print()

    print_vector("Optimal translation vector t", result["t"])
    print()

    print("Final c-RMSD = {:.10f} Angstrom".format(result["c_rmsd"]))
    print()

    print("Transformation that maps structure 1 segment onto structure 2 segment:")
    print("    x_aligned = Q x + t")
    print()

if __name__ == "__main__":
    main()