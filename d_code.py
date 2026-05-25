from pathlib import Path
from collections import defaultdict
import csv
import math

ROOT = Path(r"C:\bio")
RESULTS_DIR = ROOT / "results_b"

BLAST_FILE = RESULTS_DIR / "blastp.tsv"
FASTA_FILE = RESULTS_DIR / "fasta36.tsv"
SSEARCH_FILE = RESULTS_DIR / "ssearch36.tsv"
GGSEARCH_FILE = RESULTS_DIR / "ggsearch36.tsv"


METHODS = ["BLAST", "FASTA", "SSEARCH", "GGSEARCH"]


def load_ranked_hits_blast(path: Path):
    """
    BLAST outfmt 6:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    Ranked by file order (already sorted by BLAST).
    """
    hits = defaultdict(list)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            q = parts[0]
            s = parts[1]
            if s not in hits[q]:
                hits[q].append(s)
    return hits


def load_ranked_hits_fasta_family(path: Path):
    """
    FASTA/SSEARCH/GGSEARCH with -m 8CBL:
    query_id, query_length, subject_id, subject_length, ...
    Ranked by file order.
    """
    hits = defaultdict(list)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            q = parts[0]
            s = parts[2]
            if s not in hits[q]:
                hits[q].append(s)
    return hits


def query_group(q_id: str):
    idx = int(q_id.replace("query", ""))
    if 0 <= idx <= 2:
        return "strong_query"
    elif 3 <= idx <= 5:
        return "weak_query"
    elif 6 <= idx <= 8:
        return "random_query"
    elif 9 <= idx <= 11:
        return "fragment_query"
    elif 12 <= idx <= 14:
        return "motif_query"
    else:
        return "unknown"


def top1_summary(all_hits):
    """
    For each query, report top-1 returned by each method.
    """
    rows = []
    all_queries = sorted(
        set().union(*[set(v.keys()) for v in all_hits.values()]),
        key=lambda x: int(x.replace("query", ""))
    )

    for q in all_queries:
        row = {
            "query": q,
            "query_group": query_group(q),
            "BLAST_top1": all_hits["BLAST"].get(q, [None])[0] if all_hits["BLAST"].get(q) else "",
            "FASTA_top1": all_hits["FASTA"].get(q, [None])[0] if all_hits["FASTA"].get(q) else "",
            "SSEARCH_top1": all_hits["SSEARCH"].get(q, [None])[0] if all_hits["SSEARCH"].get(q) else "",
            "GGSEARCH_top1": all_hits["GGSEARCH"].get(q, [None])[0] if all_hits["GGSEARCH"].get(q) else "",
        }
        rows.append(row)
    return rows


def spearman_top5(list1, list2):
    """
    Compute Spearman rank correlation on top-5 sets.
    If an item is absent from a method's top-5, assign rank = 6.
    Union of items from the two top-5 lists is used.
    """
    top1 = list1[:5]
    top2 = list2[:5]
    items = sorted(set(top1) | set(top2))

    n = len(items)
    if n <= 1:
        return 1.0

    rank1 = {}
    rank2 = {}

    for i, x in enumerate(top1, start=1):
        rank1[x] = i
    for i, x in enumerate(top2, start=1):
        rank2[x] = i

    for x in items:
        if x not in rank1:
            rank1[x] = 6
        if x not in rank2:
            rank2[x] = 6

    d2_sum = 0
    for x in items:
        d = rank1[x] - rank2[x]
        d2_sum += d * d

    rho = 1 - (6 * d2_sum) / (n * (n * n - 1))
    return rho


def pairwise_spearman_matrix(all_hits):
    """
    Compute average Spearman top-5 correlation for every pair of methods,
    averaging across all queries.
    """
    all_queries = sorted(
        set().union(*[set(v.keys()) for v in all_hits.values()]),
        key=lambda x: int(x.replace("query", ""))
    )

    per_query_rows = []
    matrix = {m1: {m2: None for m2 in METHODS} for m1 in METHODS}

    for m1 in METHODS:
        for m2 in METHODS:
            values = []
            for q in all_queries:
                l1 = all_hits[m1].get(q, [])
                l2 = all_hits[m2].get(q, [])
                rho = spearman_top5(l1, l2)
                values.append(rho)

                if m1 < m2:
                    per_query_rows.append({
                        "query": q,
                        "query_group": query_group(q),
                        "method_1": m1,
                        "method_2": m2,
                        "spearman_top5": f"{rho:.6f}"
                    })

            avg_rho = sum(values) / len(values) if values else float("nan")
            matrix[m1][m2] = avg_rho

    return matrix, per_query_rows


def save_csv(path: Path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def save_matrix_csv(path: Path, matrix):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["method"] + METHODS)
        for m1 in METHODS:
            row = [m1]
            for m2 in METHODS:
                row.append(f"{matrix[m1][m2]:.6f}")
            w.writerow(row)


def print_matrix(matrix):
    print("\n=== Average Spearman top-5 correlation matrix ===")
    header = ["method"] + METHODS
    print("\t".join(header))
    for m1 in METHODS:
        row = [m1] + [f"{matrix[m1][m2]:.4f}" for m2 in METHODS]
        print("\t".join(row))


def main():
    blast_hits = load_ranked_hits_blast(BLAST_FILE)
    fasta_hits = load_ranked_hits_fasta_family(FASTA_FILE)
    ssearch_hits = load_ranked_hits_fasta_family(SSEARCH_FILE)
    ggsearch_hits = load_ranked_hits_fasta_family(GGSEARCH_FILE)

    all_hits = {
        "BLAST": blast_hits,
        "FASTA": fasta_hits,
        "SSEARCH": ssearch_hits,
        "GGSEARCH": ggsearch_hits,
    }

    # Part 1: top-1 comparison
    top1_rows = top1_summary(all_hits)
    save_csv(RESULTS_DIR / "d_top1_summary.csv", top1_rows)

    # Print top-1 summary briefly
    print("=== Top-1 summary ===")
    for row in top1_rows:
        print(
            row["query"],
            row["query_group"],
            row["BLAST_top1"],
            row["FASTA_top1"],
            row["SSEARCH_top1"],
            row["GGSEARCH_top1"],
            sep="\t"
        )

    # Part 2: pairwise top-5 Spearman
    matrix, per_query_rows = pairwise_spearman_matrix(all_hits)
    save_csv(RESULTS_DIR / "d_pairwise_spearman_per_query.csv", per_query_rows)
    save_matrix_csv(RESULTS_DIR / "d_pairwise_spearman_matrix.csv", matrix)

    print_matrix(matrix)
    print("\nFiles written in:", RESULTS_DIR)


if __name__ == "__main__":
    main()