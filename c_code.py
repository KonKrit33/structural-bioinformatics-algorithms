from pathlib import Path
import csv
from collections import defaultdict

ROOT = Path(r"C:\bio")
RESULTS_DIR = ROOT / "results_b"

BLAST_FILE = RESULTS_DIR / "blastp.tsv"
FASTA_FILE = RESULTS_DIR / "fasta36.tsv"
SSEARCH_FILE = RESULTS_DIR / "ssearch36.tsv"
GGSEARCH_FILE = RESULTS_DIR / "ggsearch36.tsv"


def load_hits_blast(path: Path):
    hits = defaultdict(set)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            qseqid = parts[0]
            sseqid = parts[1]
            hits[qseqid].add(sseqid)
    return hits


def load_hits_fasta_family(path: Path):
    hits = defaultdict(set)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            qseqid = parts[0]
            sseqid = parts[2]
            hits[qseqid].add(sseqid)
    return hits


def prot_group(prot_id: str):
    idx = int(prot_id.replace("prot", ""))
    if 0 <= idx <= 24:
        return "strong_db"
    elif 25 <= idx <= 49:
        return "weak_db"
    elif 50 <= idx <= 74:
        return "random_db"
    elif 75 <= idx <= 99:
        return "motif_db"
    else:
        return "unknown"


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


def summarize_method_losses(dp_gold, heuristic_hits, allowed_query_groups=None, allowed_db_groups=None):
    rows = []
    total_gold = 0
    total_missed = 0

    for q in sorted(dp_gold.keys(), key=lambda x: int(x.replace("query", ""))):
        qg = query_group(q)

        if allowed_query_groups is not None and qg not in allowed_query_groups:
            continue

        gold = set(dp_gold[q])
        heur = set(heuristic_hits.get(q, set()))

        if allowed_db_groups is not None:
            gold = {x for x in gold if prot_group(x) in allowed_db_groups}
            heur = {x for x in heur if prot_group(x) in allowed_db_groups}

        if len(gold) == 0:
            missed = set()
            pct = 0.0
        else:
            missed = gold - heur
            pct = 100.0 * len(missed) / len(gold)

        rows.append({
            "query": q,
            "query_group": qg,
            "gold_hits": len(gold),
            "missed_hits": len(missed),
            "missed_pct": f"{pct:.2f}",
            "missed_ids": ",".join(sorted(missed))
        })

        total_gold += len(gold)
        total_missed += len(missed)

    total_pct = 0.0 if total_gold == 0 else 100.0 * total_missed / total_gold
    return rows, total_gold, total_missed, total_pct


def save_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    blast_hits = load_hits_blast(BLAST_FILE)
    fasta_hits = load_hits_fasta_family(FASTA_FILE)
    ssearch_hits = load_hits_fasta_family(SSEARCH_FILE)
    ggsearch_hits = load_hits_fasta_family(GGSEARCH_FILE)

    # DP gold standard = union of exact methods
    all_queries = set(ssearch_hits.keys()) | set(ggsearch_hits.keys())
    dp_gold = {}
    for q in all_queries:
        dp_gold[q] = set(ssearch_hits.get(q, set())) | set(ggsearch_hits.get(q, set()))

    # 1) OVERALL analysis: all queries, all DB groups
    fasta_rows_overall, fasta_gold_overall, fasta_missed_overall, fasta_pct_overall = summarize_method_losses(
        dp_gold, fasta_hits
    )
    blast_rows_overall, blast_gold_overall, blast_missed_overall, blast_pct_overall = summarize_method_losses(
        dp_gold, blast_hits
    )

    # 2) CLEAN WEAK analysis: only weak queries AND only weak_db sequences
    fasta_rows_weak_clean, fasta_gold_weak_clean, fasta_missed_weak_clean, fasta_pct_weak_clean = summarize_method_losses(
        dp_gold,
        fasta_hits,
        allowed_query_groups={"weak_query"},
        allowed_db_groups={"weak_db"},
    )
    blast_rows_weak_clean, blast_gold_weak_clean, blast_missed_weak_clean, blast_pct_weak_clean = summarize_method_losses(
        dp_gold,
        blast_hits,
        allowed_query_groups={"weak_query"},
        allowed_db_groups={"weak_db"},
    )

    save_csv(RESULTS_DIR / "c_fasta_overall.csv", fasta_rows_overall)
    save_csv(RESULTS_DIR / "c_blast_overall.csv", blast_rows_overall)
    save_csv(RESULTS_DIR / "c_fasta_clean_weak.csv", fasta_rows_weak_clean)
    save_csv(RESULTS_DIR / "c_blast_clean_weak.csv", blast_rows_weak_clean)

    print("=== OVERALL vs DP GOLD ===")
    print(f"FASTA missed {fasta_missed_overall}/{fasta_gold_overall} = {fasta_pct_overall:.2f}%")
    print(f"BLAST missed {blast_missed_overall}/{blast_gold_overall} = {blast_pct_overall:.2f}%")

    print("\n=== CLEAN WEAK ANALYSIS (weak queries vs weak_db only) ===")
    print(f"FASTA missed weak homologs {fasta_missed_weak_clean}/{fasta_gold_weak_clean} = {fasta_pct_weak_clean:.2f}%")
    print(f"BLAST missed weak homologs {blast_missed_weak_clean}/{blast_gold_weak_clean} = {blast_pct_weak_clean:.2f}%")

    print("\nCSV files written in:", RESULTS_DIR)


if __name__ == "__main__":
    main()