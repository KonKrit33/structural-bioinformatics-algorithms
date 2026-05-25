from pathlib import Path
import subprocess
import shutil
import time
import csv
from collections import defaultdict

ROOT = Path(r"C:\bio")
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results_e"
BLASTDB_DIR = RESULTS_DIR / "blastdb"

DB_FASTA = DATA_DIR / "db.fasta"
QUERIES_FASTA = DATA_DIR / "queries.fasta"

BLAST_BIN_DIR = ROOT / "data" / "bin"


def resolve_windows_exe(name: str, bindir=None) -> str:
    candidates = []
    if bindir is not None:
        candidates.append(Path(bindir) / f"{name}.exe")
        candidates.append(Path(bindir) / name)

    for c in candidates:
        if c.exists():
            return str(c)

    p = shutil.which(name) or shutil.which(f"{name}.exe")
    if p:
        return p

    raise FileNotFoundError(f"Executable not found: {name}")


def ensure_nonempty_file(p: Path, label: str):
    if not p.exists():
        raise FileNotFoundError(f"Missing {label}: {p}")
    if p.stat().st_size == 0:
        raise ValueError(f"Empty file {label}: {p}")


def run_and_time(cmd, stdout_path=None, stderr_path=None):
    t0 = time.perf_counter()

    stdout_handle = open(stdout_path, "w", encoding="utf-8") if stdout_path else subprocess.DEVNULL
    stderr_handle = open(stderr_path, "w", encoding="utf-8") if stderr_path else subprocess.DEVNULL

    try:
        subprocess.run(cmd, check=True, stdout=stdout_handle, stderr=stderr_handle, text=True)
    except subprocess.CalledProcessError as e:
        if stdout_path and hasattr(stdout_handle, "close"):
            stdout_handle.close()
        if stderr_path and hasattr(stderr_handle, "close"):
            stderr_handle.close()
            with open(stderr_path, "r", encoding="utf-8", errors="replace") as f:
                err = f.read()
            raise RuntimeError(f"Command failed:\n{cmd}\n\nstderr:\n{err}") from e
        raise
    finally:
        if stdout_path and hasattr(stdout_handle, "close") and not stdout_handle.closed:
            stdout_handle.close()
        if stderr_path and hasattr(stderr_handle, "close") and not stderr_handle.closed:
            stderr_handle.close()

    return time.perf_counter() - t0


def load_blast_hits(path: Path):
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


def summarize_hits(hits_dict):
    rows = []
    for q in sorted(hits_dict.keys(), key=lambda x: int(x.replace("query", ""))):
        subjects = hits_dict[q]
        strong = sum(1 for s in subjects if prot_group(s) == "strong_db")
        weak = sum(1 for s in subjects if prot_group(s) == "weak_db")
        random = sum(1 for s in subjects if prot_group(s) == "random_db")
        motif = sum(1 for s in subjects if prot_group(s) == "motif_db")
        top1 = subjects[0] if subjects else ""

        rows.append({
            "query": q,
            "query_group": query_group(q),
            "n_hits": len(subjects),
            "top1": top1,
            "strong_db_hits": strong,
            "weak_db_hits": weak,
            "random_db_hits": random,
            "motif_db_hits": motif
        })
    return rows


def aggregate_by_query_group(rows):
    groups = defaultdict(lambda: {
        "queries": 0,
        "total_hits": 0,
        "strong_db_hits": 0,
        "weak_db_hits": 0,
        "random_db_hits": 0,
        "motif_db_hits": 0,
    })

    for r in rows:
        g = r["query_group"]
        groups[g]["queries"] += 1
        groups[g]["total_hits"] += r["n_hits"]
        groups[g]["strong_db_hits"] += r["strong_db_hits"]
        groups[g]["weak_db_hits"] += r["weak_db_hits"]
        groups[g]["random_db_hits"] += r["random_db_hits"]
        groups[g]["motif_db_hits"] += r["motif_db_hits"]

    out = []
    for g, vals in groups.items():
        out.append({
            "query_group": g,
            "queries": vals["queries"],
            "avg_total_hits": vals["total_hits"] / vals["queries"],
            "avg_strong_db_hits": vals["strong_db_hits"] / vals["queries"],
            "avg_weak_db_hits": vals["weak_db_hits"] / vals["queries"],
            "avg_random_db_hits": vals["random_db_hits"] / vals["queries"],
            "avg_motif_db_hits": vals["motif_db_hits"] / vals["queries"],
        })
    return out


def save_csv(path: Path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    BLASTDB_DIR.mkdir(parents=True, exist_ok=True)

    ensure_nonempty_file(DB_FASTA, "DB_FASTA")
    ensure_nonempty_file(QUERIES_FASTA, "QUERIES_FASTA")

    makeblastdb = resolve_windows_exe("makeblastdb", BLAST_BIN_DIR)
    blastp = resolve_windows_exe("blastp", BLAST_BIN_DIR)

    # Build DB once
    cmd = [
        makeblastdb,
        "-in", str(DB_FASTA),
        "-dbtype", "prot",
        "-out", str(BLASTDB_DIR / "db")
    ]
    dt_db = run_and_time(
        cmd,
        stdout_path=RESULTS_DIR / "makeblastdb.stdout.txt",
        stderr_path=RESULTS_DIR / "makeblastdb.stderr.txt"
    )

    runtimes = [{"step": "makeblastdb", "seconds": f"{dt_db:.6f}"}]

    blast_outputs = {}
    for ws in [2, 7]:
        out_file = RESULTS_DIR / f"blastp_wordsize_{ws}.tsv"

        cmd = [
            blastp,
            "-query", str(QUERIES_FASTA),
            "-db", str(BLASTDB_DIR / "db"),
            "-word_size", str(ws),
            "-out", str(out_file),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ]

        dt = run_and_time(
            cmd,
            stdout_path=RESULTS_DIR / f"blastp_wordsize_{ws}.stdout.txt",
            stderr_path=RESULTS_DIR / f"blastp_wordsize_{ws}.stderr.txt"
        )
        runtimes.append({"step": f"blastp_wordsize_{ws}", "seconds": f"{dt:.6f}"})
        blast_outputs[ws] = out_file

    save_csv(RESULTS_DIR / "e_runtimes.csv", runtimes)

    for ws, path in blast_outputs.items():
        hits = load_blast_hits(path)
        per_query = summarize_hits(hits)
        by_group = aggregate_by_query_group(per_query)

        save_csv(RESULTS_DIR / f"e_wordsize_{ws}_per_query.csv", per_query)
        save_csv(RESULTS_DIR / f"e_wordsize_{ws}_by_query_group.csv", by_group)

    print("Finished successfully.")
    print("Results written in:", RESULTS_DIR)
    print("Inspect these files:")
    print(" - e_runtimes.csv")
    print(" - e_wordsize_2_per_query.csv")
    print(" - e_wordsize_2_by_query_group.csv")
    print(" - e_wordsize_7_per_query.csv")
    print(" - e_wordsize_7_by_query_group.csv")


if __name__ == "__main__":
    main()