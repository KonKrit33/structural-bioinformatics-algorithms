from pathlib import Path
import subprocess
import shutil
import time
import csv


# =========================
# CONFIG
# =========================

ROOT = Path(r"C:\bio")
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results_b"
BLASTDB_DIR = RESULTS_DIR / "blastdb"

DB_FASTA = DATA_DIR / "db.fasta"
QUERIES_FASTA = DATA_DIR / "queries.fasta"

BLAST_BIN_DIR = ROOT / "data" / "bin"

WSL_FASTA36 = "/home/konkrit/fasta36/bin/fasta36"
WSL_SSEARCH36 = "/home/konkrit/fasta36/bin/ssearch36"
WSL_GGSEARCH36 = "/home/konkrit/fasta36/bin/ggsearch36"


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

    raise FileNotFoundError(f"Δεν βρέθηκε το executable: {name}")


def windows_to_wsl_path(p: Path) -> str:
    p = p.resolve()
    drive = p.drive[0].lower()
    rest = p.as_posix()[2:]
    return f"/mnt/{drive}{rest}"


def ensure_nonempty_file(p: Path, label: str):
    if not p.exists():
        raise FileNotFoundError(f"Δεν βρέθηκε το {label}: {p}")
    if p.stat().st_size == 0:
        raise ValueError(f"Το αρχείο {label} είναι άδειο: {p}")


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
            raise RuntimeError(f"Απέτυχε η εντολή:\n{cmd}\n\nstderr:\n{err}") from e
        raise
    finally:
        if stdout_path and hasattr(stdout_handle, "close") and not stdout_handle.closed:
            stdout_handle.close()
        if stderr_path and hasattr(stderr_handle, "close") and not stderr_handle.closed:
            stderr_handle.close()

    return time.perf_counter() - t0


def run_wsl_and_time(program_path, query_path_win: Path, db_path_win: Path, output_file: Path, stderr_file: Path):
    query_wsl = windows_to_wsl_path(query_path_win)
    db_wsl = windows_to_wsl_path(db_path_win)

    cmd = [
        "wsl",
        program_path,
        "-m", "8CBL",
        query_wsl,
        db_wsl
    ]
    return run_and_time(cmd, stdout_path=output_file, stderr_path=stderr_file)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    BLASTDB_DIR.mkdir(parents=True, exist_ok=True)

    ensure_nonempty_file(DB_FASTA, "DB_FASTA")
    ensure_nonempty_file(QUERIES_FASTA, "QUERIES_FASTA")

    makeblastdb = resolve_windows_exe("makeblastdb", BLAST_BIN_DIR)
    blastp = resolve_windows_exe("blastp", BLAST_BIN_DIR)

    runtimes = []

    cmd = [
        makeblastdb,
        "-in", str(DB_FASTA),
        "-dbtype", "prot",
        "-out", str(BLASTDB_DIR / "db")
    ]
    dt = run_and_time(
        cmd,
        stdout_path=RESULTS_DIR / "makeblastdb.stdout.txt",
        stderr_path=RESULTS_DIR / "makeblastdb.stderr.txt"
    )
    runtimes.append(("makeblastdb", dt))

    cmd = [
        blastp,
        "-query", str(QUERIES_FASTA),
        "-db", str(BLASTDB_DIR / "db"),
        "-out", str(RESULTS_DIR / "blastp.tsv"),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    dt = run_and_time(
        cmd,
        stdout_path=RESULTS_DIR / "blastp.stdout.txt",
        stderr_path=RESULTS_DIR / "blastp.stderr.txt"
    )
    runtimes.append(("blastp", dt))

    dt = run_wsl_and_time(
        WSL_FASTA36,
        QUERIES_FASTA,
        DB_FASTA,
        RESULTS_DIR / "fasta36.tsv",
        RESULTS_DIR / "fasta36.stderr.txt"
    )
    runtimes.append(("fasta36", dt))

    dt = run_wsl_and_time(
        WSL_SSEARCH36,
        QUERIES_FASTA,
        DB_FASTA,
        RESULTS_DIR / "ssearch36.tsv",
        RESULTS_DIR / "ssearch36.stderr.txt"
    )
    runtimes.append(("ssearch36", dt))

    dt = run_wsl_and_time(
        WSL_GGSEARCH36,
        QUERIES_FASTA,
        DB_FASTA,
        RESULTS_DIR / "ggsearch36.tsv",
        RESULTS_DIR / "ggsearch36.stderr.txt"
    )
    runtimes.append(("ggsearch36", dt))

    with open(RESULTS_DIR / "runtimes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["program", "seconds"])
        for name, sec in runtimes:
            w.writerow([name, f"{sec:.6f}"])

    print("Finished successfully.")
    print("Results folder:", RESULTS_DIR)


if __name__ == "__main__":
    main()