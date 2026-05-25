import os
import random

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(BASE_DIR, "..", "data"))

ALPHABET = list(
    "ACDEFGHIKLMNPQRSTVWY"
)

SEQ_LEN = 200

N_HOMOLOGS = 25
N_WEAK_HOMOLOGS = 25
N_RANDOM = 25
N_MOTIF = 25

N_QUERY_STRONG = 3
N_QUERY_WEAK = 3
N_QUERY_RANDOM = 3
N_QUERY_FRAGMENT = 3
N_QUERY_MOTIF = 3

MUT_STRONG = 0.15
MUT_WEAK = 0.50

random.seed(0)


def random_seq(n):
    return "".join(
        random.choice(ALPHABET)
        for _ in range(n)
    )


def mutate(seq, rate):
    s = list(seq)
    for i in range(len(s)):
        if random.random() < rate:
            choices = [
                a for a in ALPHABET
                if a != s[i]
            ]
            s[i] = random.choice(choices)
    return "".join(s)


BASE = random_seq(SEQ_LEN)
MOTIF = BASE[80:100]


def insert_motif(seq):
    pos = random.randint(40, 160)
    return (
        seq[:pos] +
        MOTIF +
        seq[pos + len(MOTIF):]
    )


db = []

for _ in range(N_HOMOLOGS):
    db.append(
        mutate(BASE, MUT_STRONG)
    )

for _ in range(N_WEAK_HOMOLOGS):
    db.append(
        mutate(BASE, MUT_WEAK)
    )

for _ in range(N_RANDOM):
    db.append(
        random_seq(SEQ_LEN)
    )

for _ in range(N_MOTIF):
    s = random_seq(SEQ_LEN)
    db.append(
        insert_motif(s)
    )


queries = []

for _ in range(N_QUERY_STRONG):
    queries.append(
        mutate(BASE, 0.05)
    )

for _ in range(N_QUERY_WEAK):
    queries.append(
        mutate(BASE, 0.50)
    )

for _ in range(N_QUERY_RANDOM):
    queries.append(
        random_seq(SEQ_LEN)
    )

fragment = BASE[80:120]

for _ in range(N_QUERY_FRAGMENT):
    queries.append(
        random_seq(80)
        + fragment +
        random_seq(80)
    )

for _ in range(N_QUERY_MOTIF):
    queries.append(
        random_seq(80)
        + MOTIF +
        random_seq(80)
    )


os.makedirs(OUTDIR, exist_ok=True)

with open(os.path.join(OUTDIR, "db.fasta"), "w") as f:
    for i, seq in enumerate(db):
        f.write(
            f">prot{i}\n{seq}\n"
        )

with open(os.path.join(OUTDIR, "queries.fasta"), "w") as f:
    for i, seq in enumerate(queries):
        f.write(
            f">query{i}\n{seq}\n"
        )

with open(os.path.join(OUTDIR, "reference.fasta"), "w") as f:
    f.write(f">reference\n{BASE}\n")

print("Files written to:", OUTDIR)


print("\n--- FIRST 3 DB SEQUENCES ---")
for i, seq in enumerate(db[:3]):
    print(f">prot{i}")
    print(seq)

print("\n--- ALL QUERIES ---")
for i, seq in enumerate(queries):
    print(f">query{i}")
    print(seq)

print("REFERENCE SEQUENCE:")
print(BASE)