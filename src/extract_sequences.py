import pandas as pd
from pyfaidx import Fasta
from pathlib import Path

GENOME_PATH = Path("data/hg38.fa")
VARIANTS_PATH = Path("outputs/clean_variants.csv")
OUT_PATH = Path("outputs/sequences.csv")

WINDOW = 256  # 256bp each side = 512bp total

df = pd.read_csv(VARIANTS_PATH)
genome = Fasta(GENOME_PATH)
print(f"genome loaded — {len(genome.keys())} chromosomes")
print(f"variants to process: {len(df):,}")

# extract sequence window around each variant
results = []
skipped = 0

for i, row in df.iterrows():
    chrom = row["chrom"]
    pos   = row["pos"] - 1      # pyfaidx is 0-based, VCF is 1-based
    start = pos - WINDOW
    end   = pos + WINDOW

    # skip if window goes out of chromosome bounds
    if chrom not in genome:
        skipped += 1
        continue
    chrom_len = len(genome[chrom])
    if start < 0 or end > chrom_len:
        skipped += 1
        continue

    seq = genome[chrom][start:end].seq.upper()

    # skip if sequence has too many unknown bases
    if seq.count("N") > 0.1 * len(seq):
        skipped += 1
        continue

    results.append({
        "chrom": chrom,
        "pos":   row["pos"],
        "ref":   row["ref"],
        "alt":   row["alt"],
        "label": row["label"],
        "seq":   seq,
    })

    if (i + 1) % 10000 == 0:
        print(f"  processed {i+1:,} / {len(df):,} — skipped so far: {skipped}")

out = pd.DataFrame(results)
out.to_csv(OUT_PATH, index=False)

print(f"\ndone")
print(f"saved: {len(out):,} sequences → {OUT_PATH}")
print(f"skipped: {skipped:,}")
print(out.head())