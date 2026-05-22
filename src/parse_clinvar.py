import gzip
import pandas as pd
from pathlib import Path

VCF_PATH = Path("data/clinvar.vcf.gz")
OUT_PATH = Path("outputs/clean_variants.csv")
OUT_PATH.parent.mkdir(exist_ok=True)

def parse_info(info_str, key):
    for field in info_str.split(";"):
        if field.startswith(f"{key}="):
            return field[len(key)+1:]
    return None

# parse the vcf
rows = []
with gzip.open(VCF_PATH, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 8:
            continue

        chrom, pos, _, ref, alt, _, _, info = parts[:8]

        if parse_info(info, "CLNVC") != "single_nucleotide_variant":
            continue

        rows.append({
            "chrom":  f"chr{chrom}",
            "pos":    int(pos),
            "ref":    ref,
            "alt":    alt,
            "clnsig": parse_info(info, "CLNSIG"),
        })

df = pd.DataFrame(rows)
print(f"SNVs loaded: {len(df):,}")
print(df["clnsig"].value_counts().head(15))

# keep only confident labels
label_map = {"Pathogenic": 1, "Benign": 0}
df = df[df["clnsig"].isin(label_map)].copy()
df["label"] = df["clnsig"].map(label_map)
df = df.drop(columns=["clnsig"]).reset_index(drop=True)

print(f"\nPathogenic/Benign only: {len(df):,}")
print(df["label"].value_counts())

# balance classes
n = min(len(df[df["label"] == 0]), len(df[df["label"] == 1]))
df = pd.concat([
    df[df["label"] == 0].sample(n, random_state=42),
    df[df["label"] == 1].sample(n, random_state=42),
]).sample(frac=1, random_state=42).reset_index(drop=True)

print(f"\nBalanced: {len(df):,} variants")
print(df["label"].value_counts())

df.to_csv(OUT_PATH, index=False)
print(f"\nsaved → {OUT_PATH}")
print(df.head())