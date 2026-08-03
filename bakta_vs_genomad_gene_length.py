import matplotlib
matplotlib.use("Agg")  # non-interactive backend, needed for headless clusters — must be set before pyplot import
import matplotlib.pyplot as plt
import pandas as pd
import os

from genomad_parsing import build_combined_virus_df
from sample_grouping import get_sample_groups

## Test run
#test_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/virome_project/genomad/genomad_all/gen_subset"
#sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/gen_test"

## Full run
genomad_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/virome_project/genomad/genomad_all/gen_results"
bakta_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/annotated_bact_all"
sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa"

A_group, B_group, C_group = get_sample_groups(sample_dir)
print(f"A_group = {A_group},\nB_group = {B_group},\nC_group = {C_group}")

# This figure is NOT depth-stratified -- it pools every sample across all
# three depth groups to give one overall Bakta-vs-geNomad comparison.
all_samples = A_group + B_group + C_group
print(f"\nTotal samples pooled: {len(all_samples)}")

output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/virome_project/genomad/genomad_all/gen_plots"
os.makedirs(output_dir, exist_ok=True)

# Adjust if your Bakta output uses a different per-sample filename/extension
# (e.g. "{sample}/{sample}.gff" instead of ".gff3").
BAKTA_PATTERN = "{sample}/{sample}.gff3"


# ---------------------------------------------------------------------------
# Bakta GFF3 parsing -- reused from parsing.py
# ---------------------------------------------------------------------------

def load_gff3(filepath):
    """Load a GFF3 file into a pandas DataFrame, skipping headers/comments."""
    cols = ["seqid", "source", "type", "start", "end",
            "score", "strand", "phase", "attributes"]

    df = pd.read_csv(
        filepath,
        sep="\t",
        comment="#",
        names=cols,
        header=None,
        dtype={"start": "Int64", "end": "Int64"},
        engine="python"
    )
    df = df.dropna(subset=["type"])
    return df


def parse_attributes(attr_string):
    """Turn 'ID=xxx;Name=yyy;product=zzz' into a dict."""
    attrs = {}
    for field in attr_string.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def expand_attributes(df):
    """Expand the attributes column into separate columns."""
    attr_dicts = df["attributes"].apply(parse_attributes)
    attr_df = pd.DataFrame(attr_dicts.tolist(), index=df.index)
    return pd.concat([df.drop(columns=["attributes"]), attr_df], axis=1)


def build_combined_bakta_df(sample_ids, bakta_dir, pattern=BAKTA_PATTERN):
    """
    Same loop/tag pattern as build_combined_virus_df: read one per-sample
    GFF3, tag with 'Sample Id', concat across samples.
    """
    all_dfs = []
    for sample_id in sample_ids:
        gff_file = os.path.join(bakta_dir, pattern.format(sample=sample_id))
        if not os.path.exists(gff_file):
            print(f"Missing file: {gff_file}")
            continue
        df = expand_attributes(load_gff3(gff_file))
        if df.empty:
            print(f"No features found for {sample_id} (empty file)")
            continue
        df.insert(0, "Sample Id", sample_id)
        all_dfs.append(df)
    if not all_dfs:
        raise ValueError("No dataframes were built — check your Bakta file paths.")
    return pd.concat(all_dfs, ignore_index=True, sort=False)


# ---------------------------------------------------------------------------
# geNomad virus gene parsing -- same pattern as build_combined_virus_df,
# but reads *_contigs_virus_genes.tsv (per-gene) instead of
# *_contigs_virus_summary.tsv (per-contig).
# ---------------------------------------------------------------------------

def build_combined_virus_genes_df(list_samples, genomad_dir):
    all_dfs = []
    for sample_id in list_samples:
        summary_dir = f"{genomad_dir}/{sample_id}/{sample_id}_contigs_summary"
        tsv_file = f"{summary_dir}/{sample_id}_contigs_virus_genes.tsv"
        if not os.path.exists(tsv_file):
            print(f"Missing file: {tsv_file}")
            continue
        df = pd.read_csv(tsv_file, sep="\t")
        if df.empty:
            print(f"No viral genes found for {sample_id} (empty file)")
            continue
        df.columns = [c.strip() for c in df.columns]
        df.insert(0, "Sample Id", sample_id)
        all_dfs.append(df)
    if not all_dfs:
        raise ValueError("No dataframes were built — check your geNomad file paths.")
    return pd.concat(all_dfs, ignore_index=True, sort=False)


# ---------------------------------------------------------------------------
# Build the two pooled (non-depth-stratified) gene length distributions
# ---------------------------------------------------------------------------

bakta_combined = build_combined_bakta_df(all_samples, bakta_dir)
bakta_cds = bakta_combined[bakta_combined["type"].str.lower() == "cds"].copy()
bakta_cds["length"] = bakta_cds["end"] - bakta_cds["start"]
bakta_lengths = bakta_cds["length"].dropna()

genomad_genes = build_combined_virus_genes_df(all_samples, genomad_dir)
genomad_genes["length"] = pd.to_numeric(genomad_genes["length"], errors="coerce")
genomad_lengths = genomad_genes["length"].dropna()

print(f"\nBakta (all CDS, all samples):    n={len(bakta_lengths)}, "
      f"mean={bakta_lengths.mean():.1f} bp, median={bakta_lengths.median():.1f} bp")
print(f"geNomad (viral genes, all samples): n={len(genomad_lengths)}, "
      f"mean={genomad_lengths.mean():.1f} bp, median={genomad_lengths.median():.1f} bp")

pd.DataFrame({
    "tool": ["Bakta"] * len(bakta_lengths) + ["geNomad"] * len(genomad_lengths),
    "length": pd.concat([bakta_lengths, genomad_lengths], ignore_index=True),
}).to_csv(os.path.join(output_dir, "bakta_vs_genomad_gene_length_overall.csv"), index=False)


# ---------------------------------------------------------------------------
# Plot: boxplot + overlaid histogram, side by side in one figure
# ---------------------------------------------------------------------------

colors = {"Bakta": "#4C72B0", "geNomad": "#DD8452"}

fig, (ax_box, ax_hist) = plt.subplots(1, 2, figsize=(13, 6))

bp = ax_box.boxplot(
    [bakta_lengths, genomad_lengths],
    tick_labels=[f"Bakta CDS\n(n={len(bakta_lengths)})", f"geNomad viral genes\n(n={len(genomad_lengths)})"],
    showfliers=False,
    patch_artist=True,
)
for patch, tool in zip(bp["boxes"], ["Bakta", "geNomad"]):
    patch.set_facecolor(colors[tool])
    patch.set_alpha(0.7)
ax_box.set_ylabel("Gene length (bp)")
ax_box.set_title("Gene length distribution")

ax_hist.hist(bakta_lengths, bins=50, alpha=0.5, label=f"Bakta CDS (n={len(bakta_lengths)})",
             color=colors["Bakta"], edgecolor="black", linewidth=0.3)
ax_hist.hist(genomad_lengths, bins=50, alpha=0.5, label=f"geNomad viral genes (n={len(genomad_lengths)})",
             color=colors["geNomad"], edgecolor="black", linewidth=0.3)
ax_hist.set_xlabel("Gene length (bp)")
ax_hist.set_ylabel("Count")
ax_hist.set_title("Gene length histogram")
ax_hist.legend()

fig.suptitle("Bakta (bacterial/archaeal) vs. geNomad (viral) gene length — all samples pooled")
plt.tight_layout()

out_path = os.path.join(output_dir, "bakta_vs_genomad_gene_length_overall.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print("\nSaved to", out_path)