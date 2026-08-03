#!/usr/bin/env python3
"""
bakta_functional_composition_vs_depth.py

Plot #1 from the report figure list: functional category composition of
Bakta CDS annotations, stacked by depth group.

Built directly on:
  - parsing.py               (load_gff3, parse_attributes, expand_attributes
                               -- reused verbatim below)
  - genomad_parsing.py's build_combined_virus_df -- the looping/tagging
                               pattern (loop over sample_id, read a per-
                               sample file, tag with 'Sample Id', concat) is
                               mirrored below as build_combined_bakta_df,
                               just pointed at Bakta's per-sample GFF3
                               instead of geNomad's virus summary TSV.

FUNCTIONAL CATEGORIES
----------------------
Genes are assigned to a broad functional category by keyword-matching their
Bakta `product` annotation against FUNCTION_KEYWORDS below. This is a
coarse, keyword-based heuristic -- not a formal KEGG/pathway mapping -- but
it's enough to show a first-pass depth trend for a report figure. State it
as a limitation: mis-annotated or ambiguously-worded products can be
miscategorized, and anything with no keyword hit (including "hypothetical
protein") falls into "Other / Unclassified". Edit FUNCTION_KEYWORDS freely
to add/refine categories relevant to your samples.

DEPTH GROUPING
---------------
Depth is supplied via a CSV (same pattern as the earlier group-based
script), since I don't know your sample-id naming convention:

    sample_id,depth_group
    BP_A_R1_224,0-10cm
    BP_A_R2_224,10-30cm
    ...

Depth groups are auto-ordered left-to-right on the plot by the first number
found in the label (e.g. "0-10cm" < "10-30cm" < "30-60cm"). Override with
--depth-order if your labels don't sort that way.

OUTPUT
------
  functional_category_counts.csv       (raw counts, depth group x category)
  functional_category_proportions.csv  (same, normalized to fractions)
  functional_composition_by_depth.png  (stacked bar plot)

USAGE
-----
python3 bakta_functional_composition_vs_depth.py \\
    --bakta-dir /path/to/prokka_results \\
    --depth-csv depth_groups.csv \\
    --outdir ./results \\
    --bakta-pattern "{sample}/{sample}.gff"
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")  # non-interactive backend, needed for headless clusters
import matplotlib.pyplot as plt
import pandas as pd


# ==========================================================================
# --- Bakta parsing: reused from parsing.py ---------------------------------
# ==========================================================================

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


# ==========================================================================
# --- Combine samples: mirrors build_combined_virus_df's loop/tag pattern --
# ==========================================================================

def build_combined_bakta_df(sample_ids, bakta_dir, pattern="{sample}/{sample}.gff"):
    """
    Same structure as build_combined_virus_df in genomad_parsing.py:
    loop over sample ids, read one per-sample file, tag rows with
    'Sample Id', concat. Here the per-sample file is Bakta's GFF3 instead
    of geNomad's virus summary TSV.
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
        raise ValueError("No dataframes were built — check your file paths.")

    combined_df = pd.concat(all_dfs, ignore_index=True, sort=False)
    return combined_df


# ==========================================================================
# --- Functional categorization ---------------------------------------------
# ==========================================================================

FUNCTION_KEYWORDS = {
    "Methanogenesis": [
        "methyl-coenzyme m reductase", "methanogenesis", "tetrahydromethanopterin",
        "coenzyme f420", "methanol--corrinoid", "formylmethanofuran",
    ],
    "Sulfate/sulfur cycling": [
        "sulfate reductase", "sulfite reductase", "adenylylsulfate reductase",
        "sulfur oxidation", "thiosulfate",
    ],
    "Nitrogen cycling": [
        "nitrogenase", "nitrate reductase", "nitrite reductase",
        "nitrous oxide reductase", "ammonia monooxygenase",
        "hydroxylamine oxidoreductase", "nitric oxide reductase",
    ],
    "Fermentation": [
        "pyruvate formate-lyase", "alcohol dehydrogenase", "lactate dehydrogenase",
        "acetate kinase", "phosphotransacetylase", "formate dehydrogenase",
    ],
    "Aerobic respiration": [
        "cytochrome c oxidase", "cytochrome o ubiquinol oxidase",
        "cytochrome bd", "nadh dehydrogenase", "succinate dehydrogenase",
        "ubiquinol oxidase",
    ],
    "Carbon degradation (CAZymes)": [
        "glycoside hydrolase", "cellulase", "hemicellulase", "chitinase",
        "xylanase", "glucosidase", "carbohydrate-active", "endoglucanase",
        "beta-glucosidase", "pectate lyase",
    ],
    "Stress/oxidative response": [
        "superoxide dismutase", "catalase", "peroxidase", "heat shock",
        "cold shock", "universal stress",
    ],
}

CATEGORY_ORDER = list(FUNCTION_KEYWORDS.keys()) + ["Other / Unclassified"]


def classify_product(product):
    if not isinstance(product, str) or not product.strip():
        return "Other / Unclassified"
    text = product.lower()
    for category, keywords in FUNCTION_KEYWORDS.items():
        for kw in keywords:
            if kw in text:
                return category
    return "Other / Unclassified"


# ==========================================================================
# --- Depth grouping ----------------------------------------------------------
# ==========================================================================

def load_depth_map(path):
    df = pd.read_csv(path)
    if not {"sample_id", "depth_group"}.issubset(df.columns):
        raise ValueError(f"--depth-csv must have columns sample_id,depth_group, got {df.columns.tolist()}")
    return dict(zip(df["sample_id"].astype(str), df["depth_group"].astype(str)))


def natural_depth_order(depth_labels):
    def sort_key(label):
        m = re.search(r"-?\d+(\.\d+)?", label)
        return float(m.group()) if m else float("inf")
    return sorted(set(depth_labels), key=sort_key)


# ==========================================================================
# --- Main --------------------------------------------------------------------
# ==========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bakta-dir", required=True, help="Directory containing per-sample Bakta output")
    ap.add_argument("--depth-csv", required=True, help="CSV with columns: sample_id,depth_group")
    ap.add_argument("--outdir", default="./functional_composition", help="Output directory")
    ap.add_argument("--bakta-pattern", default="{sample}/{sample}.gff",
                     help="Path template (relative to --bakta-dir) for each sample's Bakta GFF3")
    ap.add_argument("--depth-order", default=None,
                     help="Comma-separated explicit depth group order, overrides auto-sort")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    depth_map = load_depth_map(args.depth_csv)
    sample_ids = list(depth_map.keys())

    combined = build_combined_bakta_df(sample_ids, args.bakta_dir, args.bakta_pattern)
    cds_df = combined[combined["type"].str.lower() == "cds"].copy()
    cds_df["Depth Group"] = cds_df["Sample Id"].map(depth_map)

    if "product" not in cds_df.columns:
        raise ValueError("No 'product' attribute found in GFF3 attributes -- cannot classify function.")

    cds_df["Function Category"] = cds_df["product"].apply(classify_product)

    if args.depth_order:
        depth_order = [d.strip() for d in args.depth_order.split(",")]
    else:
        depth_order = natural_depth_order(cds_df["Depth Group"].dropna().unique())

    counts = (
        cds_df.groupby(["Depth Group", "Function Category"])
        .size()
        .unstack(fill_value=0)
    )
    counts = counts.reindex(index=depth_order, fill_value=0)
    for cat in CATEGORY_ORDER:
        if cat not in counts.columns:
            counts[cat] = 0
    counts = counts[CATEGORY_ORDER]

    proportions = counts.div(counts.sum(axis=1).replace(0, pd.NA), axis=0).fillna(0)

    counts.to_csv(os.path.join(args.outdir, "functional_category_counts.csv"))
    proportions.to_csv(os.path.join(args.outdir, "functional_category_proportions.csv"))

    # ---- stacked bar plot ----
    fig, ax = plt.subplots(figsize=(max(6, 1.3 * len(depth_order)), 7))
    bottom = pd.Series(0.0, index=proportions.index)
    cmap = plt.cm.tab20
    colors = [cmap(i / max(1, len(CATEGORY_ORDER) - 1)) for i in range(len(CATEGORY_ORDER))]

    for cat, color in zip(CATEGORY_ORDER, colors):
        ax.bar(proportions.index, proportions[cat], bottom=bottom, label=cat,
               color=color, edgecolor="black", linewidth=0.3)
        bottom = bottom + proportions[cat]

    ax.set_ylabel("Fraction of CDS genes")
    ax.set_xlabel("Depth group")
    ax.set_title("Functional category composition of CDS genes by depth")
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()

    out_path = os.path.join(args.outdir, "functional_composition_by_depth.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved to", out_path)

    print("\nSamples per depth group:")
    print(cds_df.groupby("Depth Group")["Sample Id"].nunique().reindex(depth_order))
    print("\nCDS counts per depth group:")
    print(counts.sum(axis=1))


if __name__ == "__main__":
    main()