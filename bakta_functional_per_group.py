import matplotlib
matplotlib.use("Agg")  # non-interactive backend, needed for headless clusters — must be set before pyplot import
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append("/home/svu/e0915697/4_neel_swip_project/2_Helene/2_code/scripts/first_scripts")

from grouping_test import build_combined_df
from sample_grouping import get_sample_groups

bakta_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/annotated_bact_all"
sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa"
#sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/test_subset"


A_group, B_group, C_group = get_sample_groups(sample_dir)
print(f"A_group = {A_group},\nB_group = {B_group},\nC_group = {C_group}")

#output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/plots"
output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/plots/all_plots"
os.makedirs(output_dir, exist_ok=True)

groups = {
    "A": A_group,
    "B": B_group,
    "C": C_group
}

N = 20  # top N products per group

for group_name, sample_ids in groups.items():
    if not sample_ids:
        print(f"Skipping group {group_name}: no samples found")
        continue

    group_df = build_combined_df(sample_ids, bakta_dir)

    # Exclude hypothetical proteins and missing products
    annotated_df = group_df[
        group_df["Product"].notna() &
        ~group_df["Product"].str.contains("hypothetical", case=False, na=False)
    ].copy()

    total_annotated = len(annotated_df)

    top_products = annotated_df["Product"].value_counts().head(N).sort_values()
    top_products_pct = top_products / total_annotated * 100

    fig, ax = plt.subplots(figsize=(11, 8))
    bars = ax.barh(top_products.index, top_products.values, color="darkorange", edgecolor="black")

    ax.set_title(f"Group {group_name}: Top {N} Functional Products (Excluding Hypothetical Proteins)")
    ax.set_xlabel("Count")
    ax.set_ylabel("Product")

    # Add percentage labels at the end of each bar
    for bar, pct in zip(bars, top_products_pct):
        width = bar.get_width()
        ax.text(
            width + (top_products.max() * 0.01),
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.1f}%",
            va="center",
            ha="left",
            fontsize=9
        )

    ax.set_xlim(0, top_products.max() * 1.15)

    plt.tight_layout()

    out_path = os.path.join(output_dir, f"top_products_group_{group_name}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)  # free memory before the next group's plot
    print("Saved to", out_path)


