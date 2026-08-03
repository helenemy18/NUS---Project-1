import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append("/home/svu/e0915697/4_neel_swip_project/2_Helene/2_code/scripts/first_scripts")

from gc_content import build_gc_df
from sample_grouping import get_sample_groups

fna_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/annotated_bact_all"
#sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/test_subset"
sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa"
#output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/plots"
output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/plots/all_plots"
os.makedirs(output_dir, exist_ok=True)

A_group, B_group, C_group = get_sample_groups(sample_dir)
print(f"A_group = {A_group},\nB_group = {B_group},\nC_group = {C_group}")

groups = {
    "A": A_group,
    "B": B_group,
    "C": C_group
}

# Build a GC dataframe per group, tagging each row with its group letter
all_gc_dfs = []
for group_name, sample_ids in groups.items():
    if not sample_ids:
        print(f"Skipping group {group_name}: no samples found")
        continue
    group_gc_df = build_gc_df(sample_ids, fna_dir)
    group_gc_df["Group"] = group_name
    all_gc_dfs.append(group_gc_df)

gc_df = pd.concat(all_gc_dfs, ignore_index=True)

# --- Overall histogram (all samples combined) ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(gc_df["GC%"], bins=50, color="seagreen", edgecolor="black")
ax.set_title("GC Content Distribution per Contig (All Samples)")
ax.set_xlabel("GC%")
ax.set_ylabel("Number of Contigs")
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "gc_content_distribution.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

# --- By group (A/B/C) ---
fig, ax = plt.subplots(figsize=(9, 6))
for group_name, sub_df in gc_df.groupby("Group"):
    sub_df["GC%"].plot(kind="density", ax=ax, label=group_name)

ax.set_title("GC Content Density by Group")
ax.set_xlabel("GC%")
ax.legend(title="Group")
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "gc_content_by_group.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

print("Plots saved to", output_dir)


