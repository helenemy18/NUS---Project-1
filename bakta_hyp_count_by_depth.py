import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append("/home/svu/e0915697/4_neel_swip_project/2_Helene/2_code/scripts/first_scripts")
from grouping_test import build_combined_df
from sample_grouping import get_sample_groups

#test_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/" # à compléter avec le répertoire qui contient les résultats
test_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/virome_project/genomad/genomad_all/gen_subset"
#sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa"
sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/gen_test"
#sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/test_subset"

A_group, B_group, C_group = get_sample_groups(sample_dir)
print(f"A_group = {A_group},\nB_group = {B_group},\nC_group = {C_group}")

output_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/working_dir_H/bact_arch_project/results/plots/all_plots"
os.makedirs(output_dir, exist_ok=True)

groups = {
    "A": A_group,
    "B": B_group,
    "C": C_group
}

## build a combined_df per group and compute hyp/char counts ##
summary_rows = {}
group_dfs = {}

for group_name, sample_ids in groups.items():
    if not sample_ids:
        print(f"Skipping group {group_name}: no samples found")
        continue

    group_df = build_combined_df(sample_ids, test_dir)
    group_dfs[group_name] = group_df

    total_count = len(group_df)
    hyp_count = group_df["Product"].str.contains("hypothetical", case=False, na=False).sum()

    summary_rows[group_name] = {"total": total_count, "hypothetical": hyp_count}

summary = pd.DataFrame.from_dict(summary_rows, orient="index")
summary["hyp_perc"] = summary["hypothetical"] / summary["total"] * 100
print(summary)

## pie chart per group ##
n_groups = len(summary)
fig, axes = plt.subplots(1, n_groups, figsize=(5 * n_groups, 5))

if n_groups == 1:
    axes = [axes]

for ax, (group_name, row) in zip(axes, summary.iterrows()):
    char_count = row["total"] - row["hypothetical"]
    ax.pie(
        [row["hypothetical"], char_count],
        labels=["Hypothetical", "Characterised"],
        autopct="%1.1f%%",
        colors=["#d95f5f", "#5f9ed9"],
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5}
    )
    ax.set_title(f"Group {group_name}")

plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hyp_vs_char_pie_per_group.png"), dpi=300, bbox_inches="tight")
print("Saved to", os.path.join(output_dir, "hyp_vs_char_pie_per_group.png"))