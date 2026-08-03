import csv
from sample_grouping import get_depth_groups

sample_dir = "/home/svu/e0915697/4_neel_swip_project/2_Helene/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa"
output_csv = "depth_groups.csv"


A_group, B_group, C_group = get_depth_groups(sample_dir)

print(f"A_group ({len(A_group)}): {A_group}")
print(f"B_group ({len(B_group)}): {B_group}")
print(f"C_group ({len(C_group)}): {C_group}")
print(f"\nTotal samples grouped: {len(A_group) + len(B_group) + len(C_group)}")

# --- Sanity check: flag any sample in C_group whose second token doesn't
# actually end in "C" -- these got there only via the `else` catch-all, not
# because they're a genuine "C" sample. Worth eyeballing before trusting
# the grouping at scale.
suspicious = [s for s in C_group if s.split("_")[1][-1] != "C"]
if suspicious:
    print(f"\n[WARNING] {len(suspicious)} sample(s) landed in group C via the "
          f"catch-all 'else' branch, NOT because their token ends in 'C'. "
          f"Double-check these are genuinely meant to be group C:")
    for s in suspicious:
        print(f"    {s}  (second token: '{s.split('_')[1]}')")

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "depth_group"])
    for s in A_group:
        writer.writerow([s, "A"])
    for s in B_group:
        writer.writerow([s, "B"])
    for s in C_group:
        writer.writerow([s, "C"])

print(f"\nSaved to {output_csv}")
