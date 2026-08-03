import os

def get_depth_groups(directory):
    sample_letters = {}

    for filename in os.listdir(directory):
        if filename.endswith("_contigs.fa"):
            sample_id = filename.replace("_contigs.fa", "")
            parts = sample_id.split("_")
            letter = parts[1][-1]
            sample_letters[sample_id] = letter

    A_group, B_group, C_group = [], [], []

    for sample_id, letter in sample_letters.items():
        if letter == "A":
            A_group.append(sample_id)
        elif letter == "B":
            B_group.append(sample_id)
        else:
            C_group.append(sample_id)

    return A_group, B_group, C_group