import csv

read_type = input("Input read type (ecoli or lambda): ")

mutated_file = f"{read_type}_mutated.csv"
calls_file = f"{read_type}_calls.csv"
reference_file = f"{read_type}.fasta"

FLANK = 10


def load_reference(path):
    seq = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                seq.append(line.upper())
    return "".join(seq)


def load_calls(path):
    calls = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 3:
                continue

            mut_type = row[0].strip().upper()
            pos = int(row[1].strip())
            base = row[2].strip().upper()

            calls.append((mut_type, pos, base))

    return calls


def possible_effects(call, reference, start, end):

    mut_type, pos, base = call
    effects = []

    if mut_type == "X":
        indices = [pos - 1]

    elif mut_type == "D":
        indices = [pos - 1, pos]

    elif mut_type == "I":
        indices = [pos - 1, pos]

    else:
        return set()

    window = reference[start:end]

    for idx in indices:
        if idx < start or idx > end:
            continue

        rel = idx - start

        if mut_type == "X":
            if 0 <= rel < len(window):
                effects.append(window[:rel] + base + window[rel + 1:])

        elif mut_type == "D":
            if 0 <= rel < len(window):
                effects.append(window[:rel] + window[rel + 1:])

        elif mut_type == "I":
            if 0 <= rel <= len(window):
                effects.append(window[:rel] + base + window[rel:])

    return set(effects)

def equivalent(call1, call2, reference):
    type1, pos1, base1 = call1
    type2, pos2, base2 = call2

    if type1 != type2:
        return False

    if abs(pos1 - pos2) > FLANK + 2:
        return False

    if type1 in {"X", "I"} and base1 != base2:
        return False

    start = max(0, min(pos1, pos2) - FLANK - 1)
    end = min(len(reference), max(pos1, pos2) + FLANK)

    effects1 = possible_effects(call1, reference, start, end)
    effects2 = possible_effects(call2, reference, start, end)

    return bool(effects1 & effects2)


def match_calls(ground_truth, found, reference):
    matched_gt = set()
    matched_found = set()

    for i, found_call in enumerate(found):
        for j, gt_call in enumerate(ground_truth):
            if j in matched_gt:
                continue

            if equivalent(found_call, gt_call, reference):
                matched_found.add(i)
                matched_gt.add(j)
                break

    tp = [found[i] for i in matched_found]
    fp = [found[i] for i in range(len(found)) if i not in matched_found]
    fn = [ground_truth[j] for j in range(len(ground_truth)) if j not in matched_gt]

    return tp, fp, fn


reference = load_reference(reference_file)
ground_truth = load_calls(mutated_file)
found = load_calls(calls_file)

tp, fp, fn = match_calls(ground_truth, found, reference)

print("-------")
print(f"Ground truth:       {len(ground_truth)}")
print(f"Found:              {len(found)}")
print(f"True positives:     {len(tp)}")
print(f"False positives:    {len(fp)}")
print(f"False negatives:    {len(fn)}")

precision = len(tp) / len(found) if found else 0
recall = len(tp) / len(ground_truth) if ground_truth else 0
f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0

print(f"Precision:          {precision:.4f}")
print(f"Recall:             {recall:.4f}")
print(f"F1:                 {f1:.4f}")

# print("\nBy mutation type")
# print("----------------")

""" for mut_type, label in [
    ("X", "substitutions"),
    ("D", "deletions"),
    ("I", "insertions"),
]:
    gt_n = sum(1 for c in ground_truth if c[0] == mut_type)
    found_n = sum(1 for c in found if c[0] == mut_type)
    tp_n = sum(1 for c in tp if c[0] == mut_type)
    fp_n = sum(1 for c in fp if c[0] == mut_type)
    fn_n = sum(1 for c in fn if c[0] == mut_type)

    print(f"{label}:")
    print(f"  My calls:           {found_n}")
    print(f"  ground truth:    {gt_n}")
    print(f"  true positives:  {tp_n}")
    print(f"  false positives: {fp_n}")
    print(f"  false negatives: {fn_n}")
 """

""" print("\nFalse negatives:")
for call in sorted(fn, key=lambda x: x[1]):
    print(call)


print("\nFalse positives:")
for call in sorted(fp, key=lambda x: x[1]):
    print(call) """