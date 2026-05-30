import csv

read_type = input("Input read type (ecoli or lambda): ")

my_file = read_type + "_calls.csv"
bcf_file = read_type + "_bcftools.vcf"


def load_my_calls(path):
    calls = set()

    with open(path, newline="", encoding="utf-8") as file:
        reader = csv.reader(file)

        for row in reader:
            calls.add((row[0], int(row[1]), row[2]))

    return calls


def load_bcftools_calls(path):
    calls = set()

    with open(path, encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            pos = int(fields[1]) - 1
            ref = fields[3]
            alt = fields[4]

            if "," in alt:
                continue

            if len(ref) == 1 and len(alt) == 1:
                calls.add(("X", pos, alt))

            elif len(ref) < len(alt):
                inserted = alt[1:]
                calls.add(("I", pos + 1, inserted))

            elif len(ref) > len(alt):
                calls.add(("D", pos + 1, "-"))

    return calls


found = load_my_calls(my_file)
ground_truth = load_bcftools_calls(bcf_file)

matched_mine = set()
matched_bcftools = set()

for i, my_call in enumerate(found):
    for j, bcf_call in enumerate(ground_truth):
        if j in matched_bcftools:
            continue

        if my_call == bcf_call:
            matched_mine.add(i)
            matched_bcftools.add(j)
            break

tp = matched_mine
fp = set(range(len(found))) - matched_mine
fn = set(range(len(ground_truth))) - matched_bcftools

precision = len(tp) / len(found) if found else 0
recall = len(tp) / len(ground_truth) if ground_truth else 0

if precision + recall > 0:
    f1 = 2 * precision * recall / (precision + recall)
else:
    f1 = 0

print(f"Ground truth:       {len(ground_truth)}")
print(f"Found:              {len(found)}")
print(f"True positives:     {len(tp)}")
print(f"False positives:    {len(fp)}")
print(f"False negatives:    {len(fn)}")
print(f"Precision:          {precision:.4f}")
print(f"Recall:             {recall:.4f}")
print(f"F1:                 {f1:.4f}")