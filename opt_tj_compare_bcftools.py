import csv

read_type = "ecoli"

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


my_calls = load_my_calls(my_file)
bcftools_calls = load_bcftools_calls(bcf_file)

same = my_calls & bcftools_calls
only_mine = my_calls - bcftools_calls
only_bcftools = bcftools_calls - my_calls

precision = len(same) / len(my_calls) if my_calls else 0
recall = len(same) / len(bcftools_calls) if bcftools_calls else 0

if precision + recall > 0:
    f1 = 2 * precision * recall / (precision + recall)
else:
    f1 = 0

print("My calls:", len(my_calls))
print("bcftools calls:", len(bcftools_calls))
print("Same:", len(same))
print("Only mine:", len(only_mine))
print("Only bcftools:", len(only_bcftools))
print("Precision:", precision)
print("Recall:", recall)
print("F1:", f1)