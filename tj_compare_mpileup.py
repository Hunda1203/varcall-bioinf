import csv
from collections import Counter

read_type = input("Input read type (ecoli or lambda): ")

my_file = read_type + "_calls.csv"
mpileup_file = read_type + "_mpileup.txt"

min_reads = 3
min_share = 0.5


def load_my_calls(path):
    calls = set()

    with open(path, newline="", encoding="utf-8") as file:
        reader = csv.reader(file)

        for row in reader:
            if len(row) == 3:
                calls.add((row[0], int(row[1]), row[2]))

    return calls


def read_mpileup_calls(path):
    calls = set()

    with open(path, encoding="utf-8") as file:
        for line in file:
            fields = line.strip().split("\t")

            if len(fields) < 5:
                continue

            pos = int(fields[1]) - 1
            depth = int(fields[3])
            bases = fields[4]

            if depth == 0:
                continue

            subs = Counter()
            inserts = Counter()
            dels = Counter()

            i = 0

            while i < len(bases):
                char = bases[i]

                if char == "^":
                    i += 2

                elif char == "$":
                    i += 1

                elif char in ".,":
                    i += 1

                elif char in "ACGTacgt":
                    subs[char.upper()] += 1
                    i += 1

                elif char == "+" or char == "-":
                    sign = char
                    i += 1

                    number = ""

                    while i < len(bases) and bases[i].isdigit():
                        number += bases[i]
                        i += 1

                    length = int(number)

                    seq = bases[i:i+length].upper()

                    if sign == "+":
                        inserts[seq] += 1
                    else:
                        dels[seq] += 1

                    i += length

                elif char == "*":
                    i += 1

                else:
                    i += 1

            if subs:
                mut = list(subs.keys())[0]

                for base in subs:
                    if subs[base] > subs[mut]:
                        mut = base

                if subs[mut] >= min_reads and subs[mut] / depth >= min_share:
                    calls.add(("X", pos, mut))

            if inserts:
                inserted = list(inserts.keys())[0]

                for seq in inserts:
                    if inserts[seq] > inserts[inserted]:
                        inserted = seq

                if inserts[inserted] >= min_reads and inserts[inserted] / depth >= min_share:
                    calls.add(("I", pos + 1, inserted))

            if dels:
                deleted = list(dels.keys())[0]

                for seq in dels:
                    if dels[seq] > dels[deleted]:
                        deleted = seq

                if len(deleted) == 1 and dels[deleted] >= min_reads and dels[deleted] / depth >= min_share:
                    calls.add(("D", pos + 1, "-"))

    return calls


my_calls = load_my_calls(my_file)
mpileup_calls = read_mpileup_calls(mpileup_file)

true_positives = my_calls & mpileup_calls
false_positives = my_calls - mpileup_calls
false_negatives = mpileup_calls - my_calls

precision = len(true_positives) / len(my_calls) if my_calls else 0
recall = len(true_positives) / len(mpileup_calls) if mpileup_calls else 0

if precision + recall > 0:
    f1 = 2 * precision * recall / (precision + recall)
else:
    f1 = 0


print("-------")
print(f"Ground truth:       {len(mpileup_calls)}")
print(f"Found:              {len(my_calls)}")
print(f"True positives:     {len(true_positives)}")
print(f"False positives:    {len(false_positives)}")
print(f"False negatives:    {len(false_negatives)}")
print(f"Precision:          {precision:.4f}")
print(f"Recall:             {recall:.4f}")
print(f"F1:                 {f1:.4f}")