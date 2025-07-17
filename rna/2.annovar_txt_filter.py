
import csv
import argparse
import os

def filter_annovar_v3(input_file, col1_index, col2_index, threshold):
    base_dir = os.path.dirname(input_file)
    output_coding = os.path.join(base_dir, "coding.txt")
    output_non_coding = os.path.join(base_dir, "noncoding-pep.txt")
    coding_framstop = os.path.join(base_dir, "coding_framstop.txt")
    coding_nonfrasnv = os.path.join(base_dir, "coding_snvnonfra.txt")

    with open(input_file, 'r') as infile,  open(output_coding, 'w', newline='') as coding_file, open(output_non_coding, 'w', newline='') as non_coding_file:

        reader = csv.reader(infile, delimiter='\t')
        coding_writer = csv.writer(coding_file, delimiter='\t')
        non_coding_writer = csv.writer(non_coding_file, delimiter='\t')

        header = next(reader)
        coding_writer.writerow(header)
        non_coding_writer.writerow(header)

        otherinfo10_idx = header.index('Otherinfo10')

        for row in reader:
            try:
                value1 = float(row[col1_index]) if row[col1_index] != "." else 0
                value2 = float(row[col2_index]) if row[col2_index] != "." else 0
            except ValueError:
                continue

            if row[otherinfo10_idx] != 'PASS':
                continue

            if value1 <= threshold and value2 <= threshold:
                if row[5] == 'exonic':
                    if row[8] != 'synonymous SNV':
                        coding_writer.writerow(row)
                else:
                    non_coding_writer.writerow(row)

    # 对编码区结果进行再次分类
    with open(output_coding, 'r') as coding_file,open(coding_framstop, 'w', newline='') as framstop_file,open(coding_nonfrasnv, 'w', newline='') as nonfrasnv_file:

        reader = csv.reader(coding_file, delimiter='\t')
        framstop_writer = csv.writer(framstop_file, delimiter='\t')
        nonfrasnv_writer = csv.writer(nonfrasnv_file, delimiter='\t')

        header = next(reader)
        framstop_writer.writerow(header)
        nonfrasnv_writer.writerow(header)

        for row in reader:
            if row[8] in ['stoploss', 'frameshift deletion', 'frameshift insertion']:
                framstop_writer.writerow(row)
            else:
                nonfrasnv_writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Filter ANNOVAR output with extended logic.")
    parser.add_argument("input_file", help="Path to the input ANNOVAR file.")
    parser.add_argument("col1", type=int, help="Column number 1 (1-based, maps to 11-24 in txt).")
    parser.add_argument("col2", type=int, help="Column number 2 (1-based, maps to 11-24 in txt).")
    parser.add_argument("threshold", type=float, help="Threshold for filtering.")

    args = parser.parse_args()

    # Convert 1-based input to 0-based Python indexing, and map to actual column index (i.e., +10)
    col1_index = args.col1 + 10 - 1
    col2_index = args.col2 + 10 - 1

    filter_annovar_v3(args.input_file, col1_index, col2_index, args.threshold)

if __name__ == "__main__":
    main()
