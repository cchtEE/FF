from collections import defaultdict
import sys

# Argv 1 - sample name
# Argv 2 - BIN SIZE

patient = sys.argv[1]
bin_size = int(sys.argv[2])

conversion = {"CM000663.2": "chr1", "CM000664.2": "chr2", "CM000665.2": "chr3", "CM000666.2": "chr4",
              "CM000667.2": "chr5", "CM000668.2": "chr6", "CM000669.2": "chr7", "CM000670.2": "chr8",
              "CM000671.2": "chr9", "CM000672.2": "chr10", "CM000673.2": "chr11", "CM000674.2": "chr12",
              "CM000675.2": "chr13", "CM000676.2": "chr14", "CM000677.2": "chr15", "CM000678.2": "chr16",
              "CM000679.2": "chr17", "CM000680.2": "chr18", "CM000681.2": "chr19", "CM000682.2": "chr20",
              "CM000683.2": "chr21", "CM000684.2": "chr22", "CM000685.2": "chrX", "CM000686.2": "chrY"}

bins_per_chromosome = list(
    map(lambda x: (x[0], int(x[1] / bin_size)),
        [["chr1", 249300000],  # chr_name <-> chr length ----> chr_name <-> bin count
         ["chr2", 243200000],
         ["chr3", 198050000],
         ["chr4", 191200000],
         ["chr5", 180950000],
         ["chr6", 171150000],
         ["chr7", 159150000],
         ["chr8", 146400000],
         ["chr9", 141250000],
         ["chr10", 135550000],
         ["chr11", 135050000],
         ["chr12", 133900000],
         ["chr13", 115200000],
         ["chr14", 107350000],
         ["chr15", 102550000],
         ["chr16", 90400000],
         ["chr17", 81200000],
         ["chr18", 78100000],
         ["chr19", 59150000],
         ["chr20", 64444167],
         ["chr21", 48150000],
         ["chr22", 51350000],
         ["chrX", 155300000],
         ["chrY", 59400000]]))


def get_writable_header():
    to_file = "patient"
    for chromosome, bins in bins_per_chromosome:
        for i in range(bins):
            to_file += ";" + (chromosome + "_" + str(i))
    return to_file


def get_writable_line(patient_name, read_counts):
    # Generate line for file
    to_file = patient_name
    for chr_, bins in bins_per_chromosome:
        for i in range(bins):
            to_file += ";" + ("0" if read_counts[chr_][i] is None else str(read_counts[chr_][i]))
    return to_file


read_count_dict = defaultdict(lambda: defaultdict(int))

for line_number, line in enumerate(sys.stdin):
    try:
        splitted = line.split("\t")
        read_chromosome, read_pos = conversion[splitted[2]], int(splitted[3])
        # Add 1 to correct bin in correct chromosome in read counts
        read_count_dict[read_chromosome][(int(read_pos / bin_size))] += 1
    except (KeyError, IndexError) as e:
        continue

sys.stdout.write(get_writable_header() + "\n" + get_writable_line(patient, read_count_dict) + "\n")
