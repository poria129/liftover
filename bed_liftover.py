from pyliftover import LiftOver

# Load liftover converter
lo = LiftOver("hg19ToHg38.over.chain.gz")

# Input and output files
input_file = "converted_output_37.bed"
output_file = "output_hg38.bed"
unmapped_file = "unmapped.bed"

with open(input_file, "r") as fin, open(output_file, "w") as fout, open(
    unmapped_file, "w"
) as funmapped:

    for line in fin:
        if line.startswith("#") or not line.strip():
            continue  # skip comments or empty lines

        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        # UCSC LiftOver expects "chr1", "chr2", ..., "chrX", "chrY", "chrM"
        if chrom == "MT":
            ucsc_chrom = "chrM"
        elif chrom == "XY":
            # XY is non-standard; no such chromosome in UCSC reference genomes
            funmapped.write(line)
            continue
        else:
            ucsc_chrom = f"chr{chrom}"

        # Try to lift over
        lifted = lo.convert_coordinate(ucsc_chrom, start)
        if lifted:
            new_chrom, new_start, strand, _ = lifted[0]
            new_end = new_start + (end - start)

            # Remove "chr" prefix from output as per your request
            new_chrom = new_chrom.replace("chr", "")
            new_fields = [new_chrom, str(int(new_start)), str(int(new_end))] + fields[
                3:
            ]
            fout.write("\t".join(new_fields) + "\n")
        else:
            funmapped.write(line)
