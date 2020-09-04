#!/usr/bin/env python
import gzip
import pysam
import sys
import argparse

# require that filename contains the specified extension
def require_extension(filename, ext):
    assert(args.input_R1[-len(ext):] == ext)

def filter_reads(filter_sequences, input_R1_fp, input_R2_fp, output_R1_fp, output_R2_fp):
    
    input_R1 = pysam.FastxFile(input_R1_fp)
    input_R2 = pysam.FastxFile(input_R2_fp)
    
    output_R1 = gzip.open(output_R1_fp, 'w')
    output_R2 = gzip.open(output_R2_fp, 'w')

    reads_filtered = 0
    reads_kept = 0
    for R1, R2 in zip(input_R1, input_R2):
        
        discard = False
        for f in filter_sequences:
            if f in R1.sequence or f in R2.sequence:
                discard = True
                break

        if discard:
            reads_filtered += 1
            continue
        else:
            reads_kept += 1

            fq1 = str(R1) + "\n"
            fq2 = str(R2) + "\n"
            output_R1.write(fq1.encode('utf-8'))
            output_R2.write(fq2.encode('utf-8'))

    print(f"reads kept: {reads_kept}, reads filtered: {reads_filtered}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Discard read pairs that contain sequencing adapters that are missed by a trimmer")
    parser.add_argument('--input_R1', required=True, default=False,
                        help="Input fasta/fastq file for first half of pair")
    
    parser.add_argument('--input_R2', required=True, default=False,
                        help="Input fasta/fastq file for second half of pair")
 
    args = parser.parse_args()
    
    # this is designed to only work in the nextflow pipeline so we put strict requirements on the input/output names
    in_ext = ".fq.gz"
    require_extension(args.input_R1, in_ext)
    require_extension(args.input_R2, in_ext)

    out_ext = "_posttrim_filter.fq.gz"
    output_R1 = args.input_R1.replace(in_ext, out_ext)
    output_R2 = args.input_R2.replace(in_ext, out_ext)

    # Illumina adapters that are occasionally leftover in reads
    S7 = "CCGAGCCCACGAGAC"
    P7 = "ATCTCGTATGCCGTCTTCTGCTTG"
    filter_sequences = [ S7, P7 ]

    filter_reads(filter_sequences, args.input_R1, args.input_R2, output_R1, output_R2)
