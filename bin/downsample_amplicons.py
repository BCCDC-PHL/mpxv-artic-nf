#!/usr/bin/env python

import argparse
import csv
import json
import re
import sys

import pysam

from collections import defaultdict
from functools import reduce


def get_idxstats(bam_path):
    idxstats = {
        "ref_name": "",
        "seq_length": 0,
        "num_mapped_segments": 0,
        "num_unmapped_segments": 0,
    }

    idxstats_out = pysam.idxstats(bam_path).strip().split('\t')
    idxstats_out[3] = idxstats_out[3].split('\n')[0]

    idxstats['ref_name'] = idxstats_out[0]
    idxstats['seq_length'] = int(idxstats_out[1])
    idxstats['num_mapped_segments'] = int(idxstats_out[2])
    idxstats['num_unmapped_segments'] = int(idxstats_out[3])

    return idxstats



def merge_primers(p1, p2):
    try:
        assert p1['direction'] == p2['direction']
    except AssertionError as e:
        err_message = "Error parsing bed file. primers " + p1['primer_id'] + " and " + p2['primer_id'] + " not in same direction, cannot be merged"
        print(err_message, file=sys.stderr)
        sys.exit(1)
    if p1['direction'] == '+':
        if p1['start'] < p2['start']:
            merged_primer = p1
        else:
            merged_primer = p2
    elif p1['direction'] == '-':
        if p1['end'] > p2['end']:
            merged_primer = p1
        else:
            merged_primer = p2

    return merged_primer


def read_bed_file(bed_file_path):
    """
          { "nCoV-2019_1_LEFT": {
                                  "chrom": "MN908947.3",
                                  "start": 30,
                                  "end": 54,
                                  "primer_id": "nCoV-2019_1_LEFT",
                                  "pool_name": "1",
                                  "direction": "+"
                                 },
            "nCoV-2019_1_RIGHT": {
                                   "chrom": "MN908947.3",
                                   "start": 1183,
                                   "end": 1205,
                                   "primer_id": "nCoV-2019_1_RIGHT",
                                   "pool_name": "1",
                                   "direction": "-"
                                 },
            ...
          }
    """
    fieldnames = [
        'chrom',
        'start',
        'end',
        'primer_id',
        'pool_name',
        'direction'
    ]
    int_fields = [
        'start',
        'end',
    ]
    primers = {}
    with open(bed_file_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect='excel-tab')
        for row in reader:
            for field in int_fields:
                row[field] = int(row[field])
            if re.search('_alt', row['primer_id']):
                primer_id = re.match('(.+)_alt', row['primer_id']).group(1)
            else:
                primer_id = row['primer_id']
            if all([not primer_id == primer_key for primer_key in primers.keys()]):
                primers[primer_id] = row
            else:
                existing_primer = primers[primer_id]
                merged_primer = merge_primers(existing_primer, row)
                merged_primer['primer_id'] = primer_id
                primers[primer_id] = merged_primer

    return primers


def midpoint(start, end):
    """
    Find the midpoint between two loci
    returns: int
    """
    mid = round(start + ((end - start) / 2))
    return mid


def get_initial_amplicon_checkpoints(parsed_bed):
    """
    input: { "nCoV-2019_1_LEFT": {
                                  "chrom": "MN908947.3",
                                  "start": 30,
                                  "end": 54,
                                  "primer_id": "nCoV-2019_1_LEFT",
                                  "pool_name": "1",
                                  "direction": "+"
                                 },
            "nCoV-2019_1_RIGHT": {
                                   "chrom": "MN908947.3",
                                   "start": 1183,
                                   "end": 1205,
                                   "primer_id": "nCoV-2019_1_RIGHT",
                                   "pool_name": "1",
                                   "direction": "-"
                                 },
            ...
          }
    returns: {
               "nCoV-2019_1": [64, 618, 1173],
               "nCoV-2019_2": [1138, 1683, 2234],
               ...
             }
    """
    amplicon_checkpoints = {}
    for primer_id, primer in parsed_bed.items():
        amplicon_id = re.match("(.+)_LEFT", primer_id)
        if amplicon_id:
            amplicon_checkpoints[amplicon_id.group(1)] = {}

    for amplicon_id in amplicon_checkpoints:
        left_primer_id = amplicon_id + "_LEFT"
        right_primer_id = amplicon_id + "_RIGHT"
        left_primer = parsed_bed[left_primer_id]
        right_primer = parsed_bed[right_primer_id]
        # consider ends of amplicon to be inner ends of primers
        amplicon_start = left_primer['end']
        amplicon_end = right_primer['start']
        amplicon_midpoint = midpoint(amplicon_start, amplicon_end)
        # add a 10-base buffer inside the ends of the primers to define our first and last checkpoints
        amplicon_checkpoints[amplicon_id] = [(amplicon_start + 10), amplicon_midpoint, (amplicon_end - 10)]
    
    return amplicon_checkpoints


def subdivide_amplicon_checkpoints(amplicon_checkpoints):
    """
    input: [64, 618, 1173]
    output: [64, 341, 618, 896, 1173]
    """
    new_amplicon_checkpoints = []
    for idx in range(len(amplicon_checkpoints)):
        try:
            new_checkpoint = midpoint(amplicon_checkpoints[idx], amplicon_checkpoints[idx + 1])
            new_amplicon_checkpoints.append(new_checkpoint)
        except IndexError as e:
            pass

        all_checkpoints = amplicon_checkpoints + new_amplicon_checkpoints
        all_checkpoints.sort()
        
    return all_checkpoints


def read_pair_generator(bam, region_string=None):
    """
    from: https://www.biostars.org/p/306041/#332022
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def main(args):
    """
    """

    idxstats = get_idxstats(args.bam)
    total_reads = idxstats['num_mapped_segments'] + idxstats['num_unmapped_segments']

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bed)
    
    amplicon_checkpoints = get_initial_amplicon_checkpoints(bed)

    for _ in range(args.amplicon_subdivisions - 1):
        for amplicon_id in amplicon_checkpoints:
            amplicon_checkpoints[amplicon_id] = subdivide_amplicon_checkpoints(amplicon_checkpoints[amplicon_id])
    
    genome_checkpoints = []
    for checkpoints in amplicon_checkpoints.values():
        genome_checkpoints.extend(checkpoints)
    genome_checkpoints.sort()

    required_depth_achieved = [False] * len(genome_checkpoints)
    
    depths = [0] * 30000

    infile = pysam.AlignmentFile(args.bam, "rb")
    bam_header = infile.header.copy().to_dict()

    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    reads_processed = 0
    reads_written = 0
    reads_discarded = 0

    for segment, mate_segment in read_pair_generator(infile):
        if segment.mapping_quality < args.mapping_quality:
            reads_discarded += 2
            reads_processed += 2
            continue
        if not segment.is_proper_pair:
            reads_discarded += 2
            reads_processed += 2
            continue
        if segment.is_unmapped or segment.is_supplementary:
            reads_discarded += 2
            reads_processed += 2
            continue
        if not segment.is_paired or segment.mate_is_unmapped:
            reads_discarded += 2
            reads_processed += 2
            continue

        checkpoints_under_segment = list(filter(lambda cp: segment.reference_start <= cp and segment.reference_end >= cp, genome_checkpoints))
        checkpoints_under_mate_segment = list(filter(lambda cp: mate_segment.reference_start <= cp and mate_segment.reference_end >= cp, genome_checkpoints))
        checkpoints_under_both_segments = checkpoints_under_segment + checkpoints_under_mate_segment
        checkpoint_depths = [depths[checkpoint] for checkpoint in checkpoints_under_both_segments]
        checkpoints_achieved_required_depth = [True if depths[checkpoint] >= args.min_depth else False for checkpoint in checkpoints_under_both_segments]
        
        if not all(checkpoints_achieved_required_depth):
            outfile.write(segment)
            reads_written += 1
            outfile.write(mate_segment)
            reads_written += 1
            for checkpoint in checkpoints_under_both_segments:
                depths[checkpoint] += 1
        else:
             reads_discarded += 2

        reads_processed += 2

        if reads_processed % 1000 == 0:
            genome_checkpoint_depths = [depths[checkpoint] for checkpoint in genome_checkpoints]
            genome_checkpoints_achieved_required_depth = [True if depth >= args.min_depth else False for depth in genome_checkpoint_depths]
            if all(genome_checkpoints_achieved_required_depth):
                break

    if total_reads > 0:
        downsampling_factor = reads_written / total_reads
    else:
        downsampling_factor = 0.0

    print(','.join(["total_input_reads","reads_processed","reads_written", "reads_discarded", "downsampling_factor"]), file=sys.stderr)
    print(','.join([str(total_reads), str(reads_processed), str(reads_written), str(reads_discarded), str(round(downsampling_factor, 4))]), file=sys.stderr)

    # close up the file handles
    infile.close()
    outfile.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Downsample alignments from an amplicon scheme.')
    parser.add_argument('bam', help='bam file containing the alignment')
    parser.add_argument('--bed', help='BED file containing the amplicon scheme')
    parser.add_argument('--min-depth', type=int, default=200, help='Subsample to n coverage')
    parser.add_argument('--mapping-quality', type=int, default=20, help='Minimum mapping quality to include read in output')
    parser.add_argument('--amplicon-subdivisions', type=int, default=3, help='How many times to divide amplicons to detect coverage')
    parser.add_argument('--verbose', action='store_true', help='Debug mode')
    args = parser.parse_args()
    main(args)
