#!/usr/bin/env python

import argparse
import csv
import json
import math
import re
import sys


def parse_gff_attribute(attribute_str):
    attribute_lst = attribute_str.split(';')
    attribute_dict = {}
    for attr in attribute_lst:
        [key, val] = attr.split('=')
        attribute_dict[key] = val
        
    return attribute_dict


def parse_gff(gff_path):
    parsed_gff = []

    gff_fieldnames = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute',
        
    ]
    with open(gff_path, 'r') as f:
        for row in f:
            if row.startswith('#'):
                pass
            else:
                break
        reader = csv.DictReader(f, fieldnames=gff_fieldnames, dialect='excel-tab')
        for row in reader:
            if row['seqname'].startswith('#'):
                pass
            else:
                row['attribute'] = parse_gff_attribute(row['attribute'])
                parsed_gff.append(row)

    return parsed_gff


def parse_variants(variants_path):
    parsed_variants = []
    with open(variants_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            parsed_variants.append(row)

    return parsed_variants



def add_ref_codon_pos(variant, cds_features_by_id):
    if variant['GFF_FEATURE'] == 'NA':
        variant['CODON_POS'] = 'NA'
    else:
        gff_feature = cds_features_by_id[variant['GFF_FEATURE']]
        variant_pos = int(variant['POS'])
        cds_start = int(gff_feature['start'])
        ref_codon_pos = int(math.floor((variant_pos - (cds_start - 3)) / 3))
        variant['CODON_POS'] = ref_codon_pos
    return variant

def add_aa_mutation_name(variant):
    feature = variant['GFF_FEATURE']
    ref_codon_pos = variant['CODON_POS']
    ref_aa = variant['REF_AA']
    alt_aa = variant['ALT_AA']
    if feature == 'NA' or ref_aa == 'NA' or ref_codon_pos == 'NA' or alt_aa == 'NA':
        variant['MUT_NAME'] = 'NA'
    else:
        variant['MUT_NAME'] = feature + ':' + ref_aa + str(ref_codon_pos) + alt_aa

    return variant

def main(args):

    gff = []
    cds_by_name = {}
    
    if args.gff:
        gff = parse_gff(args.gff)
        for gff_feature in gff:
            if gff_feature['feature'] == 'CDS':
                cds_by_name[gff_feature['attribute']['ID']] = gff_feature

    variants = parse_variants(args.variants)    

    for variant in variants:
        if args.gff:
            variant = add_ref_codon_pos(variant, cds_by_name)
            variant = add_aa_mutation_name(variant)
        else:
            variant['CODON_POS'] = 'NA'
            variant['MUT_NAME'] = 'NA'

    output_fieldnames = [
        'REGION',
        'POS',
	'REF',
	'ALT',
	'REF_DP',
	'REF_RV',
	'REF_QUAL',
        'ALT_DP',
	'ALT_RV',
	'ALT_QUAL',
	'ALT_FREQ',
	'TOTAL_DP',
	'PVAL',
	'PASS',
	'GFF_FEATURE',
	'REF_CODON',
        'REF_AA',
	'ALT_CODON',
        'ALT_AA',
        'CODON_POS',
        'MUT_NAME',
    ]

    csv.register_dialect('unix-tab', delimiter='\t', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix-tab')
    writer.writeheader()
    for variant in variants:
        writer.writerow(variant)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Add codon position column to ivar variants output')
    parser.add_argument('variants', help='bam file containing the alignment')
    parser.add_argument('-g', '--gff', help='GFF file used to label GFF_FEATURE field of ivar variants outptut')
    args = parser.parse_args()
    main(args)
