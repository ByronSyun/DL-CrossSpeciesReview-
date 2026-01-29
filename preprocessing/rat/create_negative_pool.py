import pandas as pd
import logging
import argparse
import os
from pyfaidx import Fasta
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(message)s')

SEQ_LENGTH = 8192
UNKNOWN_TISSUE_ID = 15


def get_sequence_context(chrom, pos, fasta_handle):
    start = pos - (SEQ_LENGTH // 2)
    end = pos + (SEQ_LENGTH // 2)
    try:
        return fasta_handle[str(chrom)][start - 1:end].seq.upper()
    except (KeyError, IndexError):
        return None


def build_negative_seqs(row, context_sequence):
    ref_allele_vcf = row['ref'].upper()
    alt_allele_vcf = row['alt'].upper()
    variant_offset = SEQ_LENGTH // 2
    
    ref_sequence = context_sequence
    alt_sequence = (
        context_sequence[:variant_offset - 1] +
        alt_allele_vcf +
        context_sequence[variant_offset - 1 + len(ref_allele_vcf):]
    )
    
    final_ref_seq = (ref_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
    final_alt_seq = (alt_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]

    return final_ref_seq, final_alt_seq


def create_negative_pool(args):
    logging.info("Loading positive sample locations")
    
    try:
        positive_df = pd.read_csv(args.positive_samples_file, sep='\t', header=None, usecols=[0])
        positive_locations = set()
        for variant_id in positive_df[0]:
            try:
                chrom, rest = variant_id.split(':', 1)
                pos = int(rest.split('_', 1)[0])
                positive_locations.add((str(chrom), pos))
            except (ValueError, IndexError):
                continue
        logging.info(f"Loaded {len(positive_locations)} positive locations to exclude")
    except FileNotFoundError:
        logging.error(f"Positive samples file not found: {args.positive_samples_file}")
        return

    logging.info("Processing VCF for negative candidates")
    fasta_handle = Fasta(args.fasta_file, key_function=lambda key: key.split(' ')[0])
    
    negative_samples = []
    vcf_line_count = 0
    with open(args.vcf_file, 'r') as f:
        for line in tqdm(f, desc="Parsing VCF"):
            if not line.startswith('#'):
                vcf_line_count += 1
                parts = line.strip().split('\t')
                chrom, pos, _, ref, alt_field = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
                alt = alt_field.split(',')[0]
                
                if not (len(ref) == 1 and len(alt) == 1):
                    continue
                
                if (str(chrom), pos) in positive_locations:
                    continue
                
                context = get_sequence_context(chrom, pos, fasta_handle)
                if not context:
                    continue
                
                variant_offset = SEQ_LENGTH // 2
                if context[variant_offset - 1] != ref.upper():
                    continue

                ref_seq, alt_seq = build_negative_seqs({'ref': ref, 'alt': alt}, context)
                negative_samples.append({
                    'variant_id': f"{chrom}:{pos}_{ref}/{alt}",
                    'ref_sequence': ref_seq,
                    'alt_sequence': alt_seq,
                    'label': 0,
                    'tissue_id': UNKNOWN_TISSUE_ID
                })

    logging.info(f"Scanned {vcf_line_count:,} VCF variants")
    logging.info(f"Found {len(negative_samples):,} negative samples")

    if not negative_samples:
        logging.warning("No negative samples found")
    
    final_df = pd.DataFrame(negative_samples)
    
    cols_to_save = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    final_df[cols_to_save].to_csv(args.output_file, sep='\t', index=False, header=False)
    
    logging.info(f"Saved to: {args.output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create negative samples from VCF, excluding positive locations.")
    parser.add_argument('--positive_samples_file', type=str, required=True, help='Path to processed positive samples TSV')
    parser.add_argument('--vcf_file', type=str, required=True, help='Path to Ensembl VCF file')
    parser.add_argument('--fasta_file', type=str, required=True, help='Path to reference genome FASTA')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save negative pool TSV')
    args = parser.parse_args()
    create_negative_pool(args)
