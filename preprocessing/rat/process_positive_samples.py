import pandas as pd
import logging
import argparse
import os
from pyfaidx import Fasta
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(message)s')

SEQ_LENGTH = 8192
UNKNOWN_TISSUE_ID = 15
TISSUE_NAMES = [
    'Adipose Tissue', 'Blood', 'Blood Vessel', 'Brain', 'Colon', 'Heart', 'Kidney',
    'Liver', 'Lung', 'Muscle', 'Nerve', 'Small Intestine', 'Skin', 'Spleen', 'Stomach'
]
TISSUE_MAP = {name.split(' ')[0].lower(): i for i, name in enumerate(TISSUE_NAMES)}


def get_and_force_sequences(variant_row, fasta_handle):
    chrom = str(variant_row['chrom'])
    pos = variant_row['pos']
    ref_allele_sqtl = variant_row['ref'].upper()
    alt_allele_sqtl = variant_row['alt'].upper()

    start = pos - (SEQ_LENGTH // 2)
    end = pos + (SEQ_LENGTH // 2)
    variant_offset = SEQ_LENGTH // 2

    try:
        context_sequence = fasta_handle[chrom][start - 1:end].seq.upper()
        
        ref_allele_fasta = context_sequence[variant_offset - 1 : variant_offset - 1 + len(ref_allele_sqtl)]
        is_mismatch = (ref_allele_fasta != ref_allele_sqtl)

        ref_len_diff = len(ref_allele_fasta) - len(ref_allele_sqtl)
        
        ref_seq_forced = (
            context_sequence[:variant_offset - 1] +
            ref_allele_sqtl +
            context_sequence[variant_offset - 1 + len(ref_allele_fasta):]
        )
        alt_seq_forced = (
            context_sequence[:variant_offset - 1] +
            alt_allele_sqtl +
            context_sequence[variant_offset - 1 + len(ref_allele_fasta):]
        )

        final_ref_seq = (ref_seq_forced + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
        final_alt_seq = (alt_seq_forced + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]

        return final_ref_seq, final_alt_seq, "OK", is_mismatch

    except (KeyError, IndexError):
        return None, None, "EXTRACTION_ERROR", False


def process_positive_samples(args):
    logging.info(f"Loading sQTLs from: {args.sqtl_file}")
    try:
        sqtl_df = pd.read_csv(args.sqtl_file, sep='\t')
        sqtl_df.rename(columns={
            sqtl_df.columns[4]: 'chrom',
            sqtl_df.columns[5]: 'pos',
            sqtl_df.columns[6]: 'ref',
            sqtl_df.columns[7]: 'alt',
            sqtl_df.columns[0]: 'tissue'
        }, inplace=True)
        sqtl_df['chrom'] = sqtl_df['chrom'].astype(str)
        logging.info(f"Loaded {len(sqtl_df):,} sQTL records")
    except FileNotFoundError:
        logging.error(f"sQTL file not found: {args.sqtl_file}")
        return

    try:
        fasta_handle = Fasta(args.fasta_file, key_function=lambda key: key.split(' ')[0])
    except FileNotFoundError:
        logging.error(f"FASTA file not found: {args.fasta_file}")
        return

    processed_samples = []
    mismatch_count = 0
    extraction_error_count = 0

    for _, row in tqdm(sqtl_df.iterrows(), total=len(sqtl_df), desc="Processing samples"):
        ref_seq, alt_seq, status, is_mismatch = get_and_force_sequences(row, fasta_handle)

        if is_mismatch:
            mismatch_count += 1
            
        if status == "OK":
            processed_row = {
                'variant_id': f"{row['chrom']}:{row['pos']}_{row['ref']}/{row['alt']}",
                'ref_sequence': ref_seq,
                'alt_sequence': alt_seq,
                'label': 1,
                'tissue_id': TISSUE_MAP.get(str(row['tissue']).lower(), UNKNOWN_TISSUE_ID),
                'original_tissue': row['tissue']
            }
            processed_samples.append(processed_row)
        elif status == "EXTRACTION_ERROR":
            extraction_error_count += 1

    if mismatch_count > 0:
        logging.info(f"Found {mismatch_count:,} mismatches between sQTL ref and FASTA")
    if extraction_error_count > 0:
        logging.info(f"Skipped {extraction_error_count:,} variants due to extraction errors")

    final_df = pd.DataFrame(processed_samples)
    
    if final_df.empty:
        logging.error("No samples processed successfully")
        return
        
    logging.info(f"Processed {len(final_df):,} positive samples")

    cols_to_save = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    final_df[cols_to_save].to_csv(args.output_file, sep='\t', index=False, header=False)
    
    logging.info(f"Saved to: {args.output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process positive sQTL samples, forcing sequence construction and diagnosing mismatches.")
    parser.add_argument('--sqtl_file', type=str, required=True, help='Path to the positive sQTL data file.')
    parser.add_argument('--fasta_file', type=str, required=True, help='Path to the reference genome FASTA file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the final processed TSV file.')
    args = parser.parse_args()
    process_positive_samples(args) 