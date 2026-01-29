import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

SEQ_LENGTH = 8192
UNKNOWN_TISSUE_ID = 15

TISSUE_NAMES = [
    'Adipose Tissue', 'Blood', 'Blood Vessel', 'Brain', 'Colon', 'Heart', 'Kidney',
    'Liver', 'Lung', 'Muscle', 'Nerve', 'Small Intestine', 'Skin', 'Spleen', 'Stomach'
]
TISSUE_MAP = {name.split(' ')[0].lower(): i for i, name in enumerate(TISSUE_NAMES)}


def map_tissue_to_id(tissue_string):
    if pd.isna(tissue_string):
        return UNKNOWN_TISSUE_ID
    
    tissues = str(tissue_string).split(',')
    for tissue in tissues:
        tissue_key = tissue.lower().strip().split('_')[0]
        if tissue_key in TISSUE_MAP:
            return TISSUE_MAP[tissue_key]
    
    return UNKNOWN_TISSUE_ID


def extract_perfect_match_sequences(variant_row, fasta_handle, chrom_mapping):
    simple_chrom = str(variant_row['CHR'])
    if simple_chrom not in chrom_mapping:
        return None, None, "NO_MAPPING"
    
    ncbi_chrom = chrom_mapping[simple_chrom]
    pos = int(variant_row['POS'])
    ref_allele_sqtl = str(variant_row['REF']).upper()
    alt_allele_sqtl = str(variant_row['ALT']).upper()

    start = pos - (SEQ_LENGTH // 2)
    end = pos + (SEQ_LENGTH // 2)
    variant_offset = SEQ_LENGTH // 2

    try:
        if start < 1 or end > len(fasta_handle[ncbi_chrom]):
            return None, None, "OUT_OF_BOUNDS"
            
        context_sequence = fasta_handle[ncbi_chrom][start - 1:end].seq.upper()
        
        ref_allele_fasta = context_sequence[variant_offset - 1 : variant_offset - 1 + len(ref_allele_sqtl)]
        if ref_allele_fasta != ref_allele_sqtl:
            return None, None, "MISMATCH_REJECTED"

        ref_seq = (
            context_sequence[:variant_offset - 1] +
            ref_allele_sqtl +
            context_sequence[variant_offset - 1 + len(ref_allele_fasta):]
        )
        alt_seq = (
            context_sequence[:variant_offset - 1] +
            alt_allele_sqtl +
            context_sequence[variant_offset - 1 + len(ref_allele_fasta):]
        )

        final_ref_seq = (ref_seq + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
        final_alt_seq = (alt_seq + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]

        return final_ref_seq, final_alt_seq, "PERFECT_MATCH"

    except (KeyError, IndexError):
        return None, None, "EXTRACTION_ERROR"


def create_high_quality_positive_dataset():
    positive_file = './processed_data/pig_positive_samples_INDEPENDENT_with_tissue.tsv'
    mapping_file = './processed_data/chromosome_mapping.tsv'
    fasta_file = './reference_genome/Sscrofa11.1_genomic.fna'
    output_file = './processed_data/pig_positive_samples_HIGH_QUALITY.tsv'

    mapping_df = pd.read_csv(mapping_file, sep='\t')
    chrom_mapping = dict(zip(mapping_df['simple_chrom'].astype(str), mapping_df['ncbi_chrom']))

    sqtl_df = pd.read_csv(positive_file, sep='\t')

    fasta_handle = Fasta(fasta_file, key_function=lambda key: key.split(' ')[0], sequence_always_upper=True)
    
    processed_samples = []
    perfect_matches = 0
    mismatches_rejected = 0
    out_of_bounds = 0
    other_errors = 0

    for _, row in tqdm(sqtl_df.iterrows(), total=len(sqtl_df), desc="Quality filtering"):
        ref_seq, alt_seq, status = extract_perfect_match_sequences(row, fasta_handle, chrom_mapping)

        if status == "PERFECT_MATCH":
            perfect_matches += 1
            tissue_id = map_tissue_to_id(row['Tissue'])
            
            processed_row = {
                'variant_id': f"{row['CHR']}:{row['POS']}_{row['REF']}/{row['ALT']}",
                'ref_sequence': ref_seq,
                'alt_sequence': alt_seq,
                'label': 1,
                'tissue_id': tissue_id,
                'original_tissue': row['Tissue']
            }
            processed_samples.append(processed_row)
        elif status == "MISMATCH_REJECTED":
            mismatches_rejected += 1
        elif status == "OUT_OF_BOUNDS":
            out_of_bounds += 1
        else:
            other_errors += 1

    logging.info(f"Total input: {len(sqtl_df):,}")
    logging.info(f"Perfect matches: {perfect_matches:,} ({perfect_matches/len(sqtl_df):.1%})")
    logging.info(f"Mismatches rejected: {mismatches_rejected:,}")
    logging.info(f"Out of bounds: {out_of_bounds:,}")

    if not processed_samples:
        logging.error("No high-quality samples found")
        return

    final_df = pd.DataFrame(processed_samples)

    cols_to_save = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    final_df[cols_to_save].to_csv(output_file, sep='\t', index=False, header=False)
    
    logging.info(f"Final dataset: {len(final_df):,} samples")
    logging.info(f"Saved to: {output_file}")


if __name__ == '__main__':
    create_high_quality_positive_dataset()
