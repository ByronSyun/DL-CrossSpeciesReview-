import pandas as pd
import os
import gzip
from pyfaidx import Fasta
from tqdm import tqdm
import csv

SEQ_LENGTH = 8192
UNKNOWN_TISSUE_ID = 15

def get_sequence_context(chrom, pos, fasta_handle):
    """Extract 8192bp sequence context around a given position."""
    start = pos - (SEQ_LENGTH // 2)
    end = pos + (SEQ_LENGTH // 2)
    
    try:
        chrom_key = str(chrom)
        if chrom_key not in fasta_handle.keys():
            chrom_key = f"chr{chrom}" if not chrom.startswith('chr') else chrom
            if chrom_key not in fasta_handle.keys():
                return None
        
        if start < 1 or end > len(fasta_handle[chrom_key]):
            return None
            
        return fasta_handle[chrom_key][start - 1:end].seq.upper()
    except (KeyError, IndexError):
        return None

def build_negative_sequences(ref_allele, alt_allele, context_sequence):
    """Build ref/alt sequences for a negative sample."""
    variant_offset = SEQ_LENGTH // 2
    ref_sequence = context_sequence
    alt_sequence = (
        context_sequence[:variant_offset - 1] +
        alt_allele +
        context_sequence[variant_offset - 1 + len(ref_allele):]
    )
    final_ref_seq = (ref_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
    final_alt_seq = (alt_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
    return final_ref_seq, final_alt_seq

def create_negative_pool():
    """Create high-quality negative pool for ChickenGTEx."""
    
    positive_samples_file = './processed_data/chicken_lead_sqtls.tsv'
    vcf_file = './chicken_raw_data/gallus_gallus_variation.vcf.gz'
    fasta_file = './reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa'
    output_file = './processed_data/chicken_negative_pool_HIGH_QUALITY.tsv'
    MAX_NEGATIVE_SAMPLES = 225000
    
    print("Creating high-quality negative pool for ChickenGTEx...")
    print(f"Target: {MAX_NEGATIVE_SAMPLES:,} negative samples")

    try:
        positive_df = pd.read_csv(positive_samples_file, sep='\t')
        positive_locations = set()
        
        for _, row in positive_df.iterrows():
            try:
                chrom = str(row['CHR'])
                pos = int(row['POS'])
                positive_locations.add((chrom, pos))
            except (ValueError, IndexError, KeyError):
                continue
                
        print(f"Loaded {len(positive_locations):,} positive locations to exclude")
    except FileNotFoundError:
        print(f"ERROR: Positive samples file not found: {positive_samples_file}")
        return

    print("Loading reference genome...")
    try:
        if fasta_file.endswith('.gz'):
            print(f"ERROR: Reference genome is compressed. Please uncompress: gunzip {fasta_file}")
            return
        
        fasta_handle = Fasta(fasta_file, key_function=lambda key: key.split(' ')[0], sequence_always_upper=True)
        print(f"Reference genome loaded ({len(fasta_handle.keys())} chromosomes)")
    except FileNotFoundError:
        print(f"ERROR: FASTA file not found: {fasta_file}")
        return

    print("Processing VCF for negative candidates...")
    
    negative_samples = []
    vcf_line_count = 0
    snp_count = 0
    excluded_positive_count = 0
    ref_mismatch_count = 0
    boundary_error_count = 0
    sequence_error_count = 0
    
    try:
        with gzip.open(vcf_file, 'rt') as f:
            for line in tqdm(f, desc="Processing VCF"):
                if line.startswith('#'):
                    continue
                    
                vcf_line_count += 1
                parts = line.strip().split('\t')
                
                if len(parts) < 5:
                    continue
                    
                chrom, pos_str, _, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]
                
                try:
                    pos = int(pos_str)
                except ValueError:
                    continue
                    
                alt = alt_field.split(',')[0]
                
                if not (len(ref) == 1 and len(alt) == 1):
                    continue
                snp_count += 1
                
                if (str(chrom), pos) in positive_locations:
                    excluded_positive_count += 1
                    continue
                
                context = get_sequence_context(chrom, pos, fasta_handle)
                if not context:
                    boundary_error_count += 1
                    continue
                
                variant_offset = SEQ_LENGTH // 2
                if context[variant_offset - 1] != ref.upper():
                    ref_mismatch_count += 1
                    continue
                
                try:
                    ref_seq, alt_seq = build_negative_sequences(ref.upper(), alt.upper(), context)
                    
                    negative_samples.append({
                        'variant_id': f"{chrom}:{pos}_{ref}/{alt}",
                        'ref_sequence': ref_seq,
                        'alt_sequence': alt_seq,
                        'label': 0,
                        'tissue_id': UNKNOWN_TISSUE_ID
                    })
                    
                    if len(negative_samples) % 50000 == 0 and len(negative_samples) > 0:
                        print(f"Progress: {len(negative_samples):,} negative samples collected")
                    
                    if len(negative_samples) >= MAX_NEGATIVE_SAMPLES:
                        print(f"Reached target of {MAX_NEGATIVE_SAMPLES:,} samples")
                        break
                        
                except Exception:
                    sequence_error_count += 1
                    continue

    except FileNotFoundError:
        print(f"ERROR: VCF file not found: {vcf_file}")
        return

    print(f"\nFiltering statistics:")
    print(f"VCF lines: {vcf_line_count:,}, SNPs: {snp_count:,}")
    print(f"Excluded: positive={excluded_positive_count:,}, ref_mismatch={ref_mismatch_count:,}, boundary={boundary_error_count:,}, sequence={sequence_error_count:,}")
    print(f"Final negative samples: {len(negative_samples):,}")
    
    if not negative_samples:
        print("ERROR: No negative samples found")
        return
    
    final_df = pd.DataFrame(negative_samples)
    cols_to_save = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    final_df[cols_to_save].to_csv(output_file, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
    
    print(f"Saved to: {output_file}")

if __name__ == '__main__':
    create_negative_pool()
