import pandas as pd
import os
import pysam
import csv


def create_high_quality_positive():
    input_file = './processed_data/chicken_lead_sqtls.tsv'
    output_file = './processed_data/chicken_positive_samples_HIGH_QUALITY.tsv'
    reference_file = './reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa'
    
    window_size = 8192
    half_window = window_size // 2
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return
    
    if not os.path.exists(reference_file):
        print(f"ERROR: Reference genome not found: {reference_file}")
        return
    
    try:
        df = pd.read_csv(input_file, sep='\t', engine='python')
    except Exception as e:
        print(f"ERROR loading data: {str(e)}")
        return
    
    if 'CHR' in df.columns and 'POS' in df.columns:
        required_cols = ['variant', 'CHR', 'POS', 'tissue']
        chr_col, pos_col, tissue_col = 'CHR', 'POS', 'tissue'
    else:
        required_cols = ['variant', 'chr', 'pos', 'tissue'] 
        chr_col, pos_col, tissue_col = 'chr', 'pos', 'tissue'
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        return
    
    try:
        ref_fasta = pysam.FastaFile(reference_file)
        available_chroms = ref_fasta.references
    except Exception as e:
        print(f"ERROR loading reference genome: {str(e)}")
        return
    
    high_quality_samples = []
    errors = {
        'chromosome_not_found': 0,
        'position_out_of_bounds': 0,
        'sequence_extraction_failed': 0,
        'reference_mismatch': 0,
        'invalid_sequence': 0,
        'missing_alleles': 0
    }
    
    for idx, row in df.iterrows():
        try:
            variant_id = row['variant']
            chrom = str(row[chr_col])
            pos = int(row[pos_col])
            tissue = row[tissue_col]
            
            if chrom.startswith('chr'):
                chrom_ref = chrom
                chrom_simple = chrom[3:]
            else:
                chrom_ref = f"chr{chrom}" if chrom not in ['MT', 'M'] else chrom
                chrom_simple = chrom
            
            if chrom_ref in available_chroms:
                chrom_final = chrom_ref
            elif chrom_simple in available_chroms:
                chrom_final = chrom_simple
            elif chrom in available_chroms:
                chrom_final = chrom
            else:
                errors['chromosome_not_found'] += 1
                continue
            
            chrom_length = ref_fasta.get_reference_length(chrom_final)
            
            start_pos = pos - half_window
            end_pos = pos + half_window
            
            min_buffer = 1000
            if start_pos < min_buffer or end_pos > (chrom_length - min_buffer):
                errors['position_out_of_bounds'] += 1
                continue
            
            try:
                sequence = ref_fasta.fetch(chrom_final, start_pos-1, start_pos-1+window_size)
                if len(sequence) < window_size * 0.75:
                    errors['sequence_extraction_failed'] += 1
                    errors['sequence_extraction_failed'] += 1
                    continue
            except Exception:
                errors['sequence_extraction_failed'] += 1
                continue
            
            if not sequence or len(sequence) < window_size * 0.75:
                errors['invalid_sequence'] += 1
                continue
            
            n_count = sequence.upper().count('N')
            if n_count > window_size * 0.1:
                errors['invalid_sequence'] += 1
                continue
            
            ref_pos_in_seq = pos - start_pos
            if ref_pos_in_seq < 0 or ref_pos_in_seq >= len(sequence):
                errors['sequence_extraction_failed'] += 1
                continue
            ref_allele_seq = sequence[ref_pos_in_seq].upper()
            
            ref_allele = None
            alt_allele = None
            
            ref_col_name = 'REF' if 'REF' in row else 'ref'
            alt_col_name = 'ALT' if 'ALT' in row else 'alt'
            
            if ref_col_name in row and pd.notna(row[ref_col_name]):
                ref_allele = str(row[ref_col_name]).upper()
            if alt_col_name in row and pd.notna(row[alt_col_name]):
                alt_allele = str(row[alt_col_name]).upper()
            
            if ref_allele and ref_allele.upper() != ref_allele_seq:
                errors['reference_mismatch'] += 1
                continue
            
            if not ref_allele or not alt_allele:
                errors['missing_alleles'] += 1
                continue
            
            alt_sequence = sequence.upper()
            if ref_allele and alt_allele and ref_pos_in_seq < len(sequence):
                actual_ref = sequence[ref_pos_in_seq].upper()
                if actual_ref != ref_allele.upper():
                    errors['reference_mismatch'] += 1
                    continue
                
                alt_sequence = list(sequence.upper())
                alt_sequence[ref_pos_in_seq] = alt_allele.upper()
                alt_sequence = ''.join(alt_sequence)
            
            hq_record = {
                'variant_id': variant_id,
                'ref_sequence': sequence.upper(),
                'alt_sequence': alt_sequence,
                'chr': chrom_simple,
                'pos': pos,
                'ref': ref_allele,
                'alt': alt_allele,
                'tissue': tissue,
                'start_pos': start_pos,
                'end_pos': end_pos
            }
            
            excluded_cols = ['chr', 'pos', 'CHR', 'POS', 'ref', 'alt', 'REF', 'ALT']
            for col in df.columns:
                if col not in hq_record and col not in excluded_cols:
                    hq_record[col] = row[col]
            
            high_quality_samples.append(hq_record)
            
        except Exception:
            errors['sequence_extraction_failed'] += 1
            continue
    
    ref_fasta.close()
    
    print(f"High-quality: {len(high_quality_samples):,}/{len(df):,} ({len(high_quality_samples)/len(df)*100:.1f}%)")
    
    if not high_quality_samples:
        print("ERROR: No high-quality samples generated")
        return
    
    hq_df = pd.DataFrame(high_quality_samples)
    
    TISSUE_NAMES = [
        'Adipose Tissue', 'Blood', 'Blood Vessel', 'Brain', 'Colon', 'Heart', 'Kidney',
        'Liver', 'Lung', 'Muscle', 'Nerve', 'Small Intestine', 'Skin', 'Spleen', 'Stomach'
    ]
    TISSUE_MAP = {
        'adipose': 0, 'blood': 1, 'vessel': 2, 'brain': 3, 'colon': 4, 'heart': 5,
        'kidney': 6, 'liver': 7, 'lung': 8, 'muscle': 9, 'nerve': 10, 'small': 11,
        'skin': 12, 'spleen': 13, 'stomach': 14
    }
    UNKNOWN_TISSUE_ID = 15
    
    def map_chicken_tissue_to_id(tissue_string):
        if pd.isna(tissue_string):
            return UNKNOWN_TISSUE_ID
        
        tissues = str(tissue_string).split(',')
        for tissue in tissues:
            tissue_clean = tissue.lower().strip()
            tissue_mappings = {
                'adipose': 'adipose', 'blood': 'blood', 'brain': 'brain', 'heart': 'heart',
                'kidney': 'kidney', 'liver': 'liver', 'lung': 'lung', 'muscle': 'muscle',
                'skin': 'skin', 'spleen': 'spleen', 'stomach': 'stomach'
            }
            
            mapped_tissue = tissue_mappings.get(tissue_clean, 'unknown')
            if mapped_tissue != 'unknown' and mapped_tissue in TISSUE_MAP:
                return TISSUE_MAP[mapped_tissue]
        
        return UNKNOWN_TISSUE_ID
    
    simplified_df = pd.DataFrame({
        'variant_id': hq_df['variant_id'],
        'ref_sequence': hq_df['ref_sequence'], 
        'alt_sequence': hq_df['alt_sequence'],
        'label': 1,
        'tissue_id': hq_df['tissue'].apply(map_chicken_tissue_to_id),
        'original_tissue': hq_df['tissue']
    })
    
    known_tissues_only = simplified_df[simplified_df['tissue_id'] != UNKNOWN_TISSUE_ID].copy()
    print(f"Known tissue samples: {len(known_tissues_only):,} (filtered {len(simplified_df) - len(known_tissues_only):,} unknown)")
    simplified_df = known_tissues_only
    
    tissue_counts = simplified_df['tissue_id'].value_counts().sort_index()
    total_samples = len(simplified_df)
    target_total = 225000
    
    if total_samples <= target_total:
        final_df = simplified_df
        print(f"Using all {total_samples:,} samples")
    else:
        sampled_dfs = []
        unique_tissues = simplified_df['tissue_id'].unique()
        min_samples_per_tissue = 1000
        reserved_samples = len(unique_tissues) * min_samples_per_tissue
        remaining_samples = target_total - reserved_samples
        
        if remaining_samples < 0:
            min_samples_per_tissue = target_total // len(unique_tissues)
            remaining_samples = target_total - len(unique_tissues) * min_samples_per_tissue
        
        for tissue_id in sorted(unique_tissues):
            tissue_samples = simplified_df[simplified_df['tissue_id'] == tissue_id]
            available = len(tissue_samples)
            proportion = available / total_samples
            additional = int(remaining_samples * proportion)
            allocated = min_samples_per_tissue + additional
            final_allocation = min(allocated, available)
            
            if available > 0:
                if final_allocation < available:
                    sampled = tissue_samples.sample(n=final_allocation, random_state=42)
                else:
                    sampled = tissue_samples
                sampled_dfs.append(sampled)
        
        final_df = pd.concat(sampled_dfs, ignore_index=True)
        final_df = final_df.sample(frac=1, random_state=42).reset_index(drop=True)
        print(f"Sampled: {len(final_df):,} samples")
    
    output_df = final_df[['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']]
    output_df.to_csv(output_file, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
    
    tissue_mapping_file = output_file.replace('.tsv', '_sptransformer_tissue_mapping.tsv')
    tissue_mapping_data = []
    for tissue_name, tissue_id in TISSUE_MAP.items():
        tissue_mapping_data.append({'tissue_name': tissue_name, 'tissue_id': tissue_id, 'full_name': TISSUE_NAMES[tissue_id]})
    tissue_mapping_data.append({'tissue_name': 'unknown', 'tissue_id': UNKNOWN_TISSUE_ID, 'full_name': 'Unknown/Other'})
    
    tissue_mapping_df = pd.DataFrame(tissue_mapping_data)
    tissue_mapping_df.to_csv(tissue_mapping_file, sep='\t', index=False)
    
    print(f"Saved to: {output_file}")
    print(f"Summary: {len(output_df):,} samples, {output_df['variant_id'].nunique():,} unique variants")


if __name__ == "__main__":
    create_high_quality_positive()
