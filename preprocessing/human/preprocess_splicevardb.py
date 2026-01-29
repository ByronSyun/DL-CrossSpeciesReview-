import os
import pandas as pd
from pyfaidx import Fasta
from pathlib import Path
import argparse

WINDOW_SIZE = 8192

def parse_hg38_variant(hg38_string):
    """Parse hg38 variant string in format 'chrom-pos-ref-alt'"""
    try:
        parts = hg38_string.strip('"').split('-')
        if len(parts) != 4:
            return None, None, None, None
        
        chrom, pos_str, ref, alt = parts
        
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
        
        pos_1based = int(pos_str)  # 1-based
        
        return chrom, pos_1based, ref, alt
    except (ValueError, IndexError):
        return None, None, None, None

def parse_hg38_info(hg38_info_path, positives_only=False, drop_ambiguous_location=False, dedup_by_hg38=False):
    """Load TSV, filter to Normal/Splice-altering, apply optional filters, keep valid SNVs."""
    try:
        df = pd.read_csv(
            hg38_info_path,
            sep='\t',
            header=0,
            skiprows=[1],
            on_bad_lines='warn',
            dtype=str
        )
        if df.shape[1] < 9:
            df = pd.read_csv(hg38_info_path, sep='\s+', header=0, skiprows=[1], on_bad_lines='warn', dtype=str)
    except Exception:
        df = pd.read_csv(hg38_info_path, sep='\t', header=0, skiprows=[1], on_bad_lines='warn', engine='python', dtype=str)

    df.columns = df.columns.str.strip()

    valid_classifications = ['Normal', 'Splice-altering']
    df['classification'] = df['classification'].str.strip().str.replace('"', '')
    df_filtered = df[df['classification'].isin(valid_classifications)]

    if positives_only:
        df_filtered = df_filtered[df_filtered['classification'] == 'Splice-altering']

    if drop_ambiguous_location and 'location' in df_filtered.columns:
        df_filtered = df_filtered[df_filtered['location'].str.strip().ne('Exonic,Intronic')]

    if dedup_by_hg38 and 'hg38' in df_filtered.columns:
        df_filtered = df_filtered.drop_duplicates(subset=['hg38'], keep='first')

    valid_nucleotides = {'A', 'T', 'G', 'C'}
    valid_variants = []

    for _, row in df_filtered.iterrows():
        hg38_string = row.get('hg38')
        if not isinstance(hg38_string, str):
            continue
        chrom, pos_1based, ref, alt = parse_hg38_variant(hg38_string)
        if (chrom and pos_1based and ref in valid_nucleotides and alt in valid_nucleotides and len(ref) == 1 and len(alt) == 1):
            variant_data = row.copy()
            variant_data['chrom'] = chrom
            variant_data['pos'] = pos_1based
            variant_data['ref'] = ref
            variant_data['alt'] = alt
            variant_data['location_clean'] = row.get('location', '').strip('"')
            valid_variants.append(variant_data)

    df_snv = pd.DataFrame(valid_variants)
    return df_snv

def extract_sequence_with_validation(fasta, chrom, pos, ref, window_size=WINDOW_SIZE):
    """Extract sequence and validate reference nucleotide."""
    def _resolve_chrom_name(fasta_obj, chrom_name):
        try:
            if chrom_name in fasta_obj.keys():
                return chrom_name
            if chrom_name.startswith('chr'):
                alt = chrom_name.replace('chr', '', 1)
                if alt in fasta_obj.keys():
                    return alt
            else:
                alt = f"chr{chrom_name}"
                if alt in fasta_obj.keys():
                    return alt
        except Exception:
            pass
        return None

    actual_chrom = _resolve_chrom_name(fasta, chrom)
    if actual_chrom is None:
        return None, f"Chromosome {chrom} not found in FASTA."

    pos_0based = pos - 1
    half_window = window_size // 2
    start = max(0, pos_0based - half_window)
    end = pos_0based + half_window

    try:
        sequence = str(fasta[actual_chrom][start:end]).upper()
        ref_pos_in_seq = pos_0based - start
        if ref_pos_in_seq < len(sequence):
            actual_ref = sequence[ref_pos_in_seq]
            if actual_ref != ref:
                return None, f"Reference mismatch: expected {ref}, got {actual_ref}"
        return sequence, None
    except Exception as e:
        return None, str(e)

def build_benchmark_dataset(df, fasta_path, output_path):
    """Build dataset with extracted sequences; write TSV."""
    if fasta_path and Path(fasta_path).exists():
        fasta = Fasta(fasta_path)
        results = []
        failed_extractions = 0
        for _, row in df.iterrows():
            variant_id = f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
            sequence, error = extract_sequence_with_validation(fasta, row['chrom'], row['pos'], row['ref'])
            if sequence is None:
                failed_extractions += 1
                continue
            label = 'splice-altering' if row['classification'] == 'Splice-altering' else 'normal'
            results.append({
                'variant_id': variant_id,
                'chrom': row['chrom'],
                'pos_0based': row['pos'] - 1,
                'ref': row['ref'],
                'alt': row['alt'],
                'label': label,
                'location': row['location_clean'],
                'sequence': sequence
            })
        output_df = pd.DataFrame(results)
    else:
        results = []
        for _, row in df.iterrows():
            variant_id = f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
            label = 'splice-altering' if row['classification'] == 'Splice-altering' else 'normal'
            results.append({
                'variant_id': variant_id,
                'chrom': row['chrom'],
                'pos_0based': row['pos'] - 1,
                'ref': row['ref'],
                'alt': row['alt'],
                'label': label,
                'location': row['location_clean']
            })
        output_df = pd.DataFrame(results)

    output_df.to_csv(output_path, sep='\t', index=False)
    return output_df

def verify_reference_genome(fasta_path):
    """Minimal FASTA existence check."""
    if not fasta_path or not Path(fasta_path).exists():
        return False
    try:
        _ = Fasta(fasta_path)
        return True
    except Exception:
        return False

def main():
    parser = argparse.ArgumentParser(description='Prepare SpliceVarDB dataset for Evo2 evaluation')
    parser.add_argument('--hg38_info', required=True, help='Path to the raw SpliceVarDB download TSV file.')
    parser.add_argument('--fasta', required=True, help='Path to hg38 reference genome FASTA file.')
    parser.add_argument('--output', required=True, help='Output path for the processed dataset TSV file.')
    parser.add_argument('--positives_only', action='store_true', help='Keep only Splice-altering variants')
    parser.add_argument('--drop_ambiguous_location', action='store_true', help='Drop rows with location == "Exonic,Intronic"')
    parser.add_argument('--dedup_by_hg38', action='store_true', help='Drop duplicate variants based on hg38 coordinate')
    parser.add_argument('--variant_list_out', default='', help='Optional path to write hg38/classification/location variant list')
    args = parser.parse_args()

    ok = verify_reference_genome(args.fasta)
    if not ok:
        raise FileNotFoundError(f"FASTA not available or unreadable: {args.fasta}")

    df_variants = parse_hg38_info(
        args.hg38_info,
        positives_only=args.positives_only,
        drop_ambiguous_location=args.drop_ambiguous_location,
        dedup_by_hg38=args.dedup_by_hg38,
    )

    # Write variant list (hg38, classification, location) for AlphaGenome/benchmarks
    try:
        if not args.variant_list_out:
            variant_dir = os.path.dirname(args.output) or '.'
            variant_path = os.path.join(variant_dir, 'splicevardb_filter_benchmark.tsv')
        else:
            variant_path = args.variant_list_out
        if len(df_variants) > 0:
            tmp = df_variants.copy()
            # Normalized hg38 in chr-pos-ref-alt format from parsed fields
            tmp['hg38_norm'] = tmp.apply(lambda r: f"{r['chrom']}-{int(r['pos'])}-{r['ref']}-{r['alt']}", axis=1)
            out_df = pd.DataFrame({
                'hg38': tmp['hg38_norm'],
                'classification': tmp['classification'],
                'location': tmp['location_clean'].astype(str).str.replace('"', '')
            })
            out_df.to_csv(variant_path, sep='\t', index=False)
            print(f"Saved {len(out_df)} rows -> {variant_path}")
    except Exception:
        pass

    final_dataset = build_benchmark_dataset(df_variants, args.fasta, args.output)
    print(f"Saved {len(final_dataset)} rows -> {args.output}")


if __name__ == "__main__":
    main() 