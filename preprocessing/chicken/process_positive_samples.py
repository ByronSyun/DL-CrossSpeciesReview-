import os
import pandas as pd
import gzip
from tqdm import tqdm
import subprocess


def load_rsid_mapping(vcf_file_path):
    rsid_map = {}
    
    if not os.path.exists(vcf_file_path):
        return rsid_map
    
    try:
        cmd = ['zcat' if vcf_file_path.endswith('.gz') else 'cat', vcf_file_path]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        for line in process.stdout:
            line = line.strip()
            if line.startswith('#'):
                continue
                
            parts = line.split('\t')
            if len(parts) >= 5:
                chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                
                if rsid.startswith('rs'):
                    rsid_map[rsid] = (chrom, pos, ref, alt)
        
        process.wait()
        print(f"Loaded {len(rsid_map):,} rsID mappings")
        
    except Exception as e:
        print(f"ERROR loading VCF: {e}")
        
    return rsid_map


def process_positive_samples():
    possible_dirs = [
        './chicken_raw_data/sQTL_summary_stats/sQTLs/',
        './chicken_raw_data/sQTL_summary_stats/',
        './chicken_raw_data/sQTL_data/',
        './chicken_raw_data/ChickenGTEx_pilot_phase/molQTL mapping/',
        './chicken_raw_data/molQTL_data/',
        './chicken_raw_data/',
    ]
    
    base_dir = None
    for dir_path in possible_dirs:
        if os.path.isdir(dir_path):
            files = os.listdir(dir_path)
            sqtl_keywords = ['sqtl', 'splice', 'molqtl', 'signifpairs', 'cis_independent_qtl', 'qtl']
            sqtl_files = [f for f in files if any(keyword in f.lower() for keyword in sqtl_keywords) 
                         and f.endswith(('.txt', '.tsv', '.gz')) 
                         and not f.endswith(('.tar.gz', '.zip', '.tar'))]
            
            if sqtl_files:
                base_dir = dir_path
                break
    
    output_dir = './processed_data/'
    output_file = os.path.join(output_dir, 'chicken_positive_samples_with_tissue.tsv')
    os.makedirs(output_dir, exist_ok=True)
    
    vcf_file_path = './chicken_raw_data/gallus_gallus_variation.vcf.gz'
    rsid_map = load_rsid_mapping(vcf_file_path)
    
    if base_dir is None:
        print("ERROR: No sQTL data directory found")
        return

    all_files = os.listdir(base_dir)
    sqtl_keywords = ['sqtl', 'splice', 'molqtl', 'signifpairs', 'cis_independent_qtl', 'qtl']
    sqtl_files = [f for f in all_files if any(keyword in f.lower() for keyword in sqtl_keywords)
                 and f.endswith(('.txt', '.tsv', '.gz', '.txt.gz', '.tsv.gz')) 
                 and not f.endswith(('.tar.gz', '.zip', '.tar'))]
    
    if not sqtl_files:
        print(f"ERROR: No sQTL files found in {base_dir}")
        return

    all_positive_samples_dfs = []

    for filename in tqdm(sqtl_files, desc="Processing sQTL files"):
        parts = filename.split('.')
        if len(parts) >= 3 and parts[0] == 'ChickenGTEx':
            tissue_name = parts[1]
        else:
            tissue_name = parts[0]
        
        file_path = os.path.join(base_dir, filename)

        try:
            if filename.endswith('.gz'):
                with gzip.open(file_path, 'rt') as f:
                    df = pd.read_csv(f, sep='\t', engine='python')
            else:
                df = pd.read_csv(file_path, sep='\t', engine='python')

            variant_col = None
            for col in ['variant_id', 'snp_id', 'rsid', 'SNP', 'id']:
                if col in df.columns:
                    variant_col = col
                    break
            
            if variant_col is None:
                continue

            essential_cols = {'variant': variant_col}
            
            for col in ['chr', 'chromosome', 'chrom', 'CHROM']:
                if col in df.columns:
                    essential_cols['chr'] = col
                    break
            
            for col in ['pos', 'position', 'POS', 'start']:
                if col in df.columns:
                    essential_cols['pos'] = col
                    break
            
            for col in ['gene_id', 'gene', 'ensembl_id', 'gene_name']:
                if col in df.columns:
                    essential_cols['gene'] = col
                    break
            
            for col in ['pvalue', 'p_value', 'pval', 'P', 'nominal_pval']:
                if col in df.columns:
                    essential_cols['pvalue'] = col
                    break

            if 'variant' in essential_cols and len(essential_cols) < 3:
                sample_variants = df[essential_cols['variant']].head(3).tolist()
                use_rsid_mapping = False
                use_direct_parsing = False
                
                try:
                    test_variant = sample_variants[0]
                    
                    if test_variant.startswith('rs'):
                        if len(rsid_map) > 0:
                            use_rsid_mapping = True
                        else:
                            continue
                    else:
                        variant_parts = test_variant.split('_')
                        if len(variant_parts) >= 4:
                            use_direct_parsing = True
                        else:
                            continue
                            
                except Exception:
                    continue
            else:
                use_rsid_mapping = False
                use_direct_parsing = False

            if len(essential_cols) < 3 and not use_rsid_mapping and not use_direct_parsing:
                continue

            if use_rsid_mapping:
                variant_ids = df[essential_cols['variant']]
                mapped_data = []
                
                for rsid in variant_ids:
                    if rsid in rsid_map:
                        chrom, pos, ref, alt = rsid_map[rsid]
                        mapped_data.append({
                            'CHR': chrom,
                            'POS': int(pos),
                            'REF': ref,
                            'ALT': alt,
                            'Tissue': tissue_name
                        })
                
                tissue_df = pd.DataFrame(mapped_data)
                
            elif use_direct_parsing:
                variant_parts = df[essential_cols['variant']].str.split('_', n=4, expand=True)
                
                tissue_df = pd.DataFrame({
                    'CHR': variant_parts[0],
                    'POS': variant_parts[1], 
                    'REF': variant_parts[2],
                    'ALT': variant_parts[3],
                    'Tissue': tissue_name
                })
                
                tissue_df.dropna(subset=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
                tissue_df['POS'] = pd.to_numeric(tissue_df['POS'], errors='coerce')
                tissue_df.dropna(subset=['POS'], inplace=True)
                tissue_df['POS'] = tissue_df['POS'].astype(int)
                
            else:
                selected_cols = []
                col_mapping = {}
                
                for std_name, col_name in essential_cols.items():
                    if col_name in df.columns:
                        selected_cols.append(col_name)
                        col_mapping[col_name] = std_name

                tissue_df = df[selected_cols].copy()
                tissue_df = tissue_df.rename(columns=col_mapping)
                tissue_df['tissue'] = tissue_name
                
                for ref_col in ['ref', 'REF', 'reference_allele']:
                    if ref_col in df.columns:
                        tissue_df['ref'] = df[ref_col]
                        break
                
                for alt_col in ['alt', 'ALT', 'alternative_allele']:
                    if alt_col in df.columns:
                        tissue_df['alt'] = df[alt_col]
                        break
            
            all_positive_samples_dfs.append(tissue_df)

        except Exception as e:
            print(f"ERROR processing {filename}: {str(e)}")
            continue

    if not all_positive_samples_dfs:
        print("ERROR: No sQTL data processed")
        return

    combined_df = pd.concat(all_positive_samples_dfs, ignore_index=True)
    
    if all(col in combined_df.columns for col in ['CHR', 'POS', 'REF', 'ALT', 'Tissue']):
        combined_df.dropna(subset=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
        
        agg_df = combined_df.groupby(['CHR', 'POS', 'REF', 'ALT'])['Tissue'].apply(
            lambda x: ','.join(sorted(list(set(x))))
        ).reset_index()
        
        print(f"Aggregated to {len(agg_df):,} unique sQTLs")
        agg_df.to_csv(output_file, sep='\t', index=False)
        
    else:
        combined_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Saved to: {output_file}")


if __name__ == "__main__":
    process_positive_samples()
