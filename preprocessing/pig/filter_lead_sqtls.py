import os
import pandas as pd
import gzip
import logging
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(message)s')


def filter_lead_sqtls():
    base_dir = './pig_raw_data/PigGTEx_v0.significant_sQTL/'
    output_dir = './processed_data/'
    output_file = os.path.join(output_dir, 'pig_positive_samples_INDEPENDENT_with_tissue.tsv')

    if not os.path.isdir(base_dir):
        logging.error(f"Source directory not found: {base_dir}")
        return

    all_sqtl_files = [f for f in os.listdir(base_dir) if f.endswith('.txt.gz')]

    if not all_sqtl_files:
        logging.error(f"No sQTL files (.txt.gz) found in {base_dir}")
        return

    all_lead_sqtls_dfs = []

    for filename in tqdm(all_sqtl_files, desc="Finding lead sQTLs"):
        tissue_name = filename.split('.')[0]
        file_path = os.path.join(base_dir, filename)

        try:
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', engine='python')
                
                lead_indices = df.groupby('phenotype_id')['pval_nominal'].idxmin()
                lead_sqtls_df = df.loc[lead_indices]
                
                variant_parts = lead_sqtls_df['variant_id'].str.split('_', n=4, expand=True)
                
                temp_df = pd.DataFrame({
                    'CHR': variant_parts[0],
                    'POS': variant_parts[1],
                    'REF': variant_parts[2],
                    'ALT': variant_parts[3],
                    'Tissue': tissue_name
                })
                
                all_lead_sqtls_dfs.append(temp_df)

        except Exception as e:
            logging.error(f"Could not process file {filename}: {e}")

    if not all_lead_sqtls_dfs:
        logging.error("No lead sQTLs were processed")
        return

    combined_df = pd.concat(all_lead_sqtls_dfs, ignore_index=True)
    
    combined_df.dropna(subset=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
    combined_df['POS'] = pd.to_numeric(combined_df['POS'], errors='coerce').astype('Int64')
    combined_df.dropna(subset=['POS'], inplace=True)

    agg_df = combined_df.groupby(['CHR', 'POS', 'REF', 'ALT'])['Tissue'].apply(
        lambda x: ','.join(sorted(list(set(x))))
    ).reset_index()

    logging.info(f"Combined into {len(agg_df):,} unique independent positive sQTLs")

    agg_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Saved to: {output_file}")


if __name__ == '__main__':
    filter_lead_sqtls()
