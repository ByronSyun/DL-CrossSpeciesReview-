import os
import pandas as pd
import gzip
import logging
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(message)s')


def process_positive_samples():
    base_dir = './pig_raw_data/PigGTEx_v0.significant_sQTL/'
    output_dir = './processed_data/'
    output_file = os.path.join(output_dir, 'pig_positive_samples_with_tissue.tsv')

    os.makedirs(output_dir, exist_ok=True)

    if not os.path.isdir(base_dir):
        logging.error(f"Source directory not found: {base_dir}")
        return

    all_sqtl_files = [f for f in os.listdir(base_dir) if f.endswith('.txt.gz')]

    if not all_sqtl_files:
        logging.error(f"No sQTL files (.txt.gz) found in {base_dir}")
        return

    all_positive_samples_dfs = []

    for filename in tqdm(all_sqtl_files, desc="Processing tissues"):
        tissue_name = filename.split('.')[0]
        file_path = os.path.join(base_dir, filename)

        try:
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', engine='python')

                if 'variant_id' not in df.columns:
                    logging.warning(f"'variant_id' column not found in {filename}, skipping")
                    continue

                variant_parts = df['variant_id'].str.split('_', n=4, expand=True)

                temp_df = pd.DataFrame({
                    'CHR': variant_parts[0],
                    'POS': variant_parts[1],
                    'REF': variant_parts[2],
                    'ALT': variant_parts[3],
                    'Tissue': tissue_name
                })
                
                all_positive_samples_dfs.append(temp_df)

        except Exception as e:
            logging.error(f"Could not process file {filename}: {e}")

    if not all_positive_samples_dfs:
        logging.error("No data was successfully processed")
        return

    combined_df = pd.concat(all_positive_samples_dfs, ignore_index=True)
    
    combined_df.dropna(subset=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
    combined_df['POS'] = pd.to_numeric(combined_df['POS'], errors='coerce')
    combined_df.dropna(subset=['POS'], inplace=True)
    combined_df['POS'] = combined_df['POS'].astype(int)

    agg_df = combined_df.groupby(['CHR', 'POS', 'REF', 'ALT'])['Tissue'].apply(
        lambda x: ','.join(sorted(list(set(x))))
    ).reset_index()

    logging.info(f"Combined into {len(agg_df):,} unique positive sQTLs from {len(all_sqtl_files)} tissues")

    agg_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Saved to: {output_file}")


if __name__ == '__main__':
    process_positive_samples()
