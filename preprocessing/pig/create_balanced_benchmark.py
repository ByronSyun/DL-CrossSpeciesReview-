import pandas as pd
import numpy as np
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(message)s')


def create_balanced_benchmark(positive_file, negative_pool_file, output_file, random_seed=42):
    positive_df = pd.read_csv(
        positive_file, 
        sep='\t', 
        header=None,
        names=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    )
    n_positive = len(positive_df)
    
    negative_pool_df = pd.read_csv(
        negative_pool_file,
        sep='\t',
        header=None,
        names=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    )
    n_negative_pool = len(negative_pool_df)
    
    if n_negative_pool < n_positive:
        logging.error(f"Insufficient negatives: need {n_positive:,}, have {n_negative_pool:,}")
        return
    
    np.random.seed(random_seed)
    negative_df = negative_pool_df.sample(n=n_positive, random_state=random_seed)
    
    combined_df = pd.concat([positive_df, negative_df], ignore_index=True)
    combined_df = combined_df.sample(frac=1.0, random_state=random_seed).reset_index(drop=True)
    
    n_total = len(combined_df)
    n_pos_final = (combined_df['label'] == 1).sum()
    n_neg_final = (combined_df['label'] == 0).sum()
    
    logging.info(f"Total samples: {n_total:,}")
    logging.info(f"Positive: {n_pos_final:,}, Negative: {n_neg_final:,}")
    
    combined_df.to_csv(output_file, sep='\t', index=False, header=False)
    logging.info(f"Saved to: {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create balanced PigGTEx benchmark')
    parser.add_argument(
        '--positive', 
        type=str,
        default='./processed_data/pig_positive_samples_HIGH_QUALITY.tsv',
        help='Path to positive samples file'
    )
    parser.add_argument(
        '--negative_pool',
        type=str, 
        default='./processed_data/pig_negative_pool_HIGH_QUALITY.tsv',
        help='Path to negative pool file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='./processed_data/piggtex_silver_benchmark_balanced.tsv',
        help='Path to output file'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed'
    )
    
    args = parser.parse_args()
    
    create_balanced_benchmark(
        positive_file=args.positive,
        negative_pool_file=args.negative_pool,
        output_file=args.output,
        random_seed=args.seed
    )
