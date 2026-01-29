import pandas as pd
import argparse
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')


def create_balanced_dataset(positive_file, negative_file, output_file):
    logging.info(f"Loading positive samples from: {positive_file}")
    try:
        pos_df = pd.read_csv(positive_file, sep='\t', header=None)
        num_positives = len(pos_df)
        logging.info(f"Found {num_positives:,} positive samples")

        neg_df = pd.read_csv(negative_file, sep='\t', header=None)
        logging.info(f"Loaded {len(neg_df):,} negative samples")

    except FileNotFoundError as e:
        logging.error(f"Input file not found: {e}")
        return

    if len(neg_df) < num_positives:
        logging.error(f"Negative pool ({len(neg_df)}) smaller than positive set ({num_positives})")
        return
        
    logging.info(f"Sampling {num_positives:,} negative samples")
    neg_df_sampled = neg_df.sample(n=num_positives, random_state=42)

    combined_df = pd.concat([pos_df, neg_df_sampled], ignore_index=True)
    final_df = combined_df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    final_df.columns = [0, 1, 2, 3, 4]
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    logging.info(f"Balanced benchmark created: {len(final_df):,} samples")
    logging.info(f"Positive: {len(final_df[final_df[3] == 1]):,}, Negative: {len(final_df[final_df[3] == 0]):,}")
    logging.info(f"Saved to: {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create balanced RatGTEx benchmark.")
    parser.add_argument('--positive_file', type=str, required=True, help='Path to positive samples TSV')
    parser.add_argument('--negative_file', type=str, required=True, help='Path to negative samples TSV')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save balanced benchmark TSV')
    args = parser.parse_args()
    
    create_balanced_dataset(args.positive_file, args.negative_file, args.output_file)
