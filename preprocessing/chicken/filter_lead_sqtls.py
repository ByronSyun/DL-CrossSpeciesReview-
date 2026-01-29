import pandas as pd
import os


def filter_lead_sqtls():
    input_file = './processed_data/chicken_positive_samples_with_tissue.tsv'
    output_file = './processed_data/chicken_lead_sqtls.tsv'
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return
    
    try:
        df = pd.read_csv(input_file, sep='\t', engine='python')
    except Exception as e:
        print(f"ERROR loading data: {str(e)}")
        return
    
    required_cols = ['CHR', 'POS', 'REF', 'ALT', 'Tissue']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        return
    
    df['variant'] = df['CHR'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
    df['tissue'] = df['Tissue']
    
    gene_col = None
    for col in ['gene', 'gene_id', 'ensembl_id', 'gene_name']:
        if col in df.columns:
            gene_col = col
            break
    
    pvalue_col = None
    for col in ['pvalue', 'p_value', 'pval', 'P', 'nominal_pval']:
        if col in df.columns:
            pvalue_col = col
            break
    
    if gene_col and pvalue_col:
        df_sorted = df.sort_values(pvalue_col, ascending=True)
        lead_sqtls = df_sorted.groupby([gene_col, 'tissue']).first().reset_index()
    elif pvalue_col:
        df_sorted = df.sort_values(pvalue_col, ascending=True)
        lead_sqtls = df_sorted.groupby(['variant', 'tissue']).first().reset_index()
    else:
        if gene_col:
            lead_sqtls = df.drop_duplicates(subset=[gene_col, 'tissue', 'variant'])
        else:
            lead_sqtls = df.drop_duplicates(subset=['variant', 'tissue'])
    
    lead_sqtls = lead_sqtls.dropna(subset=['variant', 'tissue'])
    
    chr_col = 'CHR' if 'CHR' in lead_sqtls.columns else 'chr'
    if chr_col in lead_sqtls.columns:
        valid_chrs = [str(i) for i in range(1, 50)] + ['X', 'Y', 'Z', 'W', 'MT']
        lead_sqtls = lead_sqtls[lead_sqtls[chr_col].astype(str).isin(valid_chrs)]
    
    if pvalue_col and pvalue_col in lead_sqtls.columns:
        lead_sqtls = lead_sqtls[
            (lead_sqtls[pvalue_col] >= 0) & 
            (lead_sqtls[pvalue_col] <= 1) & 
            (lead_sqtls[pvalue_col].notna())
        ]
    
    lead_sqtls.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered to {len(lead_sqtls):,} lead sQTLs")
    print(f"Saved to: {output_file}")


if __name__ == "__main__":
    filter_lead_sqtls()
