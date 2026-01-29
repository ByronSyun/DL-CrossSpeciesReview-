# Embedding-based MLP Training 

This directory contains  MLP training scripts for all species using embeddings from DNABERT-2, Evo2, Genos.

## Architecture

All models use the same lightweight MLP architecture:
- **Hidden layers**: [256, 128]
- **Activation**: ReLU
- **Optimizer**: Adam (lr=0.001, alpha=0.0001)
- **Training**: Early stopping with validation (10%)
- **Train/test split**: 80/20 (stratified, random_state=42)


## Input Format

All scripts expect NPZ files with the following structure:

```python
{
    'embeddings': np.array,  # shape: (N, embedding_dim)
    'labels': np.array,      # shape: (N,) - binary labels (0=normal, 1=splice-altering)
    'variant_ids': np.array  # shape: (N,) - variant identifiers
}
```
**Usage Example**:
```bash
python train_genos_mlp.py \
    --input_file <path_to_genos_embeddings.npz> \
    --model_file <path_to_save_model.joblib> \
    --scaler_file <path_to_save_scaler.joblib> \
    --metrics_file <path_to_save_metrics.json>
```

## Output

Each script outputs:
1. **Trained model** (`.joblib`) - MLP classifier
2. **Scaler** (`.joblib`) - StandardScaler for feature normalization
3. **Metrics** (`.json`) - Performance metrics:
   - AUROC, AUPRC
   - Accuracy, Precision, Recall, Specificity
   - F1-score, MCC
   - Confusion matrix
   - Dataset sizes
