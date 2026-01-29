import numpy as np
import json
import argparse
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, average_precision_score, matthews_corrcoef,
    confusion_matrix
)
import joblib


def main(args):
    data = np.load(args.input_file, allow_pickle=True)
    
    variant_ids = data['variant_ids']
    embeddings = data['embeddings']
    labels = data['labels']
    
    y = np.array([1 if str(label).lower() == 'splice-altering' else 0 for label in labels])
    X = embeddings
    
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        X, y, variant_ids, test_size=0.2, random_state=42, stratify=y
    )
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    mlp = MLPClassifier(
        hidden_layer_sizes=(256, 128),
        activation='relu',
        solver='adam',
        max_iter=1000,
        random_state=42,
        early_stopping=True,
        n_iter_no_change=20,
        verbose=True
    )
    
    mlp.fit(X_train_scaled, y_train)
    
    y_pred = mlp.predict(X_test_scaled)
    y_pred_proba = mlp.predict_proba(X_test_scaled)[:, 1]
    
    cm = confusion_matrix(y_test, y_pred)
    TN, FP, FN, TP = cm[0][0], cm[0][1], cm[1][0], cm[1][1]
    specificity = float(TN / (TN + FP)) if (TN + FP) > 0 else 0.0

    metrics = {
        'auroc': float(roc_auc_score(y_test, y_pred_proba)),
        'auprc': float(average_precision_score(y_test, y_pred_proba)),
        'accuracy': float(accuracy_score(y_test, y_pred)),
        'precision': float(precision_score(y_test, y_pred)),
        'recall': float(recall_score(y_test, y_pred)),
        'specificity': specificity,
        'f1_score': float(f1_score(y_test, y_pred)),
        'mcc': float(matthews_corrcoef(y_test, y_pred)),
        'confusion_matrix': cm.tolist(),
        'n_train': len(X_train),
        'n_test': len(X_test),
        'embedding_dim': embeddings.shape[1]
    }
    
    print(f"AUROC: {metrics['auroc']:.4f}")
    print(f"AUPRC: {metrics['auprc']:.4f}")
    print(f"Accuracy: {metrics['accuracy']:.4f}")
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Recall: {metrics['recall']:.4f}")
    print(f"Specificity: {metrics['specificity']:.4f}")
    print(f"F1: {metrics['f1_score']:.4f}")
    print(f"MCC: {metrics['mcc']:.4f}")
    
    joblib.dump(mlp, args.model_file)
    joblib.dump(scaler, args.scaler_file)
    
    with open(args.metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train MLP on Genos embeddings")
    parser.add_argument("--input_file", type=str, required=True, help="Input NPZ file")
    parser.add_argument("--model_file", type=str, required=True, help="Output model file")
    parser.add_argument("--scaler_file", type=str, required=True, help="Output scaler file")
    parser.add_argument("--metrics_file", type=str, required=True, help="Output metrics JSON")
    
    args = parser.parse_args()
    main(args)
