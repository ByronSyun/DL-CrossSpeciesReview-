import numpy as np
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, average_precision_score, accuracy_score, 
                              confusion_matrix, precision_score, recall_score, f1_score, matthews_corrcoef)
import joblib
import argparse
import json


def main(args):
    data = np.load(args.input_file, allow_pickle=True)
    variant_ids = data['variant_ids']
    embeddings = data['embeddings']
    labels = data['labels']
    
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        embeddings, labels, variant_ids,
        test_size=0.2, random_state=42, stratify=labels
    )
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = MLPClassifier(
        hidden_layer_sizes=(256, 128),
        activation='relu',
        solver='adam',
        max_iter=1000,
        random_state=42,
        early_stopping=True,
        n_iter_no_change=20,
        verbose=True,
        learning_rate_init=0.001,
        alpha=0.0001,
        batch_size=64
    )
    
    model.fit(X_train_scaled, y_train)
    
    y_pred = model.predict(X_test_scaled)
    y_proba = model.predict_proba(X_test_scaled)[:, 1]
    
    auroc = roc_auc_score(y_test, y_proba)
    auprc = average_precision_score(y_test, y_proba)
    accuracy = accuracy_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)
    
    tn, fp, fn, tp = cm.ravel()
    precision = precision_score(y_test, y_pred, zero_division=0)
    recall = recall_score(y_test, y_pred, zero_division=0)
    f1 = f1_score(y_test, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_test, y_pred)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"Specificity: {specificity:.4f}")
    print(f"F1: {f1:.4f}")
    print(f"MCC: {mcc:.4f}")
    
    joblib.dump(model, args.model_output_file)
    joblib.dump(scaler, args.model_output_file.replace('.joblib', '_scaler.joblib'))
    
    metrics = {
        'auroc': float(auroc),
        'auprc': float(auprc),
        'accuracy': float(accuracy),
        'precision': float(precision),
        'recall': float(recall),
        'specificity': float(specificity),
        'f1_score': float(f1),
        'mcc': float(mcc),
        'confusion_matrix': cm.tolist(),
        'n_train': int(len(X_train)),
        'n_test': int(len(X_test))
    }
    
    with open(args.metrics_output_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    results_df = pd.DataFrame({
        'variant_id': ids_test,
        'true_label': y_test,
        'predicted_label': y_pred,
        'predicted_proba': y_proba
    })
    results_df.to_csv(args.proba_output_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train MLP on DNABERT-2 embeddings")
    parser.add_argument("--input_file", type=str, required=True, help="Input NPZ file")
    parser.add_argument("--model_output_file", type=str, required=True, help="Output model file")
    parser.add_argument("--metrics_output_file", type=str, required=True, help="Output metrics JSON")
    parser.add_argument("--proba_output_file", type=str, required=True, help="Output predictions TSV")
    args = parser.parse_args()
    main(args)
