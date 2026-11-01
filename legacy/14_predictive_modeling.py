#!/usr/bin/env python3
"""
Predictive Modeling - GPR153 and MGMT
Build models to predict survival outcomes based on gene expression
"""

import os
import sys
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, accuracy_score, classification_report
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_data():
    """Load expression and survival data"""
    # Load survival data
    surv_path = 'data/results/survival/gpr153_survival_data.csv'
    
    if not os.path.exists(surv_path):
        log_message("ERROR: Survival data not found. Run script 10 first.")
        return None, None
    
    log_message("Loading survival data...")
    surv_data = pd.read_csv(surv_path, index_col=0)
    
    # Load expression data
    expr_path = 'data/processed/combined_gbm_normal_expression.csv'
    expr_data = pd.read_csv(expr_path, index_col=0)
    
    # Filter GBM samples
    gbm_expr = expr_data[expr_data['group'] == 'GBM'].drop('group', axis=1)
    
    # Extract patient IDs
    def extract_patient_id(sample_id):
        parts = str(sample_id).split('-')
        if len(parts) >= 3:
            return f"{parts[0]}-{parts[1]}-{parts[2]}"
        return sample_id
    
    surv_data.index = [extract_patient_id(sid) for sid in surv_data.index]
    gbm_expr.index = [extract_patient_id(sid) for sid in gbm_expr.index]
    
    # Remove duplicates
    surv_data = surv_data[~surv_data.index.duplicated(keep='first')]
    gbm_expr = gbm_expr[~gbm_expr.index.duplicated(keep='first')]
    
    # Align
    common = surv_data.index.intersection(gbm_expr.index)
    
    surv_data = surv_data.loc[common]
    gbm_expr = gbm_expr.loc[common]
    
    log_message(f"Aligned {len(common)} samples")
    
    return gbm_expr, surv_data

def prepare_features(expr_data, genes_of_interest):
    """Prepare feature matrix"""
    available_genes = [g for g in genes_of_interest if g in expr_data.columns]
    
    if len(available_genes) == 0:
        return None, []
    
    X = expr_data[available_genes].values
    
    return X, available_genes

def create_survival_labels(surv_data, threshold_days=365):
    """Create binary labels: short-term vs long-term survival"""
    # 1 = short survival (< threshold), 0 = long survival (>= threshold)
    labels = (surv_data['survival_days'] < threshold_days).astype(int)
    
    return labels.values

def train_random_forest(X_train, y_train, X_test, y_test):
    """Train Random Forest classifier"""
    rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
    rf.fit(X_train, y_train)
    
    # Predictions
    y_pred = rf.predict(X_test)
    y_pred_proba = rf.predict_proba(X_test)[:, 1]
    
    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred_proba)
    
    # Feature importance
    importance = rf.feature_importances_
    
    return rf, accuracy, auc, importance, y_pred, y_pred_proba

def train_logistic_regression(X_train, y_train, X_test, y_test):
    """Train Logistic Regression"""
    # Standardize features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    lr = LogisticRegression(random_state=42, max_iter=1000)
    lr.fit(X_train_scaled, y_train)
    
    # Predictions
    y_pred = lr.predict(X_test_scaled)
    y_pred_proba = lr.predict_proba(X_test_scaled)[:, 1]
    
    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred_proba)
    
    # Coefficients
    coefficients = lr.coef_[0]
    
    return lr, accuracy, auc, coefficients, y_pred, y_pred_proba

def main():
    """Main predictive modeling pipeline"""
    log_message("Starting predictive modeling")
    print("=" * 70)
    print("PREDICTIVE MODELING PIPELINE")
    print("=" * 70)
    
    # Load data
    expr_data, surv_data = load_data()
    if expr_data is None:
        return
    
    # Create output directory
    output_dir = 'data/results/ml'
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare features
    print("\n" + "=" * 70)
    print("FEATURE PREPARATION")
    print("=" * 70)
    
    genes_of_interest = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    
    X, feature_names = prepare_features(expr_data, genes_of_interest)
    
    if X is None:
        log_message("ERROR: No features available")
        return
    
    print(f"  Features: {', '.join(feature_names)}")
    print(f"  Samples: {X.shape[0]}")
    print(f"  Features: {X.shape[1]}")
    
    # Create labels (1-year survival threshold)
    y = create_survival_labels(surv_data, threshold_days=365)
    
    print(f"\n  Short-term survival (<1 year): {y.sum()} samples")
    print(f"  Long-term survival (>=1 year): {(1-y).sum()} samples")
    
    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    
    print(f"\n  Training set: {len(X_train)} samples")
    print(f"  Test set: {len(X_test)} samples")
    
    # Random Forest
    print("\n" + "=" * 70)
    print("RANDOM FOREST CLASSIFIER")
    print("=" * 70)
    
    rf_model, rf_acc, rf_auc, rf_importance, rf_pred, rf_proba = train_random_forest(
        X_train, y_train, X_test, y_test
    )
    
    print(f"\n  Test Accuracy: {rf_acc:.4f}")
    print(f"  Test AUC: {rf_auc:.4f}")
    
    print(f"\n  Feature Importance:")
    for gene, imp in zip(feature_names, rf_importance):
        print(f"    {gene}: {imp:.4f}")
    
    # Save feature importance
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': rf_importance
    }).sort_values('importance', ascending=False)
    importance_df.to_csv(f'{output_dir}/rf_feature_importance.csv', index=False)
    
    # Logistic Regression
    print("\n" + "=" * 70)
    print("LOGISTIC REGRESSION")
    print("=" * 70)
    
    lr_model, lr_acc, lr_auc, lr_coef, lr_pred, lr_proba = train_logistic_regression(
        X_train, y_train, X_test, y_test
    )
    
    print(f"\n  Test Accuracy: {lr_acc:.4f}")
    print(f"  Test AUC: {lr_auc:.4f}")
    
    print(f"\n  Feature Coefficients:")
    for gene, coef in zip(feature_names, lr_coef):
        print(f"    {gene}: {coef:.4f}")
    
    # Save coefficients
    coef_df = pd.DataFrame({
        'feature': feature_names,
        'coefficient': lr_coef
    }).sort_values('coefficient', ascending=False)
    coef_df.to_csv(f'{output_dir}/lr_coefficients.csv', index=False)
    
    # Model comparison
    print("\n" + "=" * 70)
    print("MODEL COMPARISON")
    print("=" * 70)
    
    comparison = pd.DataFrame({
        'model': ['Random Forest', 'Logistic Regression'],
        'accuracy': [rf_acc, lr_acc],
        'auc': [rf_auc, lr_auc]
    })
    
    print(comparison.to_string(index=False))
    
    comparison.to_csv(f'{output_dir}/model_comparison.csv', index=False)
    
    # Save predictions
    predictions_df = pd.DataFrame({
        'true_label': y_test,
        'rf_prediction': rf_pred,
        'rf_probability': rf_proba,
        'lr_prediction': lr_pred,
        'lr_probability': lr_proba
    })
    predictions_df.to_csv(f'{output_dir}/model_predictions.csv', index=False)
    
    print(f"\n✓ Saved: rf_feature_importance.csv")
    print(f"✓ Saved: lr_coefficients.csv")
    print(f"✓ Saved: model_comparison.csv")
    print(f"✓ Saved: model_predictions.csv")
    
    print("\n" + "=" * 70)
    log_message("Predictive modeling complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
