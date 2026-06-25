import os
import json
import logging
from pathlib import Path
import pandas as pd
import numpy as np

import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import (
    roc_auc_score, 
    accuracy_score, 
    precision_score, 
    recall_score, 
    f1_score,
    average_precision_score
)

# Import custom settings
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.settings import CONFIG, setup_logging

def load_data(filepath: str):
    """Load the feature CSV file and separate features, labels, and groups."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading data from {filepath}...")
    
    df = pd.read_csv(filepath)
    
    # Identify feature columns
    feature_cols = [c for c in df.columns if c.startswith('esm_dim_') or c.startswith('kmer_')]
    logger.info(f"Identified {len(feature_cols)} feature columns.")
    
    X = df[feature_cols]
    
    # For training data, map labels
    if 'label' in df.columns and set(df['label'].dropna().unique()).issubset({'Positive', 'Negative'}):
        y = df['label'].map({'Positive': 1, 'Negative': 0}).values
    else:
        y = None
        
    return df, X, y, feature_cols

def run_lopgo_cv(X, y, groups, genotypes, accessions, params, output_dir):
    """Run Leave-One-Group-Out Cross Validation."""
    logger = logging.getLogger(__name__)
    logger.info("Starting Leave-One-P-Genotype-Out (LOPGO) Cross-Validation...")
    
    logo = LeaveOneGroupOut()
    fold_metrics = []
    
    # Calculate scale_pos_weight
    pos_count = np.sum(y == 1)
    neg_count = np.sum(y == 0)
    scale_pos_weight = neg_count / pos_count if pos_count > 0 else 1.0
    logger.info(f"Class distribution: {pos_count} Positive, {neg_count} Negative. scale_pos_weight: {scale_pos_weight:.2f}")
    
    # Copy params and add scale_pos_weight
    fold_params = params.copy()
    fold_params['scale_pos_weight'] = scale_pos_weight
    
    # To store out-of-fold predictions
    oof_predictions = []

    for train_idx, test_idx in logo.split(X, y, groups=groups):
        holdout_genotype = groups.iloc[test_idx].iloc[0]
        logger.info(f"Fold: Holding out {holdout_genotype} ({len(test_idx)} samples)")
        
        # Skip if holdout has only one class and it's a problem for metrics, 
        # but for LOPGO it's expected that a holdout might only be all human or all animal.
        # We can still compute predictions.
        
        X_train, y_train = X.iloc[train_idx], y[train_idx]
        X_test, y_test = X.iloc[test_idx], y[test_idx]
        
        model = xgb.XGBClassifier(**fold_params, random_state=42)
        model.fit(X_train, y_train)
        
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test)
        
        # Save OOF predictions
        for i, idx in enumerate(test_idx):
            oof_predictions.append({
                'accession': accessions.iloc[idx],
                'genotype': holdout_genotype,
                'true_label': y_test[i],
                'pred_prob': y_prob[i],
                'pred_class': y_pred[i]
            })
        
        # Compute metrics if both classes are present in the test set, else skip AUC
        metrics = {'genotype': holdout_genotype, 'n_samples': len(y_test)}
        
        # Accuracy
        metrics['accuracy'] = accuracy_score(y_test, y_pred)
        
        # If both classes are present, calculate robust metrics
        if len(np.unique(y_test)) > 1:
            metrics['roc_auc'] = roc_auc_score(y_test, y_prob)
            metrics['pr_auc'] = average_precision_score(y_test, y_prob)
            metrics['precision'] = precision_score(y_test, y_pred, zero_division=0)
            metrics['recall'] = recall_score(y_test, y_pred, zero_division=0)
            metrics['f1'] = f1_score(y_test, y_pred, zero_division=0)
        else:
            metrics['roc_auc'] = np.nan
            metrics['pr_auc'] = np.nan
            metrics['precision'] = np.nan
            metrics['recall'] = np.nan
            metrics['f1'] = np.nan
            
            # Custom logging for single-class holdout
            true_class = np.unique(y_test)[0]
            pred_acc = np.mean(y_pred == true_class)
            logger.info(f"  Holdout has only class {true_class}. Prediction accuracy on this class: {pred_acc:.2f}")
            
        fold_metrics.append(metrics)
        logger.info(f"  Accuracy: {metrics['accuracy']:.4f}")
        
    metrics_df = pd.DataFrame(fold_metrics)
    
    # Save OOF predictions
    oof_df = pd.DataFrame(oof_predictions)
    
    reports_dir = Path(CONFIG['reports_dir'])
    reports_dir.mkdir(exist_ok=True, parents=True)
    
    metrics_file = reports_dir / "lopgo_cv_results.csv"
    oof_file = reports_dir / "lopgo_oof_predictions.csv"
    
    metrics_df.to_csv(metrics_file, index=False)
    oof_df.to_csv(oof_file, index=False)
    logger.info(f"Saved LOPGO CV results to {metrics_file}")
    
    return metrics_df

def train_final_model(X, y, params, feature_cols):
    """Train the final XGBoost model on the entire dataset."""
    logger = logging.getLogger(__name__)
    logger.info("Training final XGBoost model on all training data...")
    
    pos_count = np.sum(y == 1)
    neg_count = np.sum(y == 0)
    scale_pos_weight = neg_count / pos_count if pos_count > 0 else 1.0
    
    final_params = params.copy()
    final_params['scale_pos_weight'] = scale_pos_weight
    
    model = xgb.XGBClassifier(**final_params, random_state=42)
    model.fit(X, y)
    
    # Extract feature importances
    importances = model.feature_importances_
    feat_imp_df = pd.DataFrame({
        'feature': feature_cols,
        'importance': importances
    }).sort_values('importance', ascending=False)
    
    reports_dir = Path(CONFIG['reports_dir'])
    feat_imp_df.to_csv(reports_dir / "feature_importances.csv", index=False)
    
    # Save model
    models_dir = Path(CONFIG['output_dir']) / "models"
    models_dir.mkdir(exist_ok=True, parents=True)
    model_path = models_dir / "xgboost_vp4_model.json"
    model.save_model(model_path)
    logger.info(f"Saved final model to {model_path}")
    
    return model

def evaluate_model(model, X_eval, eval_df):
    """Predict on the evaluation dataset."""
    logger = logging.getLogger(__name__)
    logger.info("Predicting on evaluation dataset...")
    
    y_prob = model.predict_proba(X_eval)[:, 1]
    
    eval_df_out = eval_df.copy()
    eval_df_out['predicted_probability'] = y_prob
    eval_df_out['predicted_class'] = (y_prob >= 0.5).astype(int)
    
    out_file = Path(CONFIG['output_dir']) / "evaluation_predictions.csv"
    eval_df_out.to_csv(out_file, index=False)
    logger.info(f"Saved evaluation predictions to {out_file}")
    
def main():
    setup_logging()
    logger = logging.getLogger(__name__)
    logger.info("=== Starting XGBoost Modeling Pipeline ===")
    
    # 1. Load Training Data
    train_file = Path(CONFIG['output_dir']) / "features_combined_training.csv"
    df_train, X_train, y_train, feature_cols = load_data(str(train_file))
    
    # 2. XGBoost Parameters
    xgb_params = {
        'max_depth': 4,
        'learning_rate': 0.05,
        'n_estimators': 150,
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'eval_metric': 'logloss'
    }
    logger.info(f"XGBoost Parameters: {xgb_params}")
    
    # 3. Run LOPGO
    groups = df_train['genotype']
    accessions = df_train['accession']
    metrics_df = run_lopgo_cv(X_train, y_train, groups, groups, accessions, xgb_params, Path(CONFIG['output_dir']))
    
    # Highlight P[6] performance specifically as requested
    p6_metrics = metrics_df[metrics_df['genotype'] == 'P[6]']
    if not p6_metrics.empty:
        logger.info("\n--- P[6] Genotype Specific Performance ---")
        for col in p6_metrics.columns:
            logger.info(f"{col}: {p6_metrics.iloc[0][col]}")
        logger.info("------------------------------------------\n")
    else:
        logger.warning("P[6] genotype not found in training groups.")
        
    # 4. Train Final Model
    model = train_final_model(X_train, y_train, xgb_params, feature_cols)
    
    # 5. Predict on Evaluation Data
    eval_file = Path(CONFIG['output_dir']) / "features_combined_evaluation.csv"
    df_eval, X_eval, _, _ = load_data(str(eval_file))
    
    evaluate_model(model, X_eval, df_eval)
    
    logger.info("=== XGBoost Modeling Pipeline Complete ===")

if __name__ == "__main__":
    main()
