import os
import logging
from pathlib import Path
import pandas as pd
import numpy as np

import xgboost as xgb
from sklearn.model_selection import StratifiedKFold
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
    logger = logging.getLogger(__name__)
    logger.info(f"Loading data from {filepath}...")
    df = pd.read_csv(filepath)
    feature_cols = [c for c in df.columns if c.startswith('esm_dim_') or c.startswith('kmer_')]
    X = df[feature_cols]
    y = df['label'].map({'Positive': 1, 'Negative': 0}).values
    return df, X, y, feature_cols

def run_stratified_kfold(X, y, params, output_dir, n_splits=5):
    logger = logging.getLogger(__name__)
    logger.info(f"Starting {n_splits}-Fold Stratified Cross-Validation...")
    
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    fold_metrics = []
    
    pos_count = np.sum(y == 1)
    neg_count = np.sum(y == 0)
    scale_pos_weight = neg_count / pos_count if pos_count > 0 else 1.0
    
    fold_params = params.copy()
    fold_params['scale_pos_weight'] = scale_pos_weight

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
        X_train, y_train = X.iloc[train_idx], y[train_idx]
        X_test, y_test = X.iloc[test_idx], y[test_idx]
        
        model = xgb.XGBClassifier(**fold_params, random_state=42)
        model.fit(X_train, y_train)
        
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test)
        
        metrics = {'fold': fold, 'n_samples': len(y_test)}
        metrics['accuracy'] = accuracy_score(y_test, y_pred)
        metrics['roc_auc'] = roc_auc_score(y_test, y_prob)
        metrics['pr_auc'] = average_precision_score(y_test, y_prob)
        metrics['precision'] = precision_score(y_test, y_pred, zero_division=0)
        metrics['recall'] = recall_score(y_test, y_pred, zero_division=0)
        metrics['f1'] = f1_score(y_test, y_pred, zero_division=0)
            
        fold_metrics.append(metrics)
        logger.info(f"  Fold {fold} - Accuracy: {metrics['accuracy']:.4f}, AUC: {metrics['roc_auc']:.4f}")
        
    metrics_df = pd.DataFrame(fold_metrics)
    
    # Calculate means
    mean_metrics = metrics_df.mean().to_dict()
    logger.info("\n--- Overall 5-Fold CV Performance ---")
    logger.info(f"Mean Accuracy: {mean_metrics['accuracy']:.4f}")
    logger.info(f"Mean ROC AUC:  {mean_metrics['roc_auc']:.4f}")
    logger.info(f"Mean F1 Score: {mean_metrics['f1']:.4f}")
    logger.info("-------------------------------------\n")
    
    reports_dir = Path(CONFIG['reports_dir'])
    reports_dir.mkdir(exist_ok=True, parents=True)
    metrics_file = reports_dir / "stratified_kfold_cv_results.csv"
    metrics_df.to_csv(metrics_file, index=False)
    logger.info(f"Saved 5-Fold CV results to {metrics_file}")
    
    return metrics_df

def main():
    setup_logging()
    logger = logging.getLogger(__name__)
    
    train_file = Path(CONFIG['output_dir']) / "features_combined_training.csv"
    df_train, X_train, y_train, feature_cols = load_data(str(train_file))
    
    xgb_params = {
        'max_depth': 4,
        'learning_rate': 0.05,
        'n_estimators': 150,
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'eval_metric': 'logloss'
    }
    
    run_stratified_kfold(X_train, y_train, xgb_params, Path(CONFIG['output_dir']))

if __name__ == "__main__":
    main()
