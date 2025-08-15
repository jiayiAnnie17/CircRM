import os
import pandas as pd
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import roc_curve, auc, classification_report, confusion_matrix
import argparse

def plot_roc_curve(fpr, tpr, label, plot_file_path):
    # Plot ROC curve and save to file
    plt.plot(fpr, tpr, label=label)
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--', label='Random')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.savefig(plot_file_path)
    plt.close()

def plot_roc_for_binary_class(y_test, y_pred_prob, plot_folder, mod_type):
    # Compute ROC curve and AUC for binary classification
    fpr, tpr, _ = roc_curve(y_test, y_pred_prob[:, 1])
    auc_score = auc(fpr, tpr)
    roc_plot_path = os.path.join(plot_folder, f"ROC_{mod_type}.png")
    plot_roc_curve(fpr, tpr, f"ROC curve (AUC = {auc_score:.3f})", roc_plot_path)

def train_and_evaluate(X, y, param_grid, output_dir, plot_folder, mod_type):
    print(f"---- Training Binary Classification Model for {mod_type} ----")

    # Split dataset into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    # Initialize XGBoost classifier with Logloss as evaluation metric
    model = XGBClassifier(
        objective="binary:logistic",
        use_label_encoder=False,
        eval_metric="logloss"
    )

    # Perform hyperparameter tuning using RandomizedSearchCV
    grid_search = RandomizedSearchCV(
        estimator=model,
        param_distributions=param_grid,
        cv=3, n_jobs=-1, n_iter=50,
        scoring='roc_auc'
    )

    # Fit the model with hyperparameter search
    grid_search.fit(X_train, y_train)
    best_model = grid_search.best_estimator_
    print("Best parameters: ", grid_search.best_params_)

    # Fit the best model on the training set
    best_model.fit(X_train, y_train)
    y_pred_prob = best_model.predict_proba(X_test)
    y_pred = best_model.predict(X_test)

    # Save trained model to file
    model_file_path = os.path.join(output_dir, f"{mod_type}.model")
    best_model.save_model(model_file_path)

    # Set class names dynamically based on modification type
    if mod_type in ['m1A', 'm6A']:
        mod_name = mod_type
        unmod_name = 'normalA'
    elif mod_type == 'm5C':
        mod_name = mod_type
        unmod_name = 'normalC'
    else:
        mod_name = mod_type
        unmod_name = 'unmod'

    # Print confusion matrix as a table
    conf_matrix = confusion_matrix(y_test, y_pred)
    print("\nConfusion Matrix (as table):")
    print(pd.DataFrame(conf_matrix, columns=[unmod_name, mod_name], index=[unmod_name, mod_name]))

    # Plot ROC curve
    plot_roc_for_binary_class(y_test, y_pred_prob, plot_folder, mod_type)

    # Print classification report
    report = classification_report(y_test, y_pred, target_names=[unmod_name, mod_name])
    print("Classification Report:")
    print(report)

    return best_model

def main():
    parser = argparse.ArgumentParser(description="Train XGBoost binary classifier for RNA modifications")
    parser.add_argument("-m", "--mod_file", required=True, help="Path to the mod (positive) CSV file")
    parser.add_argument("-u", "--unmod_file", required=True, help="Path to the unmod (negative) CSV file")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save the trained model")
    parser.add_argument("-p", "--plot_folder", required=True, help="Directory to save plots")
    parser.add_argument("-t", "--mod_type", required=True, help="Modification type (e.g., m5C, m6A)")
    args = parser.parse_args()

    # Create output directories if not exist
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.plot_folder, exist_ok=True)

    # Read CSV files (comma-separated)
    mod = pd.read_csv(args.mod_file, header=None)
    unmod = pd.read_csv(args.unmod_file, header=None)

    # Assign labels: mod=1, unmod=0
    mod['label'] = 1
    unmod['label'] = 0

    # Combine datasets
    data = pd.concat([mod, unmod], axis=0)
    X = data.iloc[:, 6:41].astype(float).values  # Extract features columns
    y = data['label'].values

    # Define hyperparameter grid for RandomizedSearchCV
    param_grid = {
        'reg_alpha': [0, 0.1, 1],
        'reg_lambda': [0, 1, 5],
        'max_depth': [3, 6, 9],
        'min_child_weight': [1, 5, 10],
        'learning_rate': [0.01, 0.05, 0.1],
        'n_estimators': [100, 500, 1000],
        'subsample': [0.5, 0.7, 0.8],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'gamma': [0, 0.1, 1],
    }

    train_and_evaluate(X, y, param_grid, args.output_dir, args.plot_folder, args.mod_type)

if __name__ == "__main__":
    main()
