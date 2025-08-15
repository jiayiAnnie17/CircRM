import os
import pandas as pd
import xgboost as xgb
import csv
from sklearn.ensemble import IsolationForest
import numpy as np
import argparse
import json
import yaml

# -------------------------
# Default p-value mapping table
# -------------------------
default_thresholds_dict = {
    'm1A': {0.001: 0.912, 0.01: 0.805, 0.05: 0.662},
    'm5C': {0.001: 0.939, 0.01: 0.804, 0.05: 0.620},
    'm6A': {0.001: 0.919, 0.01: 0.818, 0.05: 0.678},
}

# -------------------------
# Command line interface
# -------------------------
parser = argparse.ArgumentParser(description="Flexible RNA modification prediction (pvalue & likelihood)")

parser.add_argument("-o", "--output_dir", required=True, help="Output directory for results")
parser.add_argument("--models_dir", required=True, help="Directory containing all model files")
parser.add_argument("--tasks", nargs="+", required=True,
                    help="Tasks in format: mod:mode:threshold1,threshold2:input_dir. "
                         "Example: m1A:pvalue:0.001,0.01:/Feature_A/tmp/m1A "
                         "m5C:likelihood:5,10:/Feature_C/tmp/m5C")
parser.add_argument("--pvalue_config", help="YAML or JSON file containing p-value thresholds mapping")

args = parser.parse_args()

output_folder = args.output_dir
models_dir = args.models_dir
tasks = args.tasks

os.makedirs(output_folder, exist_ok=True)

# -------------------------
# oad p-value configuration file
# -------------------------
thresholds_dict = default_thresholds_dict.copy()
if args.pvalue_config:
    if args.pvalue_config.endswith(".json"):
        with open(args.pvalue_config, "r") as f:
            loaded = json.load(f)
    elif args.pvalue_config.endswith((".yaml", ".yml")):
        with open(args.pvalue_config, "r") as f:
            loaded = yaml.safe_load(f)
    else:
        raise ValueError("Unsupported config file format. Must be .json or .yaml")

    thresholds_dict = {
        mod: {float(k): v for k, v in val.items()}
        for mod, val in loaded.items()
    }

# -------------------------
# IsolationForest
# -------------------------
iso_forest = IsolationForest(contamination='auto', random_state=42)

def calculate_likelihood(A, B):
    return A / B if B > 0 else np.inf

# -------------------------
# Find model file for each modification
# -------------------------
def find_model_path(mod):
    for fname in os.listdir(models_dir):
        if mod in fname and fname.endswith((".model", ".json")):
            return os.path.join(models_dir, fname)
    raise FileNotFoundError(f"No XGBoost model file found for {mod} in {models_dir}")

# -------------------------
# Load XGBoost model (.model or .json)
# -------------------------
def load_xgb_model(model_path):
    ext = os.path.splitext(model_path)[1].lower()
    if ext not in [".model", ".json"]:
        raise ValueError(f"Unsupported XGBoost model format: {ext}. Must be .model or .json")
    model = xgb.XGBClassifier()
    model.load_model(model_path)
    return model

# -------------------------
# Parse tasks
# -------------------------
parsed_tasks = []
for task in tasks:
    parts = task.split(":")
    if len(parts) != 4:
        raise ValueError(f"Invalid task format: {task}. Must be mod:mode:thresholds:input_dir")
    mod, mode, ths_str, task_input_dir = parts
    thresholds = [float(x) for x in ths_str.split(",")]
    if not os.path.exists(task_input_dir):
        raise FileNotFoundError(f"Input directory for {mod} does not exist: {task_input_dir}")
    parsed_tasks.append((mod, mode, thresholds, task_input_dir))

# -------------------------
# Main loop
# -------------------------
for (mod, mode, thresholds, task_input_dir) in parsed_tasks:
    if mod not in thresholds_dict:
        raise ValueError(f"Unknown modification: {mod}")
    # Load model for this modification
    model_path = find_model_path(mod)
    model = load_xgb_model(model_path)

    for threshold in thresholds:
        # Determine threshold value
        if mode == "pvalue":
            if threshold not in thresholds_dict[mod]:
                raise ValueError(f"No p-value mapping for {mod} at threshold {threshold}")
            threshold_value = thresholds_dict[mod][threshold]
        elif mode == "likelihood":
            threshold_value = threshold
        else:
            raise ValueError(f"Invalid mode: {mode}")

        all_results = []

        # Iterate through all CSV files in the input directory
        for file_name in os.listdir(task_input_dir):
            if file_name.endswith(".csv"):
                data_file = os.path.join(task_input_dir, file_name)
                data = pd.read_csv(data_file, sep=',')

                # Separate original data and features
                original_data = data.iloc[:, :6]
                X = data.iloc[:, 6:41].astype(float).values

                # Detect outliers using IsolationForest
                outliers = iso_forest.fit_predict(X)
                X_cleaned = X[outliers == 1]

                # Predict probabilities on cleaned data
                y_pred_prob_cleaned = model.predict_proba(X_cleaned)
                predicted_classes = []

                for prob in y_pred_prob_cleaned:
                    if mode == "pvalue":
                        predicted_class = np.argmax(prob)
                        if predicted_class == 1 and prob[1] >= threshold_value:
                            predicted_classes.append(1)
                        else:
                            predicted_classes.append(0)
                    elif mode == "likelihood":
                        A = prob[1]
                        B = prob[0]
                        likelihood = calculate_likelihood(A, B)
                        if likelihood >= threshold_value:
                            predicted_classes.append(1)
                        else:
                            predicted_classes.append(0)

                # Reconstruct predictions including outliers
                y_pred = []
                cleaned_idx = 0
                for outlier in outliers:
                    if outlier == -1:
                        y_pred.append(0)
                    else:
                        y_pred.append(predicted_classes[cleaned_idx])
                        cleaned_idx += 1

                # Concatenate original data with predictions
                result_data = pd.concat([original_data, pd.Series(y_pred)], axis=1)
                all_results.append(result_data)

        # Combine all results and save to output file
        final_result = pd.concat(all_results, ignore_index=True)
        mode_short = "p" if mode == "pvalue" else "l"
        output_filename = f"{mod}_{mode_short}_{threshold}.tsv"
        output_path = os.path.join(output_folder, output_filename)
        final_result.to_csv(output_path, sep='\t', quoting=csv.QUOTE_NONE,
                            escapechar=' ', index=False, header=False)
        print(f"[{mod} | {mode} | threshold={threshold}] Saved to {output_path}")
