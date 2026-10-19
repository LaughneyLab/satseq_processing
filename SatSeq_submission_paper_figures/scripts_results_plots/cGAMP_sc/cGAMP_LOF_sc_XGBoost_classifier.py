import xgboost as xgb
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import scanpy as sc

def train_and_evaluate_xgboost(adata, observation, params):
    # 1. Isolate mSTING_WT and STOP cells for training/testing
    train_cells = adata[adata.obs[observation].isin(['mSTING_WT', 'STOP'])].copy()
    
    # Map labels explicitly: mSTING_WT = 1, STOP = 0
    label_map = {'mSTING_WT': 1, 'STOP': 0}
    y_train_test = train_cells.obs[observation].map(label_map).values
    X_train_test = train_cells.X
    
    
    # 2. Split into train and test/eval sets
    X_train, X_test, y_train, y_test = train_test_split(
        X_train_test, y_train_test, stratify=y_train_test, test_size=0.2,
        random_state=42)
    X_test, X_eval, y_test, y_eval = train_test_split(
        X_test, y_test, stratify=y_test, test_size=0.5, random_state=42)
    
    # 3. Train the XGBoost model
    clf = xgb.XGBClassifier(**params, 
                            enable_categorical=False)
    clf.fit(
        X_train, y_train,
        eval_set=[(X_eval, y_eval)],
        verbose=10
    )
    
    # 4. Evaluate on the test set
    y_pred = clf.predict(X_test)
    y_pred_proba = clf.predict_proba(X_test)[:, 1]
    
    print(f"\nTest Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"Test AUC: {roc_auc_score(y_test, y_pred_proba):.4f}")

    # 5. Apply the model to the "Other" cells (not STOP or mSTING_WT)
    other_cells = adata[~adata.obs[observation].isin(['mSTING_WT',
                                                      'STOP'])].copy()
    
    if other_cells.n_obs > 0:
        X_other = other_cells.X
        other_pred = clf.predict(X_other)
        other_pred_proba = clf.predict_proba(X_other)[:, 1]
        
        # Add predictions back to the subset's obs for inspection
        other_cells.obs['xgboost_pred_label'] = other_pred
        other_cells.obs['xgboost_pred_proba'] = other_pred_proba
        
        print(f"\nApplied model to {other_cells.n_obs} 'Other' cells.")
    else:
        print("\nNo 'Other' cells found to predict on.")

    return clf, other_cells, y_test, y_pred

# --- Execution ---
# Set up XGBoost parameters
params = {"max_depth": 4, 
    "learning_rate": 0.05,
    "objective": "binary:logistic", 
    "eval_metric": "auc", 
    "n_estimators": 128}


# Collect each cell line's REF (mSTING_WT) vs. STOP test-set predictions so
# both confusion matrices can be drawn together as one PNG at the end.
confusion_results = {}

for cell_line in ['4T1','EO771']:

    # Load sc data
    adata = sc.read_h5ad("".join(["adata_",cell_line,"_cGAMP.h5ad"]))

    # Run pipeline with the correct metadata column
    trained_clf, adata, y_test, y_pred = train_and_evaluate_xgboost(
        adata, "variant_substitution", params)

    confusion_results[cell_line] = (y_test, y_pred)

    # 1. Extract the standard base variants
    adata.obs["variant"] = adata.obs["variant_substitution"].str.extract(
        r"([A-Za-z]+\d+)"
    )
    
    # 2. Fill in the special cases where the regex returns NaN
    special_cases = ["mSTING_WT", "dCTT", "STOP", "RRDD"]
    
    is_special = adata.obs["variant_substitution"].isin(special_cases)
    
    adata.obs.loc[is_special, "variant"] = adata.obs.loc[is_special,
                                                              "variant_substitution"]
    
    # 3. Create the summarized DataFrame
    # Group by base_variant and aggregate mean, SD, and binary counts
    df_summary = (
        adata.obs.groupby("variant")
        .agg(
            mean_xgboost_pred_proba=("xgboost_pred_proba", "mean"),
            sd_xgboost_pred_proba=("xgboost_pred_proba", "std"),
            count_label_1=(
                "xgboost_pred_label",
                lambda x: (x == 1).sum(),
            ),  # Count 1s
            count_label_0=(
                "xgboost_pred_label",
                lambda x: (x == 0).sum(),
            ),  # Count 0s
        )).sort_values(by="mean_xgboost_pred_proba",
                       ascending=True).reset_index()
                       
                       
    # save
    df_summary.to_csv("".join(["XGBoost_output_",cell_line,"_cGAMP.csv"]))


# --- Combined REF (mSTING_WT) vs. STOP confusion matrices, both cell lines,
#     one PNG, raw counts in all 4 quadrants (no normalization) ---
cell_lines = list(confusion_results.keys())
fig, axes = plt.subplots(1, len(cell_lines), figsize=(6 * len(cell_lines), 5))
if len(cell_lines) == 1:
    axes = [axes]

for ax, cell_line in zip(axes, cell_lines):
    y_test, y_pred = confusion_results[cell_line]
    ConfusionMatrixDisplay.from_predictions(
        y_test, y_pred,
        display_labels=["STOP", "mSTING_WT"],
        normalize=None,      # raw counts, not fractions
        values_format="d",   # display as integers
        ax=ax,
        colorbar=False,
    )
    ax.set_title(f"{cell_line} REF/STOP Test Set")

fig.suptitle("REF (mSTING_WT) vs. STOP Validation Confusion Matrices")
fig.tight_layout()
fig.savefig("confusion_matrices_REF_STOP_validation.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("\nWrote confusion_matrices_REF_STOP_validation.png")













