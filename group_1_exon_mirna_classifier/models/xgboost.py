import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import classification_report, precision_recall_curve, average_precision_score
from scipy.sparse import hstack
from xgboost import XGBClassifier

def get_gc_content(seq_array):
    """Calculates GC content percentage for an array of sequences."""
    gc_contents = []
    for seq in seq_array:
        seq_upper = str(seq).upper()
        g_count = seq_upper.count('G')
        c_count = seq_upper.count('C')
        length = len(seq_upper) if len(seq_upper) > 0 else 1
        gc_contents.append((g_count + c_count) / length)
    return np.array(gc_contents).reshape(-1, 1)

print("Loading data...")
data_path = "data/combined_exon_mirna_2.csv"
df = pd.read_csv(data_path)

df['sequence'] = df['sequence'].astype(str)
df['label'] = df['label'].map({'exon': 0, 'miRNA': 1})

X_train_seq, X_test_seq, y_train, y_test = train_test_split(
    df['sequence'].values,
    df['label'].values,
    test_size=0.2,
    random_state=42,
    stratify=df['label'].values
)

print("Feature Engineering: Extracting TF-IDF k-mers & GC Content...")
vectorizer = TfidfVectorizer(analyzer='char', ngram_range=(2, 5), max_features=10000)
X_train_tfidf = vectorizer.fit_transform(X_train_seq)
X_test_tfidf = vectorizer.transform(X_test_seq)

# Calculate GC content
X_train_gc = get_gc_content(X_train_seq)
X_test_gc = get_gc_content(X_test_seq)

# Combine sparse k-mer matrix with dense GC content feature
X_train_features = hstack([X_train_tfidf, X_train_gc])
X_test_features = hstack([X_test_tfidf, X_test_gc])

print(f"Generated {X_train_features.shape[1]} total features (k-mers + GC).")

print("Starting Hyperparameter Tuning with RandomizedSearchCV...")

xgb_base = XGBClassifier(
    random_state=42,
    eval_metric='logloss',
    n_jobs=-1
)

param_grid = {
    'n_estimators': [100, 300, 500],
    'learning_rate': [0.01, 0.05, 0.1, 0.2],
    'max_depth': [4, 6, 8, 10],
    'scale_pos_weight': [3.0, 4.0, 5.0, 6.0],  # Tuning the imbalance penalty
    'subsample': [0.7, 0.8, 0.9, 1.0],         # Row sampling
    'colsample_bytree': [0.6, 0.8, 1.0],       # Feature sampling per tree
    'gamma': [0, 1, 5]                         # Minimum loss reduction for a split
}

search = RandomizedSearchCV(
    estimator=xgb_base,
    param_distributions=param_grid,
    n_iter=15,
    scoring='f1',       
    cv=3,
    verbose=2,
    random_state=42,
    n_jobs=1
)

search.fit(X_train_features, y_train)

best_model = search.best_estimator_
print(f"\nBest Parameters found:\n{search.best_params_}")

print("\nEvaluating Best Model...")

y_probs = best_model.predict_proba(X_test_features)[:, 1]

precisions, recalls, thresholds = precision_recall_curve(y_test, y_probs)

# Find the threshold that yields the best F1 Score
f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-10)
best_i = np.argmax(f1_scores)
best_threshold = thresholds[best_i] if best_i < len(thresholds) else 0.5
best_f1 = f1_scores[best_i]

print(f"\n=> Optimal Decision Threshold: {best_threshold:.4f} (Max F1: {best_f1:.4f})")

# Apply the optimal threshold
y_pred_optimal = (y_probs >= best_threshold).astype(int)

print("\n--- Classification Report (Optimized Threshold) ---")
print(classification_report(y_test, y_pred_optimal, target_names=["Exon (0)", "miRNA (1)"]))

# --- Plotting PR Curves and Thresholds ---
ap_score = average_precision_score(y_test, y_probs)

plt.figure(figsize=(14, 5))

# Plot 1: Precision-Recall Curve (y-axis: Precision, x-axis: Recall)
plt.subplot(1, 2, 1)
plt.plot(recalls, precisions, label=f'PR-AUC = {ap_score:.4f}', color='blue', linewidth=2)
# Highlight the optimal threshold point
plt.plot(recalls[best_i], precisions[best_i], 'ro', label=f'Best F1 (Thr={best_threshold:.2f})')
plt.xlabel('Recall (Sensitivity)')
plt.ylabel('Precision (PPV)')
plt.title('Precision-Recall Curve (XGBoost)')
plt.legend(loc='lower left')
plt.grid(True)

# Plot 2: Precision & Recall vs. Decision Threshold
plt.subplot(1, 2, 2)
plt.plot(thresholds, precisions[:-1], 'b--', label='Precision', linewidth=2)
plt.plot(thresholds, recalls[:-1], 'g-', label='Recall', linewidth=2)
# Draw vertical line for the best threshold
plt.axvline(x=best_threshold, color='r', linestyle=':', label='Optimal Threshold')
plt.xlabel('Decision Threshold')
plt.ylabel('Score')
plt.title('Precision and Recall vs. Threshold')
plt.legend(loc='lower left')
plt.grid(True)

plt.tight_layout()

# Save the plot
os.makedirs("plots", exist_ok=True)
plot_path = "plots/xgboost_tuned_pr_curves.png"
plt.savefig(plot_path)
print(f"\nPlots saved successfully to '{plot_path}'")
