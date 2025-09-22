# D2S_ph
---

# Phage–Host Interaction Prediction with \$D\_{2}^{\ast}\$ Features

This repository provides a full pipeline for predicting phage–host interactions using \$D\_{2}^{\ast}\$-based feature extraction and machine learning models. The framework includes **sequence preprocessing, feature extraction, model training, evaluation, and visualization**.

---

## 📂 Project Structure

```
.
├── Executor.py          # Execution pipeline for pairwise sequence analysis
├── Similarity.py        # \$D\_{2}^{\ast}\$ similarity functions
├── Data.py              # FullSequence class for sequence handling
├── run/                 # Output directory
│   ├── data_pn/         # Input sequence pairs (JSON)
│   └── Pn_1k/           # Results directory
├── feature_utils.py     # Statistical & structural feature extraction
├── train_models.py      # Model training and evaluation
└── README.md            # Project documentation
```

---

## ⚙️ Workflow

### 1. Sequence Preprocessing

* Input: sequence pairs in JSON (`seq_pairs.json`)
* Operations:

  * **Divide sequences** into windows (`W=400`) with stride (`S=400`).
  * **Prepare k-mer (k=5) and Markov model (order=1)** representations.

```bash
# Step 1: split sequences
python main.py divide

# Step 2: compute k-mer + Markov background
python main.py prepare

# Step 3: run $D_{2}^{\ast}$ similarity analysis
python main.py batch_analysis
```

Each sequence pair produces a JSON file containing window-level \$D\_{2}^{\ast}\$ scores.

---

### 2. Feature Extraction

From the score matrix \$M\$, the pipeline extracts **59-dimensional statistical and structural descriptors**:

* **Distributional statistics**: mean, std, quantiles, skewness, kurtosis, entropy, Gini coefficient.
* **Top-k values**: statistics on local maxima.
* **Threshold coverage**: absolute + z-score based.
* **Matrix structure**: SVD eigen-spectrum, nuclear norm, Frobenius norm.
* **Positional features**: mean/variance of hotspot indices.
* **Extended descriptors**: tail quantiles, ratios, KL divergence from normality, positional entropy.

> Example:

```python
features, names = aggregate_plus_enhanced3(score_matrix)
```

---

### 3. Model Training & Evaluation

* Dataset split: **90% training / 10% testing**, stratified by class.
* Models supported:

  * XGBoost
  * Random Forest
  * Gradient Boosting
  * AdaBoost
  * Extra Trees
  * SVM
  * Logistic Regression
  * Naive Bayes
  * KNN
  * MLP

```python
from train_models import models

# Fit models
for name, model in models.items():
    model.fit(X_train, y_train)
```

* Metrics computed:

  * Accuracy
  * Precision
  * Recall
  * F1-score
  * ROC-AUC (Train & Test)

---

### 4. Visualization

The script generates comprehensive plots for model comparison:

* **Radar chart** of multi-metric performance
* **Test accuracy bar chart**
* **Precision/Recall/F1 grouped bars**
* **AUC ranking**
* **ROC curves (Train vs Test)**

```python
# Example: ROC curve visualization
fpr, tpr, _ = roc_curve(y_test, y_test_proba)
plt.plot(fpr, tpr, label=f"AUC={auc:.3f}")
```

---

## 🚀 Quick Start

```bash
# 1. Prepare environment
pip install -r requirements.txt

# 2. Run sequence preprocessing
python main.py

# 3. Train models and evaluate
python train_models.py
```

---

## 📊 Example Results

* **XGBoost** achieves the best performance:

  * Accuracy: **0.85**
  * AUC: **0.899**
  * Balanced precision & recall

* Random Forest and Gradient Boosting show competitive results, while simpler models (Naive Bayes, Logistic Regression) underperform.

---

## 📧 Contact

* Maintainers: Tianlai Huang†, Mianzhi Dai†, Zhong Li\*, Hongyu Duan\*
* † Equal contribution
* * Corresponding authors: [dhy.scut@outlook.com](mailto:duanhongyu@scut.edu.cn)

---

