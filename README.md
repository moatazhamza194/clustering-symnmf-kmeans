# clustering-symnmf-kmeans
SymNMF and K-means clustering implemented in Python and C with a Python C-extension, including analysis with silhouette scores.

Implementation of **Symmetric Non-negative Matrix Factorization (SymNMF)** and **K-means** clustering in both **Python** and **C**. Includes a Python C-extension for performance and an analysis tool comparing the two methods with silhouette scores.

## ✨ Features
- SymNMF algorithm implemented in Python and C
- Python C-extension (`symnmfmodule.c`) for integration
- Functions to compute:
  - Similarity matrix
  - Diagonal degree matrix
  - Normalized similarity matrix
- Hard clustering assignment from SymNMF association matrix
- Comparison with K-means using **silhouette score** (`sklearn.metrics`)
- Clean build system (`Makefile`, `setup.py`)

## 🛠️ Tech Stack
- **Languages:** Python, C  
- **Libraries:** NumPy, scikit-learn (for evaluation), Matplotlib (optional plots)  
- **Build tools:** Makefile, Python C API, setup.py  

## 🚀 Usage

### Build
```bash
make
python3 setup.py build_ext --inplace
