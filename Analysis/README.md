# Analysis

This section contains standalone python programs to analyze a dataset or evaluate predictions. All the programs use CSV files as inputs for reference datasets and predictions. The format for datasets is described in the DataGeneration section.

**Examples:**
- `python analyze_dataset.py` - displays help
- `python analyze_dataset.py --dataset_path=../DataGeneration/data/Stransitive/Stransitive.csv`
- `python evaluate_symmetric_consistency.py --symmetric_pairs_list=../DataGeneration/data/Ssym/Ssym_SCONES_idx.csv --results_path=../DataCollection/results/METHODNAME_Stransitive.csv --results_fieldname=ddG_METHODNAME`
- `python evaluate_transitive_consistency.py --transitive_tuples_list=../DataGeneration/data/Stransitive/Stransitive_idx.csv --results_path=../DataCollection/results/METHODNAME_Stransitive.csv --results_fieldname=ddG_METHODNAME`

**Format of results CSV:**

The reference dataset CSV file contains an `id` field that uniquely identifies each sample in the dataset. The results CSV file must contain an `id` field and a prediction field. All the programs match samples in the dataset CSV files and the results CSV file using the `id` field.

```
id,ddG,somefield1,somefield2
1F,-0.33,someval,someval
1R,0.33,someval,someval
.
.
.
```

## Installation

1. `conda create --name scones_analysis python=3.6`
2. `conda activate scones_analysis`
3. `conda install numpy pandas seaborn`
4. `conda install -c schrodinger -c conda-forge pymol-bundle` (required by `analyze_structural_changes.py` only)