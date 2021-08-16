# Dataset Generation

This directory contains tools to create datasets. We preprocess all published datasets and store them in a common consistent CSV format.

## Directory Structure

```
.
├── data
│   ├── FireProtDB
│   ├── LegacyPDB
│   ├── PDB
│   ├── RESSEQ
│   ├── S350
│   ├── S768
│   ├── Ssym
│   ├── Stransitive
│   ├── TrainingDatasets
│   └── UniProt
└── utils
```

The root directory contains standalone programs to curate datasets. These programs read/write data from the `data` directory. The `utils` directory contains utilities that are shared across multiple programs. `LegacyPDB`, `PDB`, `RESSEQ` and `UniProt` are used to cache data to improve processing speeds. The remaining subdirectories in `data` contain existing and newly curated datasets.

## Storage Format

UniProt sequences are stored in FASTA files in the `UniProt` directory. The `PDB` directory contains PBD structures in mmCIF format. The `RESSEQ` directory contains `.resseq` files which is a mapping of resseq ids to sequence positions (zero based indexing) in the amino acid sequence extracted from structure.

All the newly curated datasets are stored in a common CSV format. The column headers of the curated datasets consist of a fixed set of fields (such as `pdb_id`) plus additional dataset specific fields (mutation category in S768 for example). Every sample in the dataset is assigned a unique identifier which in combination with the dataset name can uniquely identify any sample. More information can be found in the documentation of `utils/DatasetGenerator.py`.

## Programs

### Test Datasets:

- `create_S350.py`
- `create_S768.py`
- `create_Ssym.py`
- `create_Stransitive.py`
- `create_FireProtDB.py`
- `create_TrainingDatasets.py`

`create_S350.py`, `create_S768.py`, `create_Ssym.py` and `create_Stransitive.py` create test datasets and store it in the common format.

### Base Training Datasts:

`create_FireProtDB.py` generates a variety of datasets with different augmentation options and order. These datasets serve as the base training dataset from which filtered training datasets are derived.

Data Processing Options:
- `cleaned`: removes exact duplicates (same id), samples with no structure, inconsistent samples
- `unique`: merges duplicate samples into one unique sample
- `symmetric`: augments dataset with reverse samples
- `transitive`: augments dataset with transitive samples

The order of operations in the dataset name is read left-to-right.

- `fireprotdb_cleaned.csv`
- `fireprotdb_cleaned_unique.csv`
- `fireprotdb_cleaned_unique_symmetric.csv`
- `fireprotdb_cleaned_unique_symmetric_transitive.csv`

Note that we add ddG values of two different samples while constructing transitive samples. This may lead to accumulation of errors that make the augmented training dataset unreliable. SCONES used `fireprotdb_cleaned_unique` only for these reasons.

### Filtered Training Datasets

`create_TrainingDatasets.py` creates filtered (removing samples similar to test set samples) datasets for each test set. Two training datasets are generated based on the two filtering criterions.