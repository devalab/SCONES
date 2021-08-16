# SCONES

- `prepare_datasets.py`: preprocesses datasets and cache them for fast training
- `train_SCONES_CV.py`: trains a model for a given training set using K-fold cross-validation
- `save_predictions.py`: saves predictions to a CSV file
- `build_ensemble.sh`: script to perfrom N rounds of cross-validation and save ensemble predictions

## Installation

1. `conda create --name scones_training python=3.6`
2. `conda activate scones_training`
3. https://pytorch.org/get-started/locally/
4. `conda install numpy pandas`
5. `pip install quantiprot`
6. `pip install git+https://github.com/jonbarron/robust_loss_pytorch`
7. `conda install -c conda-forge scikit-learn`
8. `conda install -c conda-forge biopython`
9. `conda install -c salilab dssp`

## Reproduction

1. Run `prepare_datasets.py`
2. Edit `build_ensemble.sh` to use the correct datasets
3. Run `build_ensemble.sh`
4. You will find the weights and logs in the `output` directory