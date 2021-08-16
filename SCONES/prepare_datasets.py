datasets_lists = [
    "../DataGeneration/data/S350/S350_SCONES.csv",
    "../DataGeneration/data/Ssym/Ssym_SCONES.csv",
    "../DataGeneration/data/S768/S768_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_easy_S350_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_easy_S768_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_easy_Ssym_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_hard_S350_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_hard_S768_SCONES.csv",
    "../DataGeneration/data/TrainingDatasets/fireprotdb_cleaned_unique_filtered_hard_Ssym_SCONES.csv",
    "../DataGeneration/data/FireProtDB/fireprotdb_cleaned_unique.csv"
]

PDB_PATH="../DataGeneration/data/PDB"

for ds_path in datasets_lists:
    import os
    base = os.path.basename(ds_path)
    name = os.path.splitext(base)[0]

    from datasets.SCONES import SCONESDataset
    dataset = SCONESDataset()
    dataset.load_dataset(ds_path, PDB_PATH)
    dataset.serialize(os.path.join("data", name + ".bin"))
