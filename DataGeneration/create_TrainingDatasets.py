def create_training_dataset(base_dataset_path, test_dataset_path, save_path):
    import pandas as pd
    csv_data = pd.read_csv(base_dataset_path)
    columns = csv_data.columns.values.tolist()
    training_base_dataset = csv_data.T.to_dict()
    training_base_dataset = { sample["id"] : sample for sample in training_base_dataset.values()}

    csv_data = pd.read_csv(test_dataset_path)
    test_dataset = csv_data.T.to_dict()
    test_dataset = { sample["id"] : sample for sample in test_dataset.values()}

    from utils.filtering import is_similar_easy, is_similar_hard, remove_similar_samples
    dataset_easy = remove_similar_samples(training_base_dataset, test_dataset, is_similar_easy)
    dataset_hard = remove_similar_samples(training_base_dataset, test_dataset, is_similar_hard)

    import os
    trainingset_basename = os.path.splitext(os.path.basename(base_dataset_path))[0]
    testset_basename = os.path.splitext(os.path.basename(test_dataset_path))[0]

    def save(records, path):
        df = pd.DataFrame.from_dict(records.values())
        df = df.reindex(columns, axis=1)
        df.to_csv(path, index=False)

    print("Training set size:", len(training_base_dataset))
    print("Training set size after easy filtering:", len(dataset_easy))
    print("Training set size after hard filtering:", len(dataset_hard))
    save(dataset_easy, os.path.join(save_path, trainingset_basename + "_filtered_easy_" +  testset_basename + ".csv"))
    save(dataset_hard, os.path.join(save_path, trainingset_basename + "_filtered_hard_" +  testset_basename + ".csv"))

TRAINSET_CSV="data/FireProtDB/fireprotdb_cleaned_unique.csv"
create_training_dataset(TRAINSET_CSV, "data/S350/S350_SCONES.csv", "data/TrainingDatasets/")
create_training_dataset(TRAINSET_CSV, "data/Ssym/Ssym_SCONES.csv", "data/TrainingDatasets/")
create_training_dataset(TRAINSET_CSV, "data/S768/S768_SCONES.csv", "data/TrainingDatasets/")

TRAINSET_TRANSITIVE_CSV="data/FireProtDB/fireprotdb_cleaned_unique_symmetric_transitive.csv"
create_training_dataset(TRAINSET_TRANSITIVE_CSV, "data/S350/S350_SCONES.csv", "data/TrainingDatasets/")
create_training_dataset(TRAINSET_TRANSITIVE_CSV, "data/Ssym/Ssym_SCONES.csv", "data/TrainingDatasets/")
create_training_dataset(TRAINSET_TRANSITIVE_CSV, "data/S768/S768_SCONES.csv", "data/TrainingDatasets/")
