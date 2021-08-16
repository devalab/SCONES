num_rounds=10
datetime=$(date '+%d-%m-%Y %H:%M:%S');
for dataset_tuple in "fireprotdb_cleaned_unique_filtered_easy_S768_SCONES S768_SCONES"
do
    set -- $dataset_tuple
    training_dataset=$1
    test_dataset=$2
    prefix="multicv-$datetime/$training_dataset"
    for r in $(seq $num_rounds)
    do
        python train_SCONES_CV.py --dir_prefix="$prefix" --train_dataset_path="data/$training_dataset.bin" &
    done
    wait
    python save_predictions.py --dataset_path="data/$test_dataset.bin" --results_path="output/SCONES/$prefix/$test_dataset.csv" --multicv_model_dir="output/SCONES/$prefix" 
done