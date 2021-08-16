import os
import numpy as np

import pandas as pd

class DataFetch:
    DRIVER = "./geckodriver"

    def fetch(self, source_dataset_csv_path, result_csv_path, fetchfn, web=True):
        source_dataset = pd.read_csv(source_dataset_csv_path)
        source_dataset = source_dataset.T.to_dict()
        source_dataset = { sample["id"] : sample for sample in source_dataset.values() }

        results = {}
        if os.path.exists(result_csv_path):
            results = pd.read_csv(result_csv_path)
            results = results.T.to_dict()
            results = { sample["id"] : sample for sample in results.values() }

        if web:
            from selenium import webdriver
            driver = webdriver.Firefox(executable_path=self.DRIVER)
        for idx, sample in source_dataset.items():
            print("Sample ID:", idx, flush=True)
            print("pdb_id: %s, chain_id: %s, ref. residue: %s, resseq_position: %d, mutated residue: %s, pH: %.2f, T: %.2f" % (sample["pdb_id"], sample["chain_id"], sample["ref_residue"], sample["resseq_position"], sample["mutated_residue"], sample.get("pH", np.nan), sample.get("T", np.nan)))
            if idx in results:
                print("[SKIP] already in results", flush=True)
                continue
            try:
                result_data = fetchfn(driver, sample) if web else fetchfn(sample)
                result = { "id" : idx, **result_data }
                results[idx] = result
                print(result, flush=True)
            except Exception as e:
                print("Error for", idx, flush=True)
                print(e, flush=True)
            finally:
                if len(results):
                    df = pd.DataFrame.from_dict(results.values())
                    df.to_csv(result_csv_path, index=False)
            print()
        if web:
            driver.close()
