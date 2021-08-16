DYNAMUT2_URL = "http://biosig.unimelb.edu.au/dynamut2/submit_prediction"

def fetch_dynamut2(driver, sample):
    driver.get(DYNAMUT2_URL)

    driver.find_element_by_id("pdb_accession_single").send_keys(sample['pdb_id'])
    driver.find_element_by_id("mutation_single").send_keys(sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"])
    driver.find_element_by_id("chain_single").send_keys(sample['chain_id'])

    import time
    time.sleep(1)
    driver.find_element_by_id("singlePredictionForm").find_element_by_xpath("//button[@type='submit']").click()

    job_url = driver.current_url

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.ID, "ppi2Prediction")),
            EC.presence_of_element_located((By.XPATH, "//span[contains(text(), 'Submission Error')]")),
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Error')]"))
        )
    )

    ddG = float(driver.find_element_by_id("ppi2Prediction").get_attribute('innerHTML').split()[0])
    return { "ddG_DynaMut2" : ddG, "url" : job_url }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/DynaMut2.csv", fetch_dynamut2)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/DynaMut2_S768.csv", fetch_dynamut2)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/DynaMut2_S350.csv", fetch_dynamut2)