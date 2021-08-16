def fetch_SDM2(driver, sample):
    driver.get("http://marid.bioc.cam.ac.uk/sdm2/prediction")

    form = driver.find_element_by_xpath("//form[@action=\"/sdm2/run_prediction\"]")
    form.find_element_by_name("pdb_code").send_keys(sample['pdb_id'])
    form.find_element_by_name("mutation").send_keys(sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"])
    form.find_element_by_name("chain").send_keys(sample['chain_id'])

    import time
    time.sleep(1)
    form.find_element_by_xpath("//form[@action=\"/sdm2/run_prediction\"]//button[@type='submit']").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Predicted pseudo')]/following-sibling::font")),
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Error')]")),
            EC.presence_of_element_located((By.XPATH, "//body[contains(text(), 'Internal Server Error')]"))
        )
    )

    ddG_element = driver.find_element_by_xpath("//h4[contains(text(), 'Predicted pseudo')]/following-sibling::font")
    ddG = float(ddG_element.get_attribute('innerHTML').split()[0])
    return { "ddG_SDM2" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/SDM2_Stransitive.csv", fetch_SDM2)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/SDM2_S768.csv", fetch_SDM2)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/SDM2_S350.csv", fetch_SDM2)