MCSM_URL = "http://biosig.unimelb.edu.au/mcsm/stability"
LEGACY_PDB_DUMP="../DataGeneration/data/LegacyPDB/"

def fetch_mCSM(driver, sample):
    driver.get(MCSM_URL)

    import os
    pdb_path = os.path.join(LEGACY_PDB_DUMP, "pdb" + sample["pdb_id"] + ".ent")

    form = driver.find_element_by_xpath("//form[@action=\"/mcsm/stability_prediction\"]")
    form.find_element_by_name("wild").send_keys(pdb_path)
    form.find_element_by_name("mutation").send_keys(sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"])
    form.find_element_by_name("chain").send_keys(sample['chain_id'])

    import time
    time.sleep(1)
    form.find_element_by_xpath("//form[@action=\"/mcsm/stability_prediction\"]//button[@type='submit']").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Predicted Stability Change')]/following-sibling::font")),
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Error')]"))
        )
    )

    job_url = driver.current_url

    ddG_element = driver.find_element_by_xpath("//h4[contains(text(), 'Predicted Stability Change')]/following-sibling::font")
    ddG = float(ddG_element.get_attribute('innerHTML').split()[0])
    return  { "ddG_mCSM" : ddG, "url" : job_url }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/mCSM_Stransitive.csv", fetch_mCSM)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/mCSM_S768.csv", fetch_mCSM)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/mCSM_S350.csv", fetch_mCSM)