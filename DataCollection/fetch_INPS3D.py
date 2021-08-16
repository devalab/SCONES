INPS3D_URL = "https://inpsmd.biocomp.unibo.it/inpsSuite/default/index3D"
LEGACY_PDB_DUMP="../DataGeneration/data/LegacyPDB/"

def fetch_INPS3D(driver, sample):
    import tempfile
    mutation_file = tempfile.NamedTemporaryFile().name

    mutation = sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"]
    with open(mutation_file, "w") as f:
        f.write(mutation)
    
    import os
    pdb_path = os.path.join(LEGACY_PDB_DUMP, "pdb" + sample["pdb_id"] + ".ent")

    driver.get(INPS3D_URL)

    driver.find_element_by_id("subtab_structure").send_keys(pdb_path)
    driver.find_element_by_id("subtab_mutations").send_keys(mutation_file)
    driver.find_element_by_id("subtab_inps_chain").send_keys(sample['chain_id'])

    import time
    time.sleep(1)
    driver.find_element_by_xpath("//input[@type='submit']").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//span[contains(text(), 'Completed')]")),
            EC.presence_of_element_located((By.XPATH, "//span[contains(text(), 'Failed')]"))
        )
    )

    driver.find_element_by_xpath("//span[contains(text(), 'Completed')]")
    jobid = driver.find_element_by_xpath("//div[contains(text(), 'Jobid:')]/following-sibling::div").get_attribute('innerHTML')
    job_url = "https://inpsmd.biocomp.unibo.it/inpsSuite/default/display_results.html?jobid=" + jobid

    driver.get(job_url)
    ddG = float(driver.find_element_by_xpath("//span[contains(text(), '%s')]/../following-sibling::div/span" % mutation).get_attribute('innerHTML'))
    return { "ddG_INPS3D" : ddG, "url" : job_url }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/INPS3D_Stransitive.csv", fetch_INPS3D)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/INPS3D_S768.csv", fetch_INPS3D)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/INPS3D_S350.csv", fetch_INPS3D)