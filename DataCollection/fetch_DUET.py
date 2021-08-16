DUET_URL = "http://biosig.unimelb.edu.au/duet/stability"

def fetch_DUET(driver, sample):
    driver.get(DUET_URL)

    if sample["id"] in ["322F", 21]:
        raise Exception("unsupported sample")

    driver.find_element_by_name("pdb_code").send_keys(sample['pdb_id'])
    driver.find_element_by_name("mutation").send_keys(sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"])
    driver.find_element_by_name("chain").send_keys(sample['chain_id'])

    import time
    time.sleep(1)
    driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/form/div[2]/div[1]/div/button").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Predicted Stability Change')]/following-sibling::font")),
            EC.presence_of_element_located((By.XPATH, "//h4[contains(text(), 'Error')]")),
            EC.presence_of_element_located((By.XPATH, "//body[contains(text(), 'Internal Server Error')]"))
        )
    )

    ddG_element = driver.find_element_by_xpath("//font[contains(text(), 'DUET')]/../following-sibling::font")
    ddG = float(ddG_element.get_attribute('innerHTML').split()[0])
    return { "ddG_DUET" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/DUET_Stransitive.csv", fetch_DUET)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/DUET_S768.csv", fetch_DUET)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/DUET_S350.csv", fetch_DUET)