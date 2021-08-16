IMUTANT3_URL = "http://gpcr.biocomp.unibo.it/cgi/predictors/I-Mutant3.0/I-Mutant3.0.cgi"

def fetch_IMutant3(driver, sample):
    driver.get(IMUTANT3_URL)

    driver.find_element_by_xpath("/html/body/center/form/p/table/tbody/tr[4]/td/div/input").click()
    driver.find_element_by_xpath("/html/body/center/form/p/table/tbody/tr[9]/td/div/input").click()
    
    import time
    time.sleep(1)

    if sample["id"] in [332]:
        raise Exception("unsupported sample")

    driver.find_element_by_name("proteina").send_keys(sample['pdb_id'])
    # driver.find_element_by_name("chain").send_keys(sample['chain_id'])
    driver.find_element_by_name("posizione").send_keys(sample['resseq_position'])
    driver.find_element_by_name("newres").send_keys(sample['mutated_residue'])

    import numpy as np
    if 'pH' in sample and not np.isnan(sample["pH"]):
        driver.find_element_by_name("ph").clear()
        driver.find_element_by_name("ph").send_keys(str(sample["pH"]))

    if 'T' in sample and not np.isnan(sample["T"]):
        driver.find_element_by_name("temp").clear()
        driver.find_element_by_name("temp").send_keys(str(sample["T"]))

    driver.find_element_by_name("submit").click()

    time.sleep(1)
    driver.find_element_by_xpath("/html/body/table[2]/tbody/tr/td[2]/table/tbody/tr/td[2]/div/font/a").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//pre[contains(text(), 'DDG Value Prediction')]")),
        )
    )

    ddG_element = driver.find_element_by_xpath("/html/body/table[2]/tbody/tr/td[2]/table/tbody/tr/td[2]/div/font/table[3]/tbody/tr[3]/td/font/b/pre")
    
    print(ddG_element.get_attribute('innerHTML').split()[3])
    ddG = float(ddG_element.get_attribute('innerHTML').split()[3])
    return { "ddG_IMutant3" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/IMutant3_Stransitive.csv", fetch_IMutant3)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/IMutant3_S768.csv", fetch_IMutant3)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/IMutant3_S350.csv", fetch_IMutant3)