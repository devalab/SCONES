DEEPDDG_URL = "http://protein.org.cn/ddg.html"

def fetch_DeepDDG(driver, sample):
    driver.get(DEEPDDG_URL)

    driver.find_element_by_name("pdbdown").send_keys(sample['pdb_id'])
    driver.find_element_by_xpath(".//input[@type='radio' and @value='selemut']").click()
    driver.find_element_by_name("mutlist").send_keys(sample["chain_id"] + ' ' +
                                                     sample['ref_residue'] + ' ' +
                                                     str(sample["resseq_position"]) + ' ' +
                                                     sample["mutated_residue"])
    
    import time
    time.sleep(1)
    driver.find_element_by_xpath("//input[@type='submit']").click()

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from utils.any_of import any_of

    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//div[@id='result']/br")),
            EC.presence_of_element_located((By.XPATH, "//div[contains(text(), 'Error')]"))
        )
    )

    job_url = driver.find_element_by_xpath("//div[@id='result']/a").get_attribute("href")
    print("URL:", job_url)

    driver.get(job_url)
    WebDriverWait(driver, 86400).until(
        any_of(
            EC.presence_of_element_located((By.XPATH, "//div[@id='waiting']/a[contains(text(), 'Download')]")),
            EC.presence_of_element_located((By.XPATH, "//div[contains(text(), 'Error')]"))
        )
    )
    link = driver.find_element_by_xpath("//div[@id='waiting']/a").click()

    ddG = float(driver.find_element_by_xpath("//pre").get_attribute("innerHTML").split()[-1])
    return { "ddG_DeepDDG" : ddG, "url" : job_url }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/DeepDDG_Stransitive.csv", fetch_DeepDDG)