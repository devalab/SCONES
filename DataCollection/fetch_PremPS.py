def fetch_PremPS(driver, sample):
    driver.get("https://lilab.jysw.suda.edu.cn/research/PremPS/")
    driver.find_element_by_id("pdb_id_input").send_keys(sample['pdb_id'])
    driver.find_element_by_xpath("//button[@id='submit_button']").click()

    if sample["id"] in [330]:
        raise Exception("unsupported sample")


    chain_item = None
    for item in driver.find_elements_by_xpath("//h5[@id='True']"):
        if item.text.strip() == sample['chain_id']:
            print("Selected option", item.text, "for chain")
            chain_item = item
            break
    if not chain_item: raise Exception("Error: chain id not found")
    chain_item.click()

    import time
    time.sleep(1)
    driver.find_element_by_xpath("//button[@type='submit']").click()

    from selenium.webdriver.support.ui import Select
    select = Select(driver.find_element_by_xpath("//select[@name='chain_select1']"))
    select.select_by_index(1) # hard-coding this because we always have exactly one chain selected
    time.sleep(1)

    residue = Select(driver.find_element_by_xpath("//select[@name='res_select1']"))
    res_select1_option = None
    for option in residue.options:
        if (sample["ref_residue"] + ' ' + str(sample['resseq_position'])) in option.text:
            print("Selected option \"%s\" for ref. residue (required: %s%d)" % (option.text, sample["ref_residue"], sample["resseq_position"]))
            res_select1_option = option
            break
    if not res_select1_option:
        raise Exception("ref. residue and position not found")
    res_select1_option.click()
    time.sleep(1)

    residue = Select(driver.find_element_by_xpath("//select[@name='mut_select1']"))
    mut_select1_option = None
    for option in residue.options:
        if sample['mutated_residue'] == option.text.split()[0]:
            print("Selected option \"%s\" for mutated residue" % (option.text))
            mut_select1_option = option
            #break
    if not mut_select1_option:
        raise Exception("mutated residue not found")
    mut_select1_option.click()
    time.sleep(1)

    driver.find_element_by_id("submit_partners").click()
    time.sleep(1)

    job_url = driver.current_url
    print("URL:", job_url)

    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    WebDriverWait(driver, 86400).until(
        EC.presence_of_element_located((By.XPATH, "/html/body/div[2]/div/div[5]/div/table/tbody/tr/td[4]"))
    )

    ddG_element = driver.find_element_by_xpath("/html/body/div[2]/div/div[5]/div/table/tbody/tr/td[4]/b")
    ddG = float(ddG_element.get_attribute('innerHTML'))
    return { "ddG_PremPS" : ddG, "url" : job_url }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/PremPS_Stransitive.csv", fetch_PremPS)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/PremPS_S768.csv", fetch_PremPS)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/PremPS_S350.csv", fetch_PremPS)