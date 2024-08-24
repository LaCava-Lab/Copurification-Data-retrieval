# %%Query to full_text PMCs
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type
from eutils import EutilsNCBIError, EutilsRequestError

import os as os
import pandas as pd 
import numpy as np 

from metapub import PubMedFetcher 
from metapub import pubmedcentral
import logging
import time

from concurrent.futures import ThreadPoolExecutor, as_completed
import csv

#configure log
file_handler = logging.FileHandler('PMC.log', mode='w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logging.getLogger().addHandler(file_handler)

# Initialize the fetcher
fetcher = PubMedFetcher()


# Retry logic for handling communication errors
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(0.4),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)

def read_query_from_file(filename):
    """Read and clean query from a file."""
    try:
        with open(filename, 'r') as file:
            query = file.read()
        return query
    except FileNotFoundError:
        logging.error(f"The file '{filename}' was not found.")
        return ""
    except Exception as e:
        logging.error(f"An error occurred while reading the query file: {e}")
        return ""

@retry_on_communication_error()
def get_list(query):
    """Retrieve all PMIDs for a given query using the PubMedFetcher."""
    num_of_articles = 500
    start_index = 0
    pmids = []
    while True:
        pmid_batch = fetcher.pmids_for_query(query,
                                            retstart=start_index,
                                            retmax=num_of_articles,
                                            pmc_only=True)
        pmids.extend(pmid_batch)
        start_index = len(pmids)
        if len(pmid_batch) < num_of_articles:
            break
    return pmids

@retry_on_communication_error()
def fetch_pmids_over_period(query_file, start="2000-01-01", stop=None):
    """Fetch PMIDs over a specified period using a query read from a file."""
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read query.")
        return np.array([])

    if stop is None:
        stop = datetime.now().strftime("%Y-%m-%d")

    start_date_str = start
    pmid_list = []

    while True:
        if date.fromisoformat(start_date_str) <= date.fromisoformat("2002-07-01"):
            month_interval = 6
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2005-11-01"):
            month_interval = 5
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2009-11-01"):
            month_interval = 4
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2011-10-01"):
            month_interval = 3
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2023-01-01"):
            month_interval = 2
        else:
            month_interval = 4

        next_start = date.fromisoformat(start_date_str) + relativedelta(months=month_interval)
        end_date = (next_start - relativedelta(days=1))
        end_date_str = end_date.strftime('%Y-%m-%d')

        date_str = f'''(("{start_date_str}"[Date - Publication] : "{end_date_str}"[Date - Publication]) '''
        pmids = get_list(date_str + query)
        pmid_list.extend(pmids)
        start_date_str = next_start.strftime('%Y-%m-%d')
        if next_start >= date.fromisoformat(stop):
            break

    # Remove duplicates by converting to a set, then back to a list
    pmid_clean_list = list(set(pmid_list))
    logging.info(f"Total PMIDs fetched: {len(pmid_clean_list)}")

    return np.array(pmid_clean_list)

def save_pmids(pmid_array, directory="PMID_lists"):
    """Save PMIDs to both a text file and a NumPy binary file."""
    os.makedirs(directory, exist_ok=True)
    date_tag = datetime.n
    txt_file_path = os.path.join(directory, f'pmids_{date_tag}.txt')

    npy_file_path = os.path.join(directory, f'pmids_{date_tag}.npy')

    np.savetxt(txt_file_path, pmid_array, fmt='%s', delimiter=",")
    logging.info(f"PMIDs saved to text file: {txt_file_path}")

    np.save(npy_file_path, pmid_array)
    logging.info(f"PMIDs saved to binary file: {npy_file_path}")


@retry_on_communication_error
def fetch_pmcid(pmid):
    try:
        pmc = pubmedcentral.get_pmcid_for_otherid(pmid)
        return pmc
    except (CommunicationError, ConnectionError) as e:
        logging.error(f"Error: API request failed for {pmid}: {e}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error for {pmid}: {e}")
        return None

def get_pmcid_for_otherid(pmid_clean_list):
    PMCIDs = []
    with ThreadPoolExecutor(max_workers=10) as executor:  # Adjust max_workers based on needs
        future_to_pmid = {executor.submit(fetch_pmcid, pmid): pmid for pmid in pmid_clean_list}
        for future in as_completed(future_to_pmid):
            pmid = future_to_pmid[future]
            try:
                pmc = future.result()
                PMCIDs.append(pmc)
            except Exception as e:
                logging.error(f"Error processing PMID {pmid}: {e}")
                PMCIDs.append(None)
    return PMCIDs

def filter_oa_database(oa_file_list, pmc_ids_filename):
    """
    Filters based on the csv database list of PMCs that are available for full_text mining and writes them  to CSV and txt.

    Parameters:
    oa_file_list (str): Filename of the CSV containing the OA file list.
    pmc_ids_filename (str): Filename of the CSV containing the PMC IDs.
    """
    # Read CSV files
    oa_file_list_df = pd.read_csv(oa_file_list)
    pmc_ids_df = pd.read_csv(pmc_ids_filename)
    # Extract PMC ID list from the DataFrame
    pmc_id_list = pmc_ids_df.iloc[:, 0].tolist()
    # Filter oa_database based on PMC ID list
    filtered_oa_database = oa_file_list_df[oa_file_list_df["Accession ID"].isin(pmc_id_list)]
    # Open a text file to write
    with open('full_text_pmc.txt', 'w') as file:
        for item in filtered_oa_database["Accession ID"]:
            file.write(str(item) + '\n')
    # Save the filtered DataFrame to a new CSV file
    filtered_oa_database.to_csv("full_text_data.csv", index=False)

    return filtered_oa_database