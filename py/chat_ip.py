# %%chat-ip

import logging
import time
import numpy as np
import pandas as pd
import requests
import json
import csv
import os as os
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type, stop_after_attempt, wait_exponential
from metapub import PubMedFetcher, PubMedArticle, pubmedcentral
from multiprocessing.pool import ThreadPool
from typing import List, Iterator, Optional
from habanero import Crossref
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
file_handler = logging.FileHandler('pipe.log', mode='w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logging.getLogger().addHandler(file_handler)

# Retry logic for handling communication errors
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(0.4),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)

# Initialize the fetcher
fetcher = PubMedFetcher()


def setup_logging(filename='log.log', level=logging.INFO):

    file_handler = logging.FileHandler(filename, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(level)
    logger.addHandler(file_handler)


def setup_retry(max_delay=10, wait_time=0.4, exception_types=(Exception,)):
    partial(
        retry,
        stop=stop_after_delay(max_delay),  # Maximum wait time in seconds.
        wait=wait_fixed(wait_time),  # Wait time between retries in seconds.
        retry=retry_if_exception_type(exception_types)
    )


@retry_on_communication_error()
def fetch_article(pmid: str) -> PubMedArticle:
    """Fetch a single article from PubMed by PMID."""
    article = fetcher.article_by_pmid(pmid)
    if article.pmid != pmid:
        logging.warning("Article with pmid=%r returned pmid=%r", pmid, article.pmid)
    return article

@retry_on_communication_error()
def fetch_articles(pmids: List[str], *, processes: Optional[int] = None) -> Iterator[PubMedArticle]:
    """Fetch multiple articles from PubMed in parallel using a thread pool."""
    with ThreadPool(processes=processes) as pool:
        for article in pool.imap_unordered(fetch_article, pmids):
            if article is not None:
                yield article

def save_articles_to_csv(pmids: List[str], csv_file: str):
    """Fetch articles and save them to a CSV file."""
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        out = csv.writer(f)
        first = True
        for article in fetch_articles(pmids, processes=5):
            row = dict(
                pmid=article.pmid,
                pmc=article.pmc,
                title=article.title,
                journal=article.journal,
                doi=article.doi,
                issn=article.issn
            )
            if first:
                out.writerow(row.keys())
                first = False
            out.writerow(row.values())


def publisher_crossref_doi(dois, issns, uids, email):
    """Fetch publishers for a list of DOIs using Crossref, fallback to ISSN if DOI is None."""
    publishers = []
    cr = Crossref(mailto=email)
    
    for doi, issn, pmid in zip(dois, issns, uids):
        try:
            if pd.isna(doi):  # Check if DOI is NaN
                if issn:  # Check if ISSN is not empty
                    logging.info(f"DOI is NaN for PMID {pmid}, attempting with ISSN {issn}.")
                    publisher = get_publisher_id_from_issn(issn, email)
                else:
                    logging.info(f"Both DOI and ISSN are empty for PMID {pmid}, skipping.")
                    publisher = None
                publishers.append(publisher)
                continue

            # Attempt to fetch publisher from DOI
            work = cr.works(ids=doi)
            publisher = work["message"].get("publisher") if work else None
            if publisher:
                publishers.append(publisher)
            else:
                raise ValueError(f"No publisher found for DOI{doi}, pmid:{pmid}")
        except Exception as e:
            logging.error(f"Failed for DOI {doi} with error {e}; attempting with ISSN {issn} if not empty.")
            # Fallback to ISSN if DOI fails and ISSN is not empty
            publisher = get_publisher_id_from_issn(issn, email) if issn else None
            publishers.append(publisher)
    
    return publishers


def identify_missing_values(df, column_name):
    """Identify rows with missing values in a specified column."""
    missing_values = df[column_name].isnull()
    return df[missing_values]

@retry_on_communication_error()
def get_publisher_id_from_issn(issn: str, email: str) -> str:
    """Query the CrossRef API for a single ISSN and return the publisher ID."""
    url = f"https://api.crossref.org/works?filter=issn:{issn}&select=publisher&mailto={email}"
    try:
        time.sleep(0.4)  # Sleep to avoid hitting API rate limits
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = json.loads(response.text)
        if "message" in data and "items" in data["message"] and data["message"]["items"]:
            first_item = data["message"]["items"][0]
            if isinstance(first_item, dict):
                publisher_id = list(first_item.values())[0]
                if publisher_id is None:
                    logging.info(f"No publisher ID found for ISSN {issn}")
                    return None
                return publisher_id
    except requests.RequestException as e:
        logging.info(f"Error: API request failed for Issn {issn}: {e}")
    return None


def get_publisher_ids_from_issn(missing_df: pd.DataFrame, email: str) -> list:
    """Read ISSNs from the missing_df dataframe, query the CrossRef API, and return a list of publisher IDs."""
    issns = missing_df['issn'].tolist()
    publisher_id_list = []
    for issn in issns:
        publisher_id = get_publisher_id_from_issn(issn, email)
        publisher_id_list.append(publisher_id)
        time.sleep(0.4)  # Sleep to avoid hitting API rate limits
    return publisher_id_list

def add_publishers_to_csv(input_csv: str, output_csv: str, email: str):
    """Add publisher information to CSV using DOI and ISSN."""
    # Load the CSV file
    df = pd.read_csv(input_csv)
    df = df.drop_duplicates(subset='title', keep='first')
    dois = df['doi'].tolist() # Prepare lists to be iterated in the /publisher_crossref_doi/ function
    issns = df['issn'].tolist()
    uids = df['pmid'].tolist()

    # Get publishers using DOIs
    publisher_list = publisher_crossref_doi(dois, issns, uids, email)
    df.loc[:, 'publisher'] = publisher_list

    # Find missing publishers
    missing_df = identify_missing_values(df, 'publisher')

    # Use ISSNs to find missing publishers
    if not missing_df.empty:
        publisher_ids_from_issn = get_publisher_ids_from_issn(missing_df, email)
        new_missing_df = missing_df.copy() #made a copy of missing_df to deal with SettingWithCopyWarning in Pandas
        new_missing_df.loc[:, 'publisher'] = publisher_ids_from_issn #adding value to the slice of this particular copy
        # Combine the original data with the new data
        df.update(new_missing_df)


    # Save the final DataFrame to a CSV file
    df.to_csv(output_csv, index=False)

def read_query_from_file(query_file):
    """Read and clean query from a file."""
    try:
        with open(query_file, 'r') as file:
            query = file.read()
        return query
    except FileNotFoundError:
        logging.error(f"The file '{query_file}' was not found.")
        return ""
    except Exception as e:
        logging.error(f"An error occurred while reading the query file: {e}")
        return ""

@retry_on_communication_error()
def get_list(query, pmc_only=False):
    """Retrieve all PMIDs for a given query using the PubMedFetcher."""
    num_of_articles = 500
    start_index = 0
    pmids = []
    while True:
        pmid_batch = fetcher.pmids_for_query(query,
                                            retstart=start_index,
                                            retmax=num_of_articles,
                                            pmc_only = pmc_only)
        pmids.extend(pmid_batch)
        start_index = len(pmids)
        if len(pmid_batch) < num_of_articles:
            break
    return pmids

@retry_on_communication_error()
def fetch_pmids_over_period(query_file, start="2000-01-01", stop_date=None):
    """Fetch PMIDs over a specified period using a query read from a file."""
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read query.")
        return np.array([])

    if stop_date is None:
        stop_date = datetime.now().strftime("%Y-%m-%d")

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
        if next_start >= date.fromisoformat(stop_date):
            break

    # Remove duplicates by converting to a set, then back to a list
    pmid_clean_list = list(set(pmid_list))
    logging.info(f"Total PMIDs fetched: {len(pmid_clean_list)}")

    return np.array(pmid_clean_list)

def save_pmids(pmid_array, directory="PMID_lists"):
    """Save PMIDs to both a text file and a NumPy binary file."""
    os.makedirs(directory, exist_ok=True)
    date_tag = datetime.now().isoformat()[:10]
    # File paths
    txt_file_path = os.path.join(directory, f'pmids_{date_tag}.txt')
    npy_file_path = os.path.join(directory, f'pmids_{date_tag}.npy')
    # Save PMIDs to txt
    np.savetxt(txt_file_path, pmid_array, fmt='%s', delimiter=",")
    logging.info(f"PMIDs saved to text file: {txt_file_path}")
    # Save PMIDs to npy
    np.save(npy_file_path, pmid_array)
    logging.info(f"PMIDs saved to binary file: {npy_file_path}")


def filter_pmc_full_text(oa_file_list, fetched_dataframe):
    """
    Filters based on the csv database list of PMCs that are available for full_text mining and writes new column to df.

    Parameters:
    oa_file_list (str): Filename of the CSV containing the OA file list.
    fetched_dataframe (str): Filename of the CSV containing the PMC IDs.
    """
    # Read CSV files
    oa_file_list_df = pd.read_csv(oa_file_list)
    pmc_df = pd.read_csv(fetched_dataframe)
    pmc_df['pmc'] = pmc_df['pmc'].apply(lambda x: 'PMC' + str(int(x)) if pd.notnull(x) else x)

    # Check oa_database based on PMC ID list and add new column statin true or false 
    pmc_df['full text available via PMC'] = pmc_df['pmc'].isin(oa_file_list_df["Accession ID"])
    # Open a text file to write
    with open('full_text_pmc.txt', 'w') as file:
        for item in pmc_df['full text available via PMC']:
            file.write(str(item) + '\n')
    # Save the filtered DataFrame to a new CSV file
    pmc_df.to_csv("final_full_text.csv", index=False)

    return pmc_df