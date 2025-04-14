import logging
import time
import numpy as np
import pandas as pd
import requests
import json
import csv
import os as os
from tqdm import tqdm  
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type, stop_after_attempt, wait_exponential
from metapub import PubMedFetcher, PubMedArticle, pubmedcentral
from multiprocessing.pool import ThreadPool
from typing import List, Iterator, Optional
from habanero import Crossref
from concurrent.futures import ThreadPoolExecutor, as_completed

# Decorator 1 
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(0.4),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)
##initialize pubmedfetcher 
fetcher = PubMedFetcher()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_query_from_file(filename):
    """Read and clean query from a file."""
    try:
        with open(filename, 'r') as file:
            query = file.read()
        return query
    except FileNotFoundError:
        return ValueError(f"The file '{filename}' was not found.")
    except Exception as e:
        return ValueError(f"An error occurred while reading the query file: {e}")

@retry_on_communication_error()
def get_list(query, full_text):
    """Retrieve all PMIDs for a given query using the PubMedFetcher.
    query : str
    full_text: bool"""
    num_of_articles = 500
    start_index = 0
    pmids = []
    while True:
        pmid_batch = fetcher.pmids_for_query(query,
                                            retstart=start_index,
                                            retmax=num_of_articles,
                                            pmc_only = full_text)
        pmids.extend(pmid_batch)
        start_index = len(pmids)
        if len(pmid_batch) < num_of_articles:
            break
    return pmids

@retry_on_communication_error()
def fetch_pmids_over_period(query_file, start="2000-01-01", stop_date=None, full_text=False):
    """Fetch PMIDs over a specified period."""
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read query.")
        return np.array([])

    if stop_date is None:
        stop_date = datetime.now().strftime("%Y-%m-%d")

    start_date = date.fromisoformat(start)
    stop_date = date.fromisoformat(stop_date)
    pmid_list = []

    # Calculate total number of iterations (2-month chunks)
    total_months = (stop_date.year - start_date.year) * 12 + (stop_date.month - start_date.month)
    total_iterations = (total_months // 2) + 1  # Approximate count

    # Initialize progress bar
    with tqdm(total=total_iterations, desc="Fetching PMIDs", mininterval=0.1) as pbar:
        while True:
            next_start = start_date + relativedelta(months=2)
            end_date = next_start - relativedelta(days=1)

            date_str = f'''(("{start_date.strftime('%Y-%m-%d')}"[Date - Publication] : "{end_date.strftime('%Y-%m-%d')}"[Date - Publication]) '''
            pmids = get_list(date_str + query, full_text=full_text)
            pmid_list.extend(pmids)

            # Update progress bar
            pbar.update(1)
            pbar.set_postfix({"PMIDs": len(pmid_list)})

            start_date = next_start
            if next_start >= stop_date:
                break

    pmid_clean_list = list(set(pmid_list))
    logging.info(f"Total PMIDs fetched: {len(pmid_clean_list)}")
    return np.array(pmid_clean_list)


@retry_on_communication_error
def fetch_pmcid(pmid):
    """converts PMID to it corresponding PMCID."""

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
    """converts a list of PMIDs to a list of corresponding PMCIDs in parallel using multithreading."""
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

def filter_oa_database(pmc_id_list):
    """
    Filters based on database list of PMCs that are available for full_text mining.

    Parameter:
    pmc_id_list (list): Filename of the list containing the PMC IDs.
    """
    # Read open access database file from NCBI server \size: ~230 mgb\
    oa_db_url = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.txt"

    # Read directly into pandas
    oa_file_list_df = pd.read_csv(oa_db_url, sep="\t", header=None, 
                    names=['File', 'Publication_info', 'PMCID', 'PMID', 'License'])

    # Filter oa_database based on PMC ID list
    filtered_oa_database = oa_file_list_df[oa_file_list_df["PMCID"].isin(pmc_id_list)]
    oa_pmcids = filtered_oa_database["PMCID"]

    return oa_pmcids



def fetch_article(pmid: str) -> PubMedArticle:
    """Fetch a single article from PubMed by PMID."""
    article = fetcher.article_by_pmid(pmid)
    if article.pmid != pmid:
        logging.warning("Article with pmid=%r returned pmid=%r", pmid, article.pmid)
    return article

@retry_on_communication_error()
def fetch_articles(pmids: List[str], *, processes: Optional[int] = None) -> Iterator[PubMedArticle]:
    """Fetch multiple articles from PubMed in parallel using a thread pool with progress bar."""
    with ThreadPool(processes=processes) as pool:
        # Wrap with tqdm for progress tracking
        for article in tqdm(
            pool.imap_unordered(fetch_article, pmids),
            total=len(pmids),
            desc="Fetching PubMed articles",
            unit="article"
        ):
            if article is not None:
                yield article

def fetch_articles_meta(pmids: List[str]) -> pd.DataFrame:
    """Fetch articles and return them as a pandas DataFrame with progress tracking."""
    articles_data = []
    
    # Initialize progress bar for fetching articles
    for article in tqdm(
        fetch_articles(pmids, processes=5),
        total=len(pmids),
        desc="Fetching articles",
        unit="row"
    ):
        articles_data.append({
            'pmid': article.pmid,
            'pmc': article.pmc,
            'title': article.title,
            'journal': article.journal,
            'doi': article.doi,
            'issn': article.issn
        })
    
    return pd.DataFrame(articles_data)

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

def process_publishers(articles_df: pd.DataFrame, email: str) -> pd.DataFrame:
    """Process publisher information and return updated DataFrame."""
    # Remove duplicates based on title
    articles_df = articles_df.drop_duplicates(subset='pmid', keep='first')

    dois = articles_df['doi'].tolist() # Prepare lists to be iterated in the /publisher_crossref_doi/ function
    issns = articles_df['issn'].tolist()
    uids = articles_df['pmid'].tolist()
    # Get publishers using DOIs
    publisher_list = publisher_crossref_doi(dois, issns, uids, email)
    articles_df.loc[:, 'publisher'] = publisher_list
    
    # Identify missing publisher values
    missing_df = articles_df[articles_df['publisher'].isnull()]  

    # Use ISSNs to find missing publishers
    if not missing_df.empty:
        publisher_ids_from_issn = get_publisher_ids_from_issn(missing_df, email)
        new_missing_df = missing_df.copy()  # Avoid SettingWithCopyWarning
        new_missing_df.loc[:, 'publisher'] = publisher_ids_from_issn
        
        # Update the original DataFrame with new publisher data
        articles_df.update(new_missing_df)
    
    return articles_df
