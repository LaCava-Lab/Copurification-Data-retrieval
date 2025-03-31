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
def get_list(query, pmc_only = True):
    """Retrieve all PMIDs for a given query using the PubMedFetcher."""
    num_of_articles = 500
    start_index = 0
    pmids = []
    while True:
        pmid_batch = fetcher.pmids_for_query(query,
                                            retstart=start_index,
                                            retmax=num_of_articles,
                                            pmc_only=pmc_only)
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

@retry_on_communication_error()
def fetch_article(pmid: str) -> dict[str, str]:
    """Fetch a single article and return its data as a dict."""
    article = fetcher.article_by_pmid(pmid)
    return {
        "pmid": article.pmid,
        "pmc": article.pmc,
        "title": article.title,
        "journal": article.journal,
        "doi": article.doi,
        "issn": article.issn,
    }

def fetch_articles_to_dataframe(pmids: list[str], workers: int = 5) -> pd.DataFrame:
    """Fetch articles in parallel and return them as a DataFrame."""
    with ThreadPoolExecutor(max_workers=workers) as executor:
        articles_data = list(executor.map(fetch_article, pmids))
    
    return pd.DataFrame(articles_data)



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

def filter_oa_database(oa_file_list, pmc_id_list):
    """
    Filters based on the csv database list of PMCs that are available for full_text mining.

    Parameters:
    oa_file_list (csv): Filename of the CSV containing the OA file list.
    pmc_id_list (list): Filename of the list containing the PMC IDs.
    """
    # Read CSV file
    oa_file_list_df = pd.read_csv(oa_file_list)
    # Filter oa_database based on PMC ID list
    filtered_oa_database = oa_file_list_df[oa_file_list_df["Accession ID"].isin(pmc_id_list)]
    oa_pmcids = filtered_oa_database["Accession ID"]

    return oa_pmcids