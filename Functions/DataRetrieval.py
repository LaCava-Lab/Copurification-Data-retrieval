import email
import os
#######################################################################
import os
import requests
import logging
import time
from dateutil.relativedelta import relativedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
import matplotlib.pyplot as plt
import csv
import pandas as pd 
import numpy as np 
from tqdm.notebook import tqdm
from Bio import Entrez
from eutils import EutilsNCBIError, EutilsRequestError
from datetime import datetime
from tqdm.auto import tqdm
#Import DR module from Functions folder
from Functions import DataRetrieval as DR

import urllib
import json
import random  
from ratelimit import limits, sleep_and_retry

#######################################################################

from Reference_files.keys import NCBI_API_KEY as api_key

import re
from Bio import Entrez, Medline
import logging
import time
import numpy as np
import pandas as pd
import requests
import json
import csv
import os as os
import urllib.parse
from time import sleep
from tqdm import tqdm  
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type, stop_after_attempt, wait_exponential
from multiprocessing.pool import ThreadPool
from typing import List, Iterator, Optional, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import xml.etree.ElementTree as ET
from matplotlib.ticker import FuncFormatter
from Bio import Entrez
from urllib.error import HTTPError
import random

# Decorator 1 
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(2),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)


# Configure max retries on failed requests (optional, default is 3)
Entrez.max_tries = 10
#####Bio.Entrez.max_tries and Bio.Entrez.sleep_between_tries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###################################################33 33##################################################################################
######################################################## Get PMIDs from Query ##########################################################

def get_pubmed_count(query, full_text=False):
    """
    Fetches the estimated count of records matching a query using Entrez ESearch.

    Args:
        query (str): The search term.
        full_text (bool): If True, searches PMC database; otherwise, searches PubMed.

    Returns:
        int: The estimated count of matching records, or 0 on failure.
    """
    db = "pmc" if full_text else "pubmed"
    if not Entrez.email:
        Entrez.email = input("Enter your email for NCBI Entrez: ").strip()
    if not Entrez.api_key and API_KEY:
        Entrez.api_key = input("Enter your NCBI API key (or press Enter to skip): ").strip() or None

    try:
        # Perform the search with retmax=0 to only get the count.
        # Entrez.esearch automatically handles rate limiting and retries.
        handle = Entrez.esearch(db=db, term=query, retmax=0, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        # Extract the count from the parsed XML record.
        count_str = record.get("Count", "0")
        
        # Safely convert the count string to an integer.
        return int(count_str) if count_str.isdigit() else 0

    except Exception as e:
        # Log any unexpected errors during the process.
        logging.error(f"Error fetching count from Entrez for query '{query[:100]}...': {e}")
        return 0



def get_list(query, full_text=False, api_key=None):
    """Fetches a list of PMIDs matching a query using Entrez ESearch """
    db = "pmc" if full_text else "pubmed"
    # Ensure Entrez.email is set before calling
    if not Entrez.email:
        logging.error("Entrez.email must be set before calling get_list.")
        return []


    start_index = 0
    pmids = []
    max_retries = 5
    retmax = 9999 
    try:
        # Use retmax=0 to just get the count and WebEnv info if needed later
        handle_count = Entrez.esearch(db=db, term=query, retmax=0, retmode="xml") # api_key included if set globally or passed
        record_count = Entrez.read(handle_count)
        handle_count.close()
        total_count_str = record_count.get("Count", "0")
        total_count = int(total_count_str) if total_count_str.isdigit() else 0

        # --- NEW ---
        # Respect the Entrez 10,000 record limit for esearch
        # ESearch can only retrieve the first 10,000 records matching the query.
        total_count = min(total_count, 9999)  # Cap at 9999 to avoid issues
        # --- END NEW ---

        if total_count == 0:
            logging.info(f"No records found for query: {query[:100]}...")
            return []
    except Exception as e:
        logging.error(f"Error getting record count for query '{query[:100]}...': {e}")
        return []

    # Main fetching loop
    start_index = 0
    pmids = []
    # Loop while haven't fetched all records up to the capped total_count
    while start_index < total_count: 
        # Determine how many records to fetch in this batch
        # Ensure we don't try to fetch more than total_count allows
        current_retmax = min(retmax, total_count - start_index)

        for attempt in range(max_retries):
            try:
                handle = Entrez.esearch(
                    db=db,  
                    term=query,
                    retstart=start_index,
                    retmax=current_retmax, # Use the calculated batch size
                    retmode="xml",
                    api_key=api_key # Pass if needed
                )
                record = Entrez.read(handle)
                handle.close()

                current_batch = record.get("IdList", [])
                batch_size = len(current_batch)
                
                # --- IMPROVED LOGIC ---
                if batch_size == 0:
                    # It's safer to log and potentially break to avoid infinite loops.
                    logging.warning(f"Received empty batch at start_index {start_index} (expected up to {current_retmax} more). Stopping fetch.")
                    break # Exit the retry loop and the while loop
                
                pmids.extend(current_batch)
                start_index += batch_size # Increment by actual batch size fetched
                
                logging.debug(f"Fetched batch: {batch_size} PMIDs (Total so far: {len(pmids)})")

                # --- TERMINATION CONDITION ---
                # Check if we have fetched all expected records 
                # (either we hit the 10k limit or the actual total)
                if len(pmids) >= total_count: 
                    logging.debug("Reached estimated total count. Fetching complete.")
                    break # Exit retry loop, condition will also exit while loop

                break # Success, exit retry loop for this batch

            except Exception as e:
                wait_time = 2 ** attempt
                logging.warning(f"Retrying fetch (attempt {attempt + 1}/{max_retries}) in {wait_time}s due to error: {e}")
                time.sleep(wait_time)
        else:
            # This block executes if the for loop completes without 'break' (all retries failed)
            logging.error(f"Failed to fetch PMIDs after {max_retries} retries for query '{query[:100]}...' at index {start_index}. Fetched {len(pmids)} so far.")
            # Break the main while loop to prevent getting stuck
            break
    return pmids


def get_fixed_month_interval(start_date):
    "for old dates of universal query"
    try:
        if start_date <= date.fromisoformat("1828-01-01"):
            return 564
        elif start_date <= date.fromisoformat("1843-01-01"):
            return 168
        elif start_date <= date.fromisoformat("1849-01-01"):
            return 72
        elif start_date <= date.fromisoformat("1854-01-01"):
            return 60
        elif start_date <= date.fromisoformat("1859-01-01"):
            return 59
        elif start_date <= date.fromisoformat("1865-01-01"):
            return 72
        else:
            return None
    except ValueError as e: # Handle potential invalid date strings if inputs are strings
        logging.error(f"Invalid date in get_fixed_month_interval: {e}")
        return None

def generate_date_batches(query, start_date, stop_date, target_papers_per_batch=10000, window_sizes=[365, 180, 90, 60, 40, 20, 10, 1], max_workers=2, full_text=False):
    """
    Generate date batches that are forced to stay within calendar years,
    but use adaptive windowing *within* each year to respect paper count limits.
    """

    print(f"Generating year-aligned batches with internal windowing from {start_date} to {stop_date}...\n")

    date_ranges = []
    total_papers = 0

    first_year = start_date.year
    last_year = stop_date.year

    def get_count_for_window(start, days):
        end = min(start + timedelta(days=days), stop_date)
        # Use consistent date formatting
        q = f'("{start.strftime("%Y/%m/%d")}"[PDat] : "{end.strftime("%Y/%m/%d")}"[PDat]) {query}'
        count = get_pubmed_count(q, full_text=full_text)
        return (start, end, count)

    for year in range(first_year, last_year + 1):
        year_start = date(year, 1, 1)
        year_end = date(year, 12, 31)

        batch_start = max(start_date, year_start)
        batch_end = min(stop_date, year_end)

        if batch_start > batch_end:
            continue

        current = batch_start

        while current <= batch_end:
            # Ensure candidates don't exceed batch_end
            candidates = [min(current + timedelta(days=w), batch_end) for w in window_sizes]

            results = []
            # Using ThreadPoolExecutor for parallel counting
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(get_count_for_window, current, (end - current).days)
                    for end in candidates
                ]
                # Collect results as they complete
                for f in as_completed(futures):
                    try:
                        results.append(f.result())
                    except Exception as e:
                        logging.error(f"Error while evaluating date window in year {year}: {e}")

            best = None
            # Sort results by window size (descending) to try larger windows first
            for start, end, count in sorted(results, key=lambda x: (x[1] - x[0]).days, reverse=True):

                if count <= target_papers_per_batch:
                    best = (start, end, count)
                    total_papers += count
                    break # Found an acceptable window, break inner loop

            if best:
                date_ranges.append(best)
                # Move to the day after the end of the selected window
                current = best[1] + timedelta(days=1)
            else:
                # Fallback if no window meets the target (e.g., even 1 day > target)
                fallback_end = min(current + timedelta(days=1), batch_end)
                # Recalculate count for the fallback window
                fallback_query = f'("{current.strftime("%Y/%m/%d")}"[PDat] : "{fallback_end.strftime("%Y/%m/%d")}"[PDat]) {query}'
                count = get_pubmed_count(fallback_query, full_text=full_text)
                if count > target_papers_per_batch:
                    logging.warning(f"1-day batch from {current} contains {count} papers, exceeding {target_papers_per_batch} threshold.")
                total_papers += count
                date_ranges.append((current, fallback_end, count))
                current = fallback_end + timedelta(days=1) # Move to next day


    print(f"\n Total batches generated: {len(date_ranges)}")
    print(f"Estimated total papers across all batches margin is 30%: {total_papers}")

    return date_ranges # Returns list of (start_date, end_date, count) tuples


def fetch_pmids_parallel(date_batches, query, full_text, max_workers):
    """Fetches PMIDs for a list of date batches in parallel."""
    results = []
    metadata = []

    def fetch_single_batch(batch_tuple):
        """Fetches PMIDs for a single batch tuple (start_date, end_date, count)."""
        # Unpack the tuple correctly (assuming 3 elements now)
        start_date, end_date, _ = batch_tuple # We don't use the pre-calculated count here
        # Format dates consistently for the Entrez query
        date_query = f'("{start_date.strftime("%Y/%m/%d")}"[PDat] : "{end_date.strftime("%Y/%m/%d")}"[PDat])'
        full_query = f"{date_query} {query}"
        pmids = get_list(full_query, full_text) # This returns a list of PMIDs
        # Create metadata using the actual count of fetched PMIDs
        meta = {"start": start_date, "end": end_date, "count": len(pmids)}
        return pmids, meta

    # Use ThreadPoolExecutor for parallel fetching
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks for each batch
        future_to_batch = {executor.submit(fetch_single_batch, batch): batch for batch in date_batches}
        # Process completed tasks
        for future in as_completed(future_to_batch):
            batch = future_to_batch[future]
            try:
                # Get the result (pmids list and metadata dict)
                pmids, meta = future.result()
                results.extend(pmids) # Add PMIDs to the main list
                metadata.append(meta) # Add metadata to the list
            except Exception as e:
                # Log errors for individual batches
                logging.error(f"Failed to fetch PMIDs for batch {batch}: {e}")

    return results, metadata # Return combined list of PMIDs and list of metadata dicts


def read_query_from_file(filename):
    """Read and clean query from a file."""
    try:
        with open(filename, 'r', encoding='utf-8') as file: # Specify encoding
            query = file.read().strip() # Strip whitespace/newlines
        if not query:
             logging.warning(f"Query file {filename} is empty.")
        return query
    except FileNotFoundError:
        logging.error(f"Query file not found: {filename}")
        return None
    except Exception as e:
        logging.error(f"Error reading query file {filename}: {e}")
        return None

def plot_density_over_time(batch_metadata):
    """
    Plots the number of papers published per year based on batch metadata.
    Expects batch_metadata to be a list of dictionaries with keys 'start', 'end', 'count'.
    'start' and 'end' should be datetime.date objects.
    """
    import matplotlib.pyplot as plt
    from collections import defaultdict

    if not batch_metadata:
        print("No batch metadata to plot.")
        return

    # Dictionary to accumulate paper counts per year
    yearly_counts = defaultdict(int)

    # Iterate over the list of dictionaries
    for item in batch_metadata:
        try:
            start_date = item["start"]  # Access the date object
            count = item["count"]
            # Ensure start_date is a date object before accessing .year
            if isinstance(start_date, date):
                year = start_date.year
                yearly_counts[year] += count
            else:
                 logging.warning(f"Invalid start date type in metadata item: {item}")
        except KeyError as e:
             logging.error(f"Missing key in batch metadata item: {e} in {item}")
        except Exception as e:
             logging.error(f"Error processing batch metadata item {item}: {e}")

    if not yearly_counts:
        print("No valid yearly counts to plot.")
        return

    # Sort years for plotting
    years = sorted(yearly_counts.keys())
    counts = [yearly_counts[year] for year in years]

    # Create the plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(years, counts, color='skyblue', edgecolor='black')
    plt.title("estimated number of Papers Published Per Year")
    plt.xlabel("Year")
    plt.ylabel("Number of Papers")
    plt.grid(axis='y', linestyle='--', alpha=0.7)


    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{int(height)}', ha='center', va='bottom', fontsize=8)

    plt.xticks(years, rotation=45)
    plt.tight_layout()
    plt.savefig(f"pmid_trend_{date.today().strftime('%Y%m%d')}.png", dpi=300, bbox_inches='tight')



def fetch_pmids_over_period(query_file, start=None, stop=None, full_text=False, plot=True, max_workers=2):
    """Main function to orchestrate fetching PMIDs over a period."""
    if not getattr(Entrez, 'email', None): # Check if Entrez.email is set
         Entrez.email = input("Please enter your email for NCBI Entrez: ").strip()
         if not Entrez.email:
             logging.error("Entrez.email is required. Exiting.")
             return []
    # Get API Key
    if not getattr(Entrez, 'api_key', None): # Check if Entrez.api_key is set
        api_input = input('Please enter NCBI API key *recommended* (or press Enter to skip): ').strip()
        API_KEY = api_input if api_input else None
        if API_KEY:
            Entrez.api_key = API_KEY # Set Entrez.api_key if using Biopython's Entrez
            logging.info("API key set for Entrez.")
        else:
            logging.warning("No API key provided. Rate limits may apply.")


    # Read query
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read or invalid query.")
        return []

    # Parse dates
    try:
        if start is None:
            start_date = date(1828, 1, 1)
        if stop is None:
            stop_date = date.today()
        else:
            stop_date = date.fromisoformat(stop)
    except ValueError as e:
        logging.error(f"Invalid date format: {e}")
        return []

    # Generate date batches
    print("Generating date batches...")
    date_batches = generate_date_batches(
        query=query,
        start_date=start_date,
        stop_date=stop_date,
        full_text=full_text,
        max_workers=max_workers,
    )

    # Plot if requested
    if plot and date_batches:
        # Prepare metadata for plotting (list of dicts)
        batch_metadata_preview = [
            {
                "start": start_date_obj,  # Keep as date object
                "end": end_date_obj,      # Keep as date object
                "count": count            # Use the precomputed count
            }
            for start_date_obj, end_date_obj, count in date_batches 
        ]
        print('Plotting the density of papers over time....')
        plot_density_over_time(batch_metadata_preview)
    elif plot:
        print("No batches generated, skipping plot.")

    # Fetch PMIDs in parallel
    print('Fetching the list of PMIDs....')
    if not date_batches:
        logging.warning("No date batches generated, returning empty list.")
        return []

    pmid_list, batch_metadata = fetch_pmids_parallel(
        date_batches=date_batches,
        query=query,
        full_text=full_text,
        max_workers=max_workers
    )

    # Remove duplicates from the final list
    unique_pmid_list = list(set(pmid_list))
    print(f'Total unique PMIDs fetched: {len(unique_pmid_list)}')

    # Print batch summary
    print('Batch metadata summary:')
    for meta in batch_metadata:
        print(f"  {meta['start']} to {meta['end']}: {meta['count']} papers")

    # Save PMIDs to file
    today_str = date.today().strftime("%d_%m_%Y")
    filename = f'pmid_list_{today_str}.txt'
    try:
        with open(filename, 'w') as file:
            for pmid in unique_pmid_list:
                file.write(f"{pmid}\n")
        print(f'PMID list saved to {filename}')
    except Exception as e:
        logging.error(f"Error saving PMIDs to file {filename}: {e}")

    return unique_pmid_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




################################################ ##################################################################################################
################################################ Fetch PubMed Metadata ############################################################# #################################################################################################################

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


def fetch_article(pmid: str):
    """Fetch a single article from PubMed by PMID."""
    article = fetcher.article_by_pmid(pmid)
    if article.pmid != pmid:
        logging.warning("Article with pmid=%r returned pmid=%r", pmid, article.pmid)
    return article

@retry_on_communication_error()
def fetch_articles(pmids: List[str], *, processes: Optional[int] = 1) -> Iterator:
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
        fetch_articles(pmids, processes=1),
        total=len(pmids),
        desc="Adding them",
        unit="row"
    ):
        articles_data.append({
            'PMID': article.pmid,
            'PMCID': article.pmc,
            'DOI': article.doi,
            'Title': article.title,
            'Authors': ', '.join(article.authors),
            'Year': article.year,
            'Journal': article.journal,
            'Volume': article.volume,
            'Issue': article.issue,
            'Pages': article.pages,
            'Abstract': article.abstract,
        })

    df = pd.DataFrame(articles_data)
    df['PMCID'] = df['PMCID'].apply(lambda x: f"PMC{x}" if pd.notnull(x) else x)
    return df







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pmid2pmcid(pmids, email):
    """Convert a list of PMIDs to PMCIDs using the NCBI ID Converter API."""

    base_url = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
    batch_size = 200
    results = []

    headers = {
        'User-Agent': f'{"pmid2pmcid"}/1.0 ({email})'
    }

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        params = {
            'tool': "pmid2pmcid",
            'email': email,
            'ids': ','.join(str(pmid) for pmid in batch),
            'format': 'json'
        }

        response = requests.get(base_url, params=params, headers=headers)
        if response.status_code != 200:
            print(f"Error fetching batch starting at index {i}: {response.status_code}")
            continue

        data = response.json()
        for record in data.get('records', []):
            pmid = record.get('pmid')
            pmcid = record.get('pmcid')
            results.append({'PMID': pmid, 'PMCID': pmcid})

        sleep(0.5)  # be polite to the API
    df = pd.DataFrame(results)    
    df = df.dropna(subset=['PMCID']).reset_index(drop=True)   # remove rows with None PMCIDs
    return df

def filter_oa_database(pmc_id_list):
    """
    Filters based on database list of PMCs that are available for full_text mining.

    Parameter:
    pmc_id_list (list): Filename of the list containing the PMC IDs.
    """
    # Read open access database file from NCBI server \size: ~230 mgb\
    oa_db_url = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.txt"

    pmc_id_list = pmc_id_list['PMCID'].tolist()

    # Read directly into pandas
    oa_file_list_df = pd.read_csv(oa_db_url, sep="\t", header=None, 
                    names=['File', 'Publication_info', 'PMCID', 'PMID', 'License'])

    # Filter oa_database based on PMC ID list
    filtered_oa_database = oa_file_list_df[oa_file_list_df["PMCID"].isin(pmc_id_list)]
    
    # Return specified columns
    return filtered_oa_database[['PMID', 'PMCID', 'Publication_info']]


def filter_df_oa_database(df):
    """
    Filters based on database list of PMCs that are available for full_text mining.
    Adds a column 'OA-noncomm' indicating if the paper is in the OA database.

    Parameter:
    df (DataFrame): DataFrame containing a 'PMCID' column
    
    Returns:
    DataFrame: Original dataframe with added 'OA-noncomm' column
    """
    # Read open access database file from NCBI server (~230 MB)
    oa_db_url = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.txt"
    
    try:
        # Download the OA database file
        response = requests.get(oa_db_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Read directly into pandas
        oa_file_list_df = pd.read_csv(oa_db_url, sep="\t", header=None, 
        names=['File', 'Publication_info', 'PMCID', 'PMID', 'License'])

        oa_pmc_ids = set(oa_file_list_df['PMCID'])

    
        
        # Add the 'OA-noncomm' column to the dataframe
        df['OA-noncomm'] = df['PMCID'].apply(
            lambda x: True if str(x).strip() in oa_pmc_ids else False
        )
        # Save the updated dataframe to a new CSV file
        df.to_csv("with_oa_column.csv", index=False)
        return df
    except Exception as e:
        print(f"Error processing OA database: {e}")
        return df



################################################################  ##############################################################################################
################################################################ ADD PUBLISHER from CRossref API ################################################################
# --- Configuration ---
# Number of DOIs per batch (keep small for URL length)
DOIS_PER_BATCH = 8
# Number of concurrent worker threads for DOI processing
MAX_WORKERS = 8 # Adjust based on testing/rate limits
# Number of concurrent worker threads for ISSN processing
MAX_ISSN_WORKERS = 8 # Adjust based on testing/rate limits, maybe lower than DOI workers
# Number of ISSNs per batch (adjust as needed, maybe keep small)
ISSN_BATCH_SIZE = 8# Adjust based on testing/rate limits

# Number of times to retry on a 500 error
MAX_RETRIES = 3
# Seconds to wait between retries
RETRY_DELAY = 0.2

# --- Rate Limiting Configuration ---
# Crossref API limit is 50 requests per second
CROSSREF_CALLS_PER_SECOND = 50
# Calculate the period based on the rate limit
CROSSREF_PERIOD_SECONDS = 1 # Per second

# --- Define the Rate-Limited Request Function ---
@sleep_and_retry # Automatically sleeps if the rate limit is exceeded before the call
@limits(calls=CROSSREF_CALLS_PER_SECOND, period=CROSSREF_PERIOD_SECONDS)
def rate_limited_get(url, **kwargs):
    """Wrapper for requests.get that enforces rate limiting."""
    return requests.get(url, **kwargs)



def query_crossref_doi_batch(doi_batch: list[str], mailto: str) -> dict[str, str | None]:
    """
    Queries Crossref for a batch of DOIs and returns a mapping of DOI to Publisher.
    Includes retry logic for HTTP 500 errors.
    Uses rate-limited requests.
    """
    WORKS_API_URL = 'https://api.crossref.org/works' # Fixed trailing spaces
    doi_publisher_map = {doi: None for doi in doi_batch if pd.notna(doi)}
    if not doi_batch:
        return doi_publisher_map

    valid_dois = [doi for doi in doi_batch if pd.notna(doi)]
    if not valid_dois:
        return doi_publisher_map

    doi_filter_str = ",".join([f"doi:{doi}" for doi in valid_dois])

    params = {
        'filter': doi_filter_str,
        'mailto': mailto,
        'select': 'DOI,publisher',
    }

    # retry
    for attempt in range(MAX_RETRIES + 1): # Try up to MAX_RETRIES times, plus the initial attempt
        try:
            # --- Use rate_limited_get instead of requests.get ---
            response = rate_limited_get(WORKS_API_URL, params=params, timeout=30)

            # Check for 500 error and retry if within retry limit
            if response.status_code == 500:
                if attempt < MAX_RETRIES:
                    print(f"    DOI Batch {valid_dois[:2]}...: Received 500 error (attempt {attempt + 1}/{MAX_RETRIES + 1}). Retrying in {RETRY_DELAY}s...")
                    time.sleep(RETRY_DELAY)
                    continue # Go to the next iteration of the retry loop
                else:
                    print(f"    DOI Batch {valid_dois[:2]}...: Received 500 error. Max retries ({MAX_RETRIES}) exceeded.")
                    # The map remains initialized with None for these DOIs
                    break # Exit the retry loop

            # If successful or other non-500 error, break out of retry loop
            response.raise_for_status() # Raise for other HTTP errors (4xx)
            data = response.json()

            items = data.get('message', {}).get('items', [])
            for item in items:
                doi = item.get('DOI')
                publisher = item.get('publisher')
                if doi and doi in doi_publisher_map:
                    doi_publisher_map[doi] = publisher
            break # Success, exit retry loop

        except requests.exceptions.RequestException as e:
            print(f"    DOI Batch {valid_dois[:2]}...: Request error (attempt {attempt + 1}): {e}")
            if attempt < MAX_RETRIES:
                 print(f"        Retrying in {RETRY_DELAY}s...")
                 time.sleep(RETRY_DELAY)
                 continue # Retry
            # If it's the last attempt or a non-retryable error, log and return current state (Nones)
            print(f"    DOI Batch {valid_dois[:2]}...: Failed after {MAX_RETRIES + 1} attempts.")
            break # Exit retry loop
        except Exception as e: # Catch other potential errors during JSON parsing etc.
            print(f"    DOI Batch {valid_dois[:2]}...: Unexpected error (attempt {attempt + 1}): {e}")
            if attempt < MAX_RETRIES:
                 print(f"        Retrying in {RETRY_DELAY}s...")
                 time.sleep(RETRY_DELAY)
                 continue
            print(f"    DOI Batch {valid_dois[:2]}...: Failed after {MAX_RETRIES + 1} attempts due to unexpected error.")
            break # Exit retry loop

    return doi_publisher_map

def query_crossref_issn_batch(issn_batch: list[str], mailto: str) -> dict[str, str | None]:
    """
    Queries Crossref /works endpoint for a BATCH of ISSNs and returns a mapping
    of ISSN to Publisher (from the first work found for that ISSN).
    """
    # Initialize map for all ISSNs in the batch
    issn_publisher_map = {issn: None for issn in issn_batch if pd.notna(issn) and str(issn).lower() not in ('nan', '')}

    if not issn_publisher_map:
        return issn_publisher_map

    # Implement retry logic for the entire batch
    for attempt in range(MAX_RETRIES + 1):
        # Track if any request in the batch failed with a retryable error (500)
        batch_has_retryable_error = False
        batch_error_details = []

        # Process each ISSN in the current batch attempt
        for issn in list(issn_publisher_map.keys()): # Iterate over a copy of keys
             if not issn or pd.isna(issn) or (isinstance(issn, str) and issn.lower() in ('nan', '')):
                 # Should ideally be filtered out already, but double-check
                 issn_publisher_map[issn] = None
                 continue

             # Construct the URL for the specific ISSN
             # Note: Using 'rows=1' to get only the first result and not load everything
             url_params = {
                 'filter': f'issn:{issn}',
                 'select': 'publisher',
                 'mailto': mailto,
                 'rows': '1'
             }
             # Build the query string properly
             from urllib.parse import urlencode
             query_string = urlencode(url_params)
             url = f"https://api.crossref.org/works?{query_string}"

             try:
                 # --- Using rate_limited_get  ---
                 response = rate_limited_get(url, timeout=30)

                 if response.status_code == 200:
                     data = response.json()
                     items = data.get('message', {}).get('items', [])
                     if items:
                         publisher = items[0].get('publisher')
                         issn_publisher_map[issn] = publisher
                     # If no items, publisher remains None (already set)
                 elif response.status_code == 404:
                     # ISSN not found, publisher remains None
                     pass
                 elif response.status_code == 500:
                     # Mark for retry, don't break inner loop yet to try others in batch
                     batch_has_retryable_error = True
                     batch_error_details.append(f"500 for ISSN {issn}")
                 else:
                     # Log other HTTP errors, but don't necessarily retry the whole batch
                     print(f"        ISSN {issn}: HTTP {response.status_code}")
                     # issn_publisher_map[issn] remains None

             except requests.exceptions.Timeout:
                 print(f"        ISSN {issn}: Timeout error (attempt {attempt + 1})")
                 batch_has_retryable_error = True # Consider timeouts retryable
                 batch_error_details.append(f"Timeout for ISSN {issn}")
             except requests.exceptions.RequestException as e:
                 print(f"        ISSN {issn}: Request error (attempt {attempt + 1}): {e}")
                 # Decide if other request errors should trigger a batch retry
                 # Let's retry for now, similar to 500s
                 batch_has_retryable_error = True
                 batch_error_details.append(f"Request error for ISSN {issn}: {e}")
             except ValueError as e: # JSON decode error
                 print(f"        ISSN {issn}: JSON decode error (attempt {attempt + 1}): {e}")
                 # Not necessarily retryable for the whole batch, but let's be cautious
                 batch_has_retryable_error = True
                 batch_error_details.append(f"JSON error for ISSN {issn}: {e}")
             except Exception as e:
                 print(f"        ISSN {issn}: Unexpected error (attempt {attempt + 1}): {e}")
                 # Not necessarily retryable for the whole batch, but let's be cautious
                 batch_has_retryable_error = True
                 batch_error_details.append(f"Unexpected error for ISSN {issn}: {e}")

        # After processing all ISSNs in the batch for this attempt
        if not batch_has_retryable_error:
            # Success or non-retryable errors for all, exit retry loop
            break
        else:
            # There were retryable errors
            if attempt < MAX_RETRIES:
                print(f"    ISSN Batch {issn_batch[:2]}...: Retryable errors occurred (attempt {attempt + 1}/{MAX_RETRIES + 1}). Details: {batch_error_details[:3]}... Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
                # Continue to the next attempt in the retry loop
            else:
                print(f"    ISSN Batch {issn_batch[:2]}...: Max retries ({MAX_RETRIES}) exceeded. Errors: {batch_error_details[:3]}...")
                # Exit retry loop, return map with Nones for failed ISSNs
                break # Exit retry loop

    return issn_publisher_map # Return the final map for this batch


def get_publishers_concurrent_from_issn(missing_df: pd.DataFrame, email: str, max_workers: int = MAX_ISSN_WORKERS) -> List[Optional[str]]:
    """
    Retrieves publishers for a DataFrame missing publisher info using concurrent ISSN lookups (batched).

    Args:
        missing_df (pd.DataFrame): DataFrame containing rows missing publisher info.
        email (str): Your email for polite API usage.
        max_workers (int): Maximum number of concurrent threads for ISSN lookups.

    Returns:
        List[Optional[str]]: A list of publisher names corresponding to the rows in missing_df.
    """
    print(f"Attempting to retrieve publishers for {len(missing_df)} records using concurrent ISSN lookup...")
    publishers = [None] * len(missing_df)

    # Prepare list of valid ISSNs
    issn_lookup_tasks = []

    for idx_in_missing_df, (original_index, row) in enumerate(missing_df.iterrows()):
        issn = row.get('ISSN')
        if pd.notna(issn) and not (isinstance(issn, str) and issn.lower() in ('nan', '')):
            if isinstance(issn, str):
                issn_clean = issn.split(',')[0].strip()
            else:
                issn_clean = str(issn).strip()

            if issn_clean and issn_clean.lower() not in ('nan', ''):
                issn_lookup_tasks.append((idx_in_missing_df, issn_clean))

    if not issn_lookup_tasks:
         print("No valid ISSNs found for lookup.")
         return publishers

    # Extract unique ISSNs for querying
    unique_issns = list(set(issn for _, issn in issn_lookup_tasks))
    print(f"  Found {len(issn_lookup_tasks)} ISSNs, {len(unique_issns)} unique.")

    # --- Create Batches of Unique ISSNs ---
    # This is the key change: group unique ISSNs into batches
    batches = [unique_issns[i:i + ISSN_BATCH_SIZE] for i in range(0, len(unique_issns), ISSN_BATCH_SIZE)]
    print(f"  Processing {len(unique_issns)} unique ISSNs in {len(batches)} batches (size {ISSN_BATCH_SIZE}) using {max_workers} workers...")
    # --- End of Batching ---

    final_issn_publisher_map = {}

    # Use tqdm for progress bar on unique ISSNs:
    pbar = tqdm(total=len(batches), desc="Processing ISSN batches", unit="Batches")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit batch jobs (each job processes a list of ISSNs)
        future_to_batch = {executor.submit(query_crossref_issn_batch, batch, email): batch for batch in batches}

        # Collect results as they complete
        for future in as_completed(future_to_batch):
            batch = future_to_batch[future]
            try:
                batch_result = future.result() # This is the dict {issn: publisher, ...} for the batch
                # Update the final map with results from this batch
                final_issn_publisher_map.update(batch_result)
                # Update progress bar - by number of batches or number of ISSNs in the batch?
                # Per ISSN (uncomment below, comment above)
                pbar.update(len(batch))
            except Exception as exc:
                print(f'\nISSN Batch {batch[:3]}... generated an exception: {exc}')
                # Update progress bar for failed batch
                pbar.update(len(batch)) # Or len(batch) if using Option 2 above
                # Initialize ISSNs in this failed batch with None in final map
                for issn in batch:
                    if pd.notna(issn):
                        final_issn_publisher_map[issn] = None

        pbar.close()

    # Map the results back to the original positions in the 'publishers' list
    # Create mapping from clean ISSN to list of positions (as before)
    issn_to_positions = {}
    for idx_in_result_list, clean_issn in issn_lookup_tasks:
        if clean_issn not in issn_to_positions:
            issn_to_positions[clean_issn] = []
        issn_to_positions[clean_issn].append(idx_in_result_list)

    # Assign publishers based on the final map
    for issn_key, position_list in issn_to_positions.items():
        publisher = final_issn_publisher_map.get(issn_key, None)
        for position in position_list:
            publishers[position] = publisher

    print("ISSN lookup processing complete.")
    return publishers

def get_publishers_concurrent(doi_list: List[str], mailto: str, max_workers: int = MAX_WORKERS) -> Dict[str, Optional[str]]:
    """
    Retrieves publishers for a list of DOIs using concurrent batch queries.

    Args:
        doi_list (List[str]): The list of DOIs.
        mailto (str): Your email for polite API usage.
        max_workers (int): Maximum number of concurrent threads.

    Returns:
        Dict[str, Optional[str]]: A dictionary mapping input DOIs to publisher names.
    """
    if not doi_list:
        return {}

    # For mapping, a set is fine
    unique_dois = list(set(doi for doi in doi_list if pd.notna(doi)))
    if not unique_dois:
        return {doi: None for doi in doi_list}

    # Create batches
    batches = [unique_dois[i:i + DOIS_PER_BATCH] for i in range(0, len(unique_dois), DOIS_PER_BATCH)]
    print(f"Processing {len(unique_dois)} unique DOIs in {len(batches)} batches using {max_workers} workers...")

    final_publisher_map = {}
    
    # tqdm
    # Wrap the batches iterator and set the total for the progress bar
    pbar = tqdm(total=len(unique_dois), desc="Processing DOIs", unit="DOI")

    # Use ThreadPoolExecutor for concurrency
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all batch jobs
        future_to_batch = {executor.submit(query_crossref_doi_batch, batch, mailto): batch for batch in batches}

        # Collect results as they complete
        for future in as_completed(future_to_batch):
            batch = future_to_batch[future]
            try:
                batch_result = future.result()
                # Update the final map with results from this batch
                final_publisher_map.update(batch_result)
                # Update the progress bar by the number of DOIs in the completed batch
                pbar.update(len(batch))
            except Exception as exc:
                print(f'\nBatch {batch[:3]}... generated an exception: {exc}')
                # Update progress bar even for failed batches
                pbar.update(len(batch))
                # Initialize DOIs in this failed batch with None
                for doi in batch:
                    if pd.notna(doi):
                        final_publisher_map[doi] = None

        # Close the progress bar
        pbar.close()

    # Ensure all original DOIs (including potential duplicates/NaNs) have an entry
    # This map will primarily be used via .map() on the DataFrame 'DOI' column
    full_map = {doi: final_publisher_map.get(doi, None) if pd.notna(doi) else None for doi in doi_list}
    print("DOI batch processing complete.")
    return full_map







# --- Main Processing Function ---
def process_publishers_concurrent(articles_df: pd.DataFrame, email: str) -> pd.DataFrame:
    """Process publisher information using concurrent Crossref API queries and return updated DataFrame."""
    print("Starting concurrent publisher processing...")
    articles_df = articles_df.copy()
    articles_df['publisher'] = None # Initialize publisher column

    # --- Step 1: Get publishers using DOIs via Concurrent Crossref batch API ---
    if 'DOI' not in articles_df.columns:
        print("Warning: 'DOI' column not found. Skipping Crossref DOI lookup.")
        missing_df = articles_df
    else:
        dois = articles_df['DOI'].tolist()

        # Get the publisher mapping using concurrency
        doi_to_publisher_full_map = get_publishers_concurrent(dois, email, MAX_WORKERS)

        # Map publishers back to the DataFrame using the original DOI list order
        # .map() is efficient for this
        articles_df['publisher'] = articles_df['DOI'].map(doi_to_publisher_full_map)

        # Identify rows where publisher is still missing after DOI lookup
        missing_df = articles_df[articles_df['publisher'].isnull() | (articles_df['publisher'] == '')]
        # Also consider rows where DOI was NaN
        missing_df = pd.concat([missing_df, articles_df[articles_df['DOI'].isnull()]], ignore_index=False).drop_duplicates()

    # --- Step 2: Use (ISSN) for missing publishers ---
    if not missing_df.empty:
        print(f"Found {len(missing_df)} records missing publisher info after DOI lookup. Attempting concurrent ISSN lookup...")
        # Use the NEW concurrent ISSN-based lookup function
        publisher_ids_from_issn = get_publishers_concurrent_from_issn(missing_df, email, MAX_ISSN_WORKERS)

        # Create a temporary DataFrame to hold the new publisher data for missing entries
        temp_missing_df = missing_df.reset_index(drop=False) # Keep original index as a column
        original_index_name = temp_missing_df.columns[0] if len(temp_missing_df.columns) > 0 else 'index'
        temp_missing_df['publisher_from_issn'] = publisher_ids_from_issn

        # Update the original DataFrame using the original index
        articles_df = articles_df.merge(
            temp_missing_df[['publisher_from_issn']], # Select only the new publisher column and the index
            left_index=True,
            right_index=True, # Merge on the original index (which was kept in the temp df)
            how='left'
        )
        # Fill the original 'publisher' column with values from 'publisher_from_issn' where 'publisher' is null
        articles_df['publisher'] = articles_df['publisher'].fillna(articles_df['publisher_from_issn'])
        # Drop the temporary column
        articles_df.drop(columns=['publisher_from_issn'], inplace=True)

        # --- Generate Filename with Date ---
        # Get today's date
        today_str = datetime.today().strftime('%Y-%m-%d') # Formats date as YYYY-MM-DD
        # Create the output filename
        output_filename = f"{today_str}_metawithPublisher.csv"


        # Save the result df to the dated CSV file
        print(f"Saving results to '{output_filename}'...")
        articles_df.to_csv(output_filename, index=False)


    print("Concurrent publisher processing complete.")
    return articles_df

############################################################################### ##########################################################################################################
############################################################################### Doenload PMC ##########################################################################################################

def download_pmc_articles(oa_pmcids, output_dir='./Full_text_jsons'):
    """
    Download full-text BioC JSON articles from PubMed Central.
    
    Args:
        oa_pmcids (list): List of PMCID strings to download
        output_dir (str): Output directory for JSON files (default: './Full_text_jsons')
    
    Returns:
        tuple: (success_count, failed_fulltext) where:
            - success_count: Number of successfully downloaded articles
            - failed_fulltext: List of PMCIDs that failed to download
    """
    # Create directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize lists to track failed downloads
    failed_fulltext = []

    # Initialize progress bar
    with tqdm(oa_pmcids, desc="Downloading articles", unit="article") as pbar:
        for pmcid in pbar:
            # Update description to show current PMCID
            pbar.set_postfix_str(f"PMCID: {pmcid}")
            
            # 1. Download full-text BioC JSON
            fulltext_url = f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmcid}/unicode'
            try:
                response = requests.get(fulltext_url)
                if response.status_code == 200:
                    with open(f'{output_dir}/{pmcid}.json', 'w', encoding='utf-8') as f:
                        f.write(response.text)
                else:
                    tqdm.write(f"Failed to download full-text for {pmcid} (Status: {response.status_code})")
                    failed_fulltext.append(pmcid)
            except Exception as e:
                tqdm.write(f"Error downloading {pmcid}: {str(e)}")
                failed_fulltext.append(pmcid)
    
    success_count = len(oa_pmcids) - len(failed_fulltext)
    
    # Print summary of failed downloads
    print("\nDownload Summary:")
    print(f"Successfully processed {success_count}/{len(oa_pmcids)} full-text files")
    if failed_fulltext:
        print("\nPMCIDs with full-text download failures:")
        print(failed_fulltext)
    
    return success_count, failed_fulltext



def format_author_name(author_str: str) -> str:
    """Convert 'surname:surname;given-names:firstname' to 'firstname surname'."""
    parts = author_str.split(';')
    surname = ""
    given_names = ""
    
    for part in parts:
        if part.startswith('surname:'):
            surname = part.replace('surname:', '').strip()
        elif part.startswith('given-names:'):
            given_names = part.replace('given-names:', '').strip()
    
    return f"{given_names} {surname}" if given_names and surname else author_str

def extract_text_from_json_to_dataframe1(directory: str, section_types: List[str]) -> pd.DataFrame:
    """
    Extracts content from JSON files with formatted author names and subtitle types.
    
    Args:
        directory: Path to directory containing JSON files
        section_types: List of section types to extract (e.g., ['ABSTRACT', 'INTRO'])
    
    Returns:
        DataFrame with columns: filename, section_type, text, subtitle
    """
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info(f"Extracting from JSON files in {directory}")

    data_rows = []
    processed_files = 0
    error_files = 0

    for filename in os.listdir(directory):
        if not filename.endswith(".json"):
            continue

        filepath = os.path.join(directory, filename)
        base_filename = os.path.splitext(filename)[0]
        
        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                try:
                    data = json.load(file)
                except json.JSONDecodeError as e:
                    logger.error(f"Invalid JSON in {filename}: {str(e)}")
                    error_files += 1
                    continue

                try:
                    documents = data[0].get('documents', [])
                    if not documents:
                        logger.warning(f"No documents in {filename}")
                        continue

                    passages = documents[0].get('passages', [])
                    if not passages:
                        logger.warning(f"No passages in {filename}")
                        continue

                    # Extract metadata from first passage
                    first_passage = passages[0]
                    first_infons = first_passage.get('infons', {})

                    
                    # Add PMID
                    if "article-id_pmid" in first_infons:
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "PMID",
                            "subtitle": None,
                            "text": first_infons['article-id_pmid']
                        })


                    # Add DOI
                    if 'article-id_doi' in first_infons:
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "DOI",
                            "subtitle": None,
                            "text": first_infons['article-id_doi']
                        })

                    # Add Year
                    if 'year' in first_infons:
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "YEAR",
                            "subtitle": None,
                            "text": first_infons['year']
                        })

                    # Add issue
                    if 'issue' in first_infons:
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "ISSUE",
                            "subtitle": None,
                            "text": first_infons['issue']
                        })

                    # Add Volume
                    if 'volume' in first_infons:
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "VOLUME",
                            "subtitle": None,
                            "text": first_infons['volume']
                        })

                    # Add formatted authors
                    i = 0
                    while f'name_{i}' in first_infons:
                        original_name = first_infons[f'name_{i}']
                        formatted_name = format_author_name(original_name)
                        data_rows.append({
                            "filename": base_filename,
                            "section_type": "AUTHOR",
                            "subtitle": None,
                            "text": formatted_name,
                        })
                        i += 1

                    # Add requested sections with subtitle
                    for passage in passages:
                        infons = passage.get('infons', {})
                        current_section = infons.get('section_type')
                        subtitle = infons.get('type')
                        
                        if current_section in section_types:
                            text = passage.get('text', '').strip()
                            if text:
                                data_rows.append({
                                    "filename": base_filename,
                                    "section_type": current_section,
                                    "subtitle": subtitle,
                                    "text": text
                                })

                    processed_files += 1

                except Exception as e:
                    logger.error(f"Error processing {filename}: {str(e)}")
                    error_files += 1

        except IOError as e:
            logger.error(f"Error reading {filename}: {str(e)}")
            error_files += 1

    # Create DataFrame
    df = pd.DataFrame(data_rows)
    
    # Log summary
    logger.info(f"Processed {processed_files} files, {error_files} errors")
    logger.info(f"Extracted {len(df)} entries")    
    return df






def add_metadata_to_dataframe(df):
    """
    Adds PubMed metadata to a dataframe.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe containing columns: filename, section_type, subtitle, text
    Returns:
    --------
    pd.DataFrame
        New dataframe with metadata inserted before each document's first occurrence
    """
    
    # 1. Create filename metadata map
    file_meta = {}
    for idx, row in df.iterrows():
        filename = row['filename']
        if filename not in file_meta:
            file_meta[filename] = {
                'pmid': row['text'] if row['section_type'] == 'PMID' else None,
                'first_idx': idx
            }
    
    # 2. Fetch metadata for each document
    metadata = {}
    for filename, info in tqdm(file_meta.items(), desc='Fetching metadata'):
        if not info['pmid']:
            continue
            
        try:
            meta_df = fetch_articles_meta([info['pmid']])
            if not meta_df.empty:
                metadata[filename] = [
                    {
                        'filename': filename,
                        'section_type': field,
                        'subtitle': None,
                        'text': str(value)
                    }
                    for _, record in meta_df.iterrows()
                    for field, value in record.items()
                ]
        except Exception as e:
            print(f'Metadata error for {filename}: {str(e)[:100]}...')
    
    # 3. Reconstruct dataframe with inserted metadata
    new_rows = []
    processed_files = set()
    
    for idx, row in df.iterrows():
        filename = row['filename']
        
        # Insert metadata before first occurrence
        if filename not in processed_files:
            processed_files.add(filename)
            if filename in metadata:
                new_rows.extend(metadata[filename])
        
        # Add original row
        new_rows.append(row.to_dict())
    
    # 4. Return new dataframe
    result_df = pd.DataFrame(new_rows)
    
    return result_df

##### Download Sypplemantary XML files from BioC API #####

def download_supplementary_materials(oa_pmcids, output_dir="supplementary_materials", format="bioc_xml", delay=1.0):
    """
    Downloads all supplementary materials for a list of PMCIDs using the NCBI BioC Supplementary Materials API.

    Parameters:
        oa_pmcids (list): List of PMCIDs (e.g., ['PMC1234567', 'PMC2345678']).
        output_dir (str): Directory to save the downloaded supplementary materials.
        format (str): Format of the supplementary materials ('bioc_xml' or 'bioc_json').
        delay (float): Delay in seconds between requests to respect NCBI's rate limits.

    Returns:
        list: List of PMCIDs that had no supplementary materials or failed to download.
    """
    base_url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi"
    os.makedirs(output_dir, exist_ok=True)

    unsaved_pmcs = []

    for pmcid in oa_pmcids:
        url = f"{base_url}/{format}/{pmcid}/all"
        print(f"Downloading supplementary materials for {pmcid} from {url}")

        try:
            response = requests.get(url)
            response.raise_for_status()
        except requests.RequestException as e:
            print(f"Failed to download supplementary materials for {pmcid}: {e}")
            unsaved_pmcs.append(pmcid)
            continue

        content_text = response.text.strip()
        if "No result can be found" in content_text:
            print(f"No supplementary materials found for {pmcid}. Skipping save.")
            unsaved_pmcs.append(pmcid)
        else:
            file_extension = "xml" if format == "bioc_xml" else "json"
            file_path = os.path.join(output_dir, f"{pmcid}.{file_extension}")
            with open(file_path, "wb") as file:
                file.write(response.content)
            print(f"Saved supplementary materials for {pmcid} to {file_path}")

        sleep(delay)

    print(f"\nSummary: {len(oa_pmcids)} total PMCIDs processed.")
    print(f"{len(unsaved_pmcs)} had no supplementary materials or failed to download.")
    print(f"{len(oa_pmcids) - len(unsaved_pmcs)} successfully saved.")

    return unsaved_pmcs


#############  full_text creat + add xml of supplementary but u need to download supplementary yourself  #############

def extract_text_from_bioc_xml(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    records = []
    filename = os.path.basename(xml_path)
    for doc in root.findall(".//document"):
        document_id = doc.find("id").text if doc.find("id") is not None else None
        passages = doc.findall("passage")

        for passage in passages:
            passage_text = passage.findtext("text", default="")
            infons = {infon.attrib["key"]: infon.text for infon in passage.findall("infon")}
            record = {
                "filename": filename,
                "source": infons.get("source", None),
                "document_id": document_id,
                "type": infons.get("type", None),
                "text": passage_text.strip()
            }
            records.append(record)

    return pd.DataFrame(records)

def extract_text_from_json_to_dataframe(directory: str, section_types: List[str], xml_directory: str = None) -> pd.DataFrame:
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info(f"Extracting from JSON files in {directory}")

    data_rows = []
    processed_files = 0
    error_files = 0

    for filename in os.listdir(directory):
        if not filename.endswith(".json"):
            continue

        filepath = os.path.join(directory, filename)
        base_filename = os.path.splitext(filename)[0]

        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                try:
                    data = json.load(file)
                except json.JSONDecodeError as e:
                    logger.error(f"Invalid JSON in {filename}: {str(e)}")
                    error_files += 1
                    continue

                try:
                    documents = data[0].get('documents', [])
                    if not documents:
                        logger.warning(f"No documents in {filename}")
                        continue

                    passages = documents[0].get('passages', [])
                    if not passages:
                        logger.warning(f"No passages in {filename}")
                        continue

                    first_passage = passages[0]
                    first_infons = first_passage.get('infons', {})

                    pmid = first_infons.get("article-id_pmid", "NOPMID")
                    entry_counter = 0
                    order_counter = 1
                    
                    # Sections
                    for passage in passages:
                        infons = passage.get('infons', {})
                        current_section = infons.get('section_type')
                        subtitle = infons.get('type')
                        if current_section in section_types:
                            text = passage.get('text', '').strip()
                            if text:
                                data_rows.append({
                                    "EntryID": f"{pmid}_{entry_counter:03}",
                                    "section_type": current_section,
                                    "subtitle": subtitle,
                                    "text": text,
                                    "Order": order_counter
                                })
                                entry_counter += 1
                                order_counter += 1

                    # Supplementary from XML (if available)
                    if xml_directory:
                        xml_path = os.path.join(xml_directory, base_filename + ".xml")
                        if os.path.isfile(xml_path):
                            xml_df = extract_text_from_bioc_xml(xml_path)
                            for _, row in xml_df.iterrows():
                                text = row['text']
                                if text:
                                    data_rows.append({
                                        "EntryID": f"{pmid}_{entry_counter:03}",
                                        "section_type": "SUPPLEMENT",
                                        "subtitle": row.get('type', None),
                                        "text": text,
                                        "Order": order_counter
                                    })
                                    entry_counter += 1
                                    order_counter += 1

                    processed_files += 1

                except Exception as e:
                    logger.error(f"Error processing {filename}: {str(e)}")
                    error_files += 1

        except IOError as e:
            logger.error(f"Error reading {filename}: {str(e)}")
            error_files += 1

    df = pd.DataFrame(data_rows)
    logger.info(f"Processed {processed_files} files, {error_files} errors")
    logger.info(f"Extracted {len(df)} entries")
    return df

###################################################################### BIOPYTHON:ENTREZ ######################################################################
def fetch_chunk_metadata(pmid_chunk: List[str]) -> pd.DataFrame:
    """Fetch and parse metadata for a single chunk of PMIDs."""
    if not pmid_chunk:
        logging.debug("Received empty PMID chunk.")
        return pd.DataFrame()

    records = []
    # Define batch size for efetch within the chunk (must be <= 10,000, 500 is safe)
    efetch_batch_size = 9999

    try:
        # 1. EPost this specific chunk
        # Ensure all PMIDs are strings and stripped
        pmid_string = ",".join(str(pmid).strip() for pmid in pmid_chunk if pmid)
        if not pmid_string:
            logging.warning("PMID string is empty after processing chunk. Skipping.")
            return pd.DataFrame() # Return empty DataFrame for this chunk

        epost_handle = Entrez.epost(db="pubmed", id=pmid_string)
        epost_result = Entrez.read(epost_handle)
        epost_handle.close()
        
        webenv = epost_result.get("WebEnv")
        query_key = epost_result.get("QueryKey")
        
        if not webenv or not query_key:
            logging.error(f"Failed to get WebEnv or QueryKey from epost for chunk. Result: {epost_result}")
            return pd.DataFrame() # Return empty DataFrame for this chunk
            
        total_in_chunk = len(pmid_chunk) # Number of PMIDs in this chunk
        logging.debug(f"epost successful for chunk (WebEnv: ...{webenv[-10:]}..., QueryKey: {query_key}). Total PMIDs: {total_in_chunk}")

        # 2. EFetch records for this chunk in sub-batches
        for start in range(0, total_in_chunk, efetch_batch_size):
            current_batch_size = min(efetch_batch_size, total_in_chunk - start)
            stream = None # Initialize stream for this sub-batch
            logging.debug(f"Fetching sub-batch: records {start} to {start + current_batch_size - 1} (size: {current_batch_size}) for chunk...")
            
            try:
                # Fetch the sub-batch using the WebEnv session for this chunk
                stream = Entrez.efetch(
                    db="pubmed",
                    rettype="medline",        # Correct for Medline parser
                    retmode="text",           # Correct for Medline parser
                    retstart=start,           # Start index within this chunk's WebEnv session
                    retmax=current_batch_size, # Number of records to fetch in this call
                    webenv=webenv,            # WebEnv from epost for this chunk
                    query_key=query_key,      # QueryKey from epost for this chunk
                )

                # Parse the sub-batch
                batch_iterator = Medline.parse(stream)

                # --- Iterate through articles in the sub-batch and extract data ---
                for article in batch_iterator:
                    # --- Extract and clean fields inline ---

                    # PMID, PMCID, Title, Abstract, Authors
                    pmid = article.get("PMID", None)
                    pmcid = article.get("PMC", None)
                    title = article.get("TI", None)
                    abstract = article.get("AB", None)
                    authors = "; ".join(article.get("AU", []))

                    # --- DOI from AID (inline) ---
                    doi = ""
                    aid_list = article.get("AID", [])
                    for aid in aid_list:
                        if "[doi]" in aid.lower(): # Case-insensitive check
                            doi = aid.split(" ")[0]  # Extract DOI part before ' [doi]'
                            break

                    # Extract dates from PHST field
                    phst_entries = article.get("PHST", [])
                    pubmed_date = None
                    pmc_release_date = None

                    for entry in phst_entries:
                        if "[pubmed]" in entry:
                            # Extract date part and remove time
                            date_part = entry.split()[0]  # Gets YYYY/MM/DD part
                            pubmed_date = date_part.replace("/", "-")  # Convert to YYYY-MM-DD
                        elif "[pmc-release]" in entry:
                            # Extract date part and remove time
                            date_part = entry.split()[0]  # Gets YYYY/MM/DD part
                            pmc_release_date = date_part.replace("/", "-")  # Convert to YYYY-MM-DD

                    date_publication = article.get("DP", None) 

                    # --- ISSN: Extract first xxxx-xxxx (inline) ---
                    raw_issn = article.get("IS", "")
                    issn_match = re.search(r"\d{4}-\d{4}", raw_issn)
                    issn = issn_match.group(0) if issn_match else None

                    # Journal, Issue, Volume
                    issue = article.get("IP", None)
                    volume = article.get("VI", None)
                    journal_abbrev = article.get("TA", None)
                    journal_full = article.get("JT", None)

                    # --- Separate Fields for MH, OT, NM, RN ---
                    mesh_terms = article.get("MH", None)
                    other_terms = article.get("OT", None)
                    substance_names = article.get("NM", None)

                    # Clean RN: extract text inside parentheses
                    registry_numbers_cleaned = []
                    for rn_entry in article.get("RN", []):
                        match = re.search(r"\(([^)]+)\)", rn_entry)
                        if match:
                            registry_numbers_cleaned.append(match.group(1))
                        else:
                            registry_numbers_cleaned.append(rn_entry)

                    # --- Build final record ---
                    record = {
                        "PMID": pmid,
                        "PMCID": pmcid,
                        "Title": title,
                        "Abstract": abstract,
                        "Journal Title": journal_full,
                        "Journal Title Abbreviation": journal_abbrev,
                        "Authors": authors,
                        "DOI": doi,
                        "Pubmed date": pubmed_date,
                        "PMC Release Date": pmc_release_date,
                        "Date of Publication": date_publication,
                        "ISSN": issn,
                        "Issue": issue,
                        "Volume": volume,
                        # Individual term lists
                        "MeSH Terms (MH)": mesh_terms,
                        "Other Terms (OT)": other_terms,
                        "Substance Names (NM)": substance_names,
                        "Registry Numbers (RN)": registry_numbers_cleaned 
                    }

                    records.append(record)
                # --- End of article iteration for this sub-batch ---

            except Exception as e: # Catch errors from efetch or parsing this sub-batch
                logging.error(f"Error fetching or parsing sub-batch starting at index {start} for chunk (PMIDs ~{pmid_chunk[start:start+3]}...): {e}")
            finally:
                # Ensure the stream for this sub-batch is closed
                if stream is not None:
                    try:
                        stream.close()
                        logging.debug(f"Stream for sub-batch {start}-{start + current_batch_size - 1} closed.")
                    except Exception as close_error:
                        logging.warning(f"Error closing stream for sub-batch {start}-{start + current_batch_size - 1}: {close_error}")


    except Exception as e: # Catch errors from the overall chunk processing (e.g., epost failure)
        logging.error(f"Failed to process chunk starting with PMID {pmid_chunk[0] if pmid_chunk else 'N/A'}: {e}")
        # Return any records collected so far for this chunk, or an empty DataFrame
        return pd.DataFrame(records)

    return pd.DataFrame(records)


# --- main function ---
def fetch_parse_pubmed_metadata(pmid_list: List[str]) -> pd.DataFrame:
    """Fetches and parses PubMed metadata for a large list of PMIDs using concurrency."""
    # Ensure Entrez is configured
    if not getattr(Entrez, 'email', None):
        Entrez.email = input("Enter your email (required by NCBI): ").strip()
    if not Entrez.email:
        logging.error("Entrez.email is required.")
        return pd.DataFrame()
    if not getattr(Entrez, 'api_key', None):
        Entrez.api_key = input("Enter your NCBI API key (optional but recommended): ").strip() or None
        logging.info("API key set.")

    # 1. Clean and deduplicate the main list
    cleaned_pmids = [str(pmid).strip() for pmid in pmid_list if pmid]
    unique_pmid_list = list(set(pmid for pmid in cleaned_pmids if pmid))
    if not unique_pmid_list:
        logging.warning("No valid PMIDs provided.")
        return pd.DataFrame()

    # 2. Divide into chunks (e.g., 9000 PMIDs per chunk)
    chunk_size_for_epost = 9999
    pmid_chunks = list(chunk_list(unique_pmid_list, chunk_size_for_epost))
    logging.info(f"Divided {len(unique_pmid_list)} PMIDs into {len(pmid_chunks)} chunks of max {chunk_size_for_epost}.")

    all_results = []
    max_workers = 8 # Adjust based on your API key limit


    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all chunks
        future_to_chunk_index = {
            executor.submit(fetch_chunk_metadata, chunk): i
            for i, chunk in enumerate(pmid_chunks)
            }

        total_chunks = len(pmid_chunks)

        with tqdm(as_completed(future_to_chunk_index), 
                total=total_chunks, 
                desc="Processing PMID Chunks", 
                unit="chunk") as progress_bar:
            
            # Collect results as they complete
            for future in progress_bar: # Iterate over the tqdm-wrapped iterator
                chunk_index = future_to_chunk_index[future]
                try:
                    chunk_df = future.result()
                    all_results.append(chunk_df)
                    # Optional: Update the progress bar description with live stats
                    collected_records = sum(len(df) for df in all_results)
                    progress_bar.set_postfix({"Records": f"{collected_records}/{len(unique_pmid_list)}"
            })
            
                    logging.debug(f"Chunk {chunk_index} processed successfully. Records: {len(chunk_df)}")
                except Exception as e:
                    logging.error(f"Chunk {chunk_index} generated an exception: {e}")
                    # Still update the counter even on error
                    collected_records = sum(len(df) for df in all_results)
                    progress_bar.set_postfix({
                    "Records": f"{collected_records}/{len(unique_pmid_list)}"
                    })


    # 4. Concatenate all results
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        logging.info(f"All chunks processed. Total records collected: {len(final_df)}")
        return final_df
    else:
        logging.error("No results collected from any chunk.")
        return pd.DataFrame()

# --- Helper function ---
def chunk_list(lst, chunk_size):
    """Yield successive chunk_size chunks from lst."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]
#############################################################################################################################
def save_pubmed_xml(pmids: List[str], output_dir: str):
    '''Save PubMed records as individual XML files using Biopython Entrez.'''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for pmid in tqdm(pmids, desc="Saving XML files"):
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        xml_data = handle.read()
        handle.close()

        if isinstance(xml_data, bytes):
            xml_data = xml_data.decode("utf-8")

        file_path = os.path.join(output_dir, f"{pmid}.xml")
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(xml_data)



####################################################            ################################################################
#################################################### Intact and Corrum data processing  ####################################################
def add_intact_corrum(df, pos_file="intact.txt", neg_file="intact_negative.txt"):
    """
    Process IntAct and Corrum database files to identify positive and negative papers,
    and classify papers in metadata.

    IntAct database zipped is ~1 gb "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip" 

    
    Parameters:
    pos_file (str): Path to positive interactions file
    neg_file (str): Path to negative interactions file
    
    Returns:
    dict: Dictionary containing processed data and statistics
    """
    
    # Load data
    pos = pd.read_csv(pos_file, sep="\t")
    neg = pd.read_csv(neg_file, sep="\t")

    # Extract PubMed IDs
    p = set([entry.split("pubmed:")[1].split('|')[0].strip() for entry in pos["Publication Identifier(s)"]])
    n = set([entry.split("pubmed:")[1].split('|')[0].strip() for entry in neg["Publication Identifier(s)"]])

    # Remove None values
    p.discard(None)
    n.discard(None)

    # Print statistics
    shared_count = len([i for i in n if i in p])
    print(f"The Intact DB has {len(p)} positive papers and {len(n)} negative papers.")
    print(f"Of the negative papers, {shared_count} are shared in the positives")
    # Convert sets to lists for isin()
    all_intact_pmids = list(n.union(p))
    df["intact"] = df["PMID"].isin(all_intact_pmids).astype(int)
    # Classify papers in metadata
    int_positives = df[df['intact'] == 1]
    int_unknowns = df[df['intact'] == 0]
    print(f"Intact Positive papers: {len(int_positives)}")
    print(f"Intact Unknown papers: {len(int_unknowns)}")

    # Load Corum data
    c=[str(i) for i in set(cor["pmid"])]
    df["corrum"] = df["PMID"].isin(c).astype(int)
    # Classify papers in metadata
    cor_positives = df[df['corrum'] == 1]
    cor_unknowns = df[df['corrum'] == 0]
    print(f"Corrum Positive papers: {len(cor_positives)}")
    print(f"Corrum Unknown papers: {len(cor_unknowns)}")
    
    # Return results
    return df


def add_complex_portal(input_dir: str, data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract PMIDs from all .tsv files in a directory and compare with metadata DataFrame.
    as there is a current directory with tsv files downloaded from complex portal "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/"
    Adds a 'Complex_portal' column to the DataFrame indicating presence (1) or

    Args:
        input_dir (str): Path to folder with .tsv files.
        data (pd.DataFrame): Metadata DataFrame with 'PMID' column.

    Returns:
        pd.DataFrame: Updated DataFrame with 'Complex_portal' column (1/0).
    """
    # --- Extract PMIDs ---
    complex_pmids = set()
    pattern = re.compile(r"pubmed:(\d+)")

    for filename in os.listdir(input_dir):
        if filename.endswith(".tsv"):
            filepath = os.path.join(input_dir, filename)
            df = pd.read_csv(filepath, sep="\t")

            # Search all columns for PMIDs
            for col in df.columns:
                for entry in df[col].dropna().astype(str):
                    matches = pattern.findall(entry.lower())
                    complex_pmids.update(matches)

    print(f"Papers in the Complex Portal: {len(complex_pmids)}")

    # --- Add Boolean Column to Data ---
    data["Complex_portal"] = data["PMID"].astype(str).isin(complex_pmids).astype(int)

    # --- Classify Papers ---
    complex_positives = data[data["Complex_portal"] == 1]
    complex_unknowns = data[data["Complex_portal"] == 0]

    positives = len(complex_positives)
    total_complex = len(complex_pmids)
    percentage = (positives / total_complex) * 100 if total_complex > 0 else 0

    print(f"Complex_portal Positive papers: {positives}")
    print(f"Complex_portal Unknown papers: {len(complex_unknowns)}")
    print(f"Coverage: {positives}/{total_complex} ({percentage:.2f}%)")
    data.to_csv("metadata_with_complex_portal.csv", index=False)

    return data



############################################################################################################
############################################# flag predatory pubs and jours ################################

def get_google_sheet_csv_url(google_url: str) -> str:
    """
    Convert Google Sheets URL to direct CSV download URL.
    """
    # Handle different Google Sheets URL formats
    if '/d/' in google_url:
        # Extract document ID
        doc_id = google_url.split('/d/')[1].split('/')[0]
        # Construct CSV export URL
        csv_url = f"https://docs.google.com/spreadsheets/d/{doc_id}/export?format=csv"
        return csv_url
    else:
        # Try to extract ID from URL parameters
        parsed = urlparse(google_url)
        query_params = parse_qs(parsed.query)
        if 'id' in query_params:
            doc_id = query_params['id'][0]
            csv_url = f"https://docs.google.com/spreadsheets/d/{doc_id}/export?format=csv"
            return csv_url
    
    raise ValueError(f"Could not parse Google Sheets URL: {google_url}")

def download_google_sheet_as_df(google_url: str, local_path: str = None) -> pd.DataFrame:
    """
    Download Google Sheet as pandas DataFrame and optionally save locally.
    """
    try:
        csv_url = get_google_sheet_csv_url(google_url)
        # Download without headers since the files don't have them
        df = pd.read_csv(csv_url, header=None)
        
        if local_path:
            df.to_csv(local_path, index=False, header=False)
            print(f"Downloaded and saved to: {local_path}")
            
        return df
    except Exception as e:
        print(f"Error downloading Google Sheet from {google_url}: {e}")
        return None

def flag_predatory_publications(
    data: pd.DataFrame,
    journal_col: str = "Journal Title",
    publisher_col: str = "publisher",
    bealls_journals_path: str = "./Predatory_publishers/bealls_journals.txt",
    bealls_publishers_path: str = "./Predatory_publishers/bealls_publishers.txt",
    list_download: str = 'y',
    output_col: str = "predatory",
    drop_clean_cols: bool = True
) -> pd.DataFrame:
    """
    Flags rows in a DataFrame as predatory (1) if the journal or publisher appears in 
    Beall's List or The Predatory Journals/Publishers List 2025.
    """
    
    # --- Clean function (removes parentheses (xxxx)) ---
    def clean_name(name: str) -> str:
        if pd.isna(name):
            return ""
        return re.sub(r"\s*$$.*?$$\s*", "", str(name)).strip()

    journal_set = set()
    publisher_set = set()

    # Google Sheets URLs
    JOURNALS_SHEET_URL = "https://drive.google.com/open?id=1Qa1lAlSbl7iiKddYINNsDB4wxI7uUA4IVseeLnCc5U4"
    PUBLISHERS_SHEET_URL = "https://drive.google.com/open?id=1BHM4aJljhbOAzSpkX1kXDUEvy6vxREZu5WJaDH6M1Vk"

    # --- Load Beall's Journals ---
    try:
        with open(bealls_journals_path, "r", encoding="utf-8") as f:
            journal_list = [clean_name(line.strip()) for line in f if line.strip()]
        journal_set.update(journal_list)
        print(f"Total unique predatory journals loaded from Beall's list: {len(journal_set)}")
    except FileNotFoundError:
        print(f"Warning: Beall's journals file not found at {bealls_journals_path}. Skipping.")

    # --- Load The Predatory Journals List ---
    if list_download.lower() != 'n':
        print("Downloading predatory journals list from Google Sheets...")
        try:
            # Download journals list (no headers)
            preda_jous = download_google_sheet_as_df(JOURNALS_SHEET_URL)
            if preda_jous is not None and len(preda_jous.columns) >= 2:
                # Second column contains journal names
                journal_names = preda_jous.iloc[:, 1].apply(clean_name).tolist()
                journal_set.update(set(journal_names))
                print(f"Total unique predatory journals loaded: {len(journal_set)}")
            else:
                print("Warning: Could not load predatory journals from Google Sheets.")
        except Exception as e:
            print(f"Warning: Error loading predatory journals from Google Sheets: {e}")
    else:
        print("Skipping Google Sheets journals download.")

    # --- Load Beall's Publishers ---
    try:
        with open(bealls_publishers_path, "r", encoding="utf-8") as f:
            publisher_list = [clean_name(line.strip()) for line in f if line.strip()]
        publisher_set.update(publisher_list)
        print(f"Total unique predatory publishers loaded from Beall's list: {len(publisher_set)}")
    except FileNotFoundError:
        print(f"Warning: Beall's publishers file not found at {bealls_publishers_path}. Skipping.")

    # --- Load The Predatory Publishers List ---
    if list_download.lower() != 'n':
        print("Downloading predatory publishers list from Google Sheets...")
        try:
            # Download publishers list (no headers)
            preda_pub = download_google_sheet_as_df(PUBLISHERS_SHEET_URL)
            if preda_pub is not None and len(preda_pub.columns) >= 2:
                # Second column contains publisher names
                publisher_names = preda_pub.iloc[:, 1].apply(clean_name).tolist()
                publisher_set.update(set(publisher_names))
                print(f"Total unique predatory publishers loaded: {len(publisher_set)}")
            else:
                print("Warning: Could not load predatory publishers from Google Sheets.")
        except Exception as e:
            print(f"Warning: Error loading predatory publishers from Google Sheets: {e}")
    else:
        print("Skipping Google Sheets publishers download.")

    # --- Check if specified columns exist ---
    if journal_col not in data.columns:
        print(f"Warning: Journal column '{journal_col}' not found in DataFrame. Available columns: {list(data.columns)}")
    if publisher_col not in data.columns:
        print(f"Warning: Publisher column '{publisher_col}' not found in DataFrame. Available columns: {list(data.columns)}")

    # --- Clean both columns ---
    journal_clean_col = f"{journal_col}_clean"
    publisher_clean_col = f"{publisher_col}_clean"

    # Only process columns that exist
    if journal_col in data.columns:
        data[journal_clean_col] = data[journal_col].apply(clean_name)
    if publisher_col in data.columns:
        data[publisher_clean_col] = data[publisher_col].apply(clean_name)

    # --- Flag if journal OR publisher is in their lists ---
    journal_flag = pd.Series([False] * len(data), index=data.index)
    publisher_flag = pd.Series([False] * len(data), index=data.index)
    
    if journal_col in data.columns:
        journal_flag = data[journal_clean_col].isin(journal_set)
    if publisher_col in data.columns:
        publisher_flag = data[publisher_clean_col].isin(publisher_set)
    
    data[output_col] = (journal_flag | publisher_flag).astype(int)

    if drop_clean_cols:
        cols_to_drop = []
        if journal_col in data.columns:
            cols_to_drop.append(journal_clean_col)
        if publisher_col in data.columns:
            cols_to_drop.append(publisher_clean_col)
        if cols_to_drop:
            data.drop(cols_to_drop, axis=1, inplace=True)

    return data

#########################################################################################################################
#########################################################################################################################