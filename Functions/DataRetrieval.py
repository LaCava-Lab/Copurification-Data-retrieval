import logging
import time
import numpy as np
import pandas as pd
import requests
import json
import csv
import os as os
import urllib.parse
from tqdm import tqdm  
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta
from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type, stop_after_attempt, wait_exponential
from metapub import PubMedFetcher, PubMedArticle, pubmedcentral
from multiprocessing.pool import ThreadPool
from typing import List, Iterator, Optional
from habanero import Crossref
from concurrent.futures import ThreadPoolExecutor, as_completed
from Reference_files.keys import API_KEY as API_KEY

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
        desc="Adding them",
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


import os
import requests
from tqdm import tqdm

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


from time import sleep


def get_pubmed_count(query, full_text=False, api_key=None):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    db = "pmc" if full_text else "pubmed"
    data = {
        "db": db,
        "term": query,
        "retmax": 0,
        "retmode": "json"
    }
    if api_key:
        data["api_key"] = api_key

    rate_limit = 10 if api_key else 3
    delay = 1 / rate_limit
    time.sleep(delay)  # rate limiting

    max_retries = 5
    for attempt in range(max_retries):
        try:
            response = requests.post(url, data=data)
            response.raise_for_status()
            return int(response.json()['esearchresult'].get('count', 0))
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                wait_time = 2 ** attempt
                logging.warning(f"Rate limit hit (HTTP 429). Retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            else:
                raise
        except Exception as e:
            logging.error(f"Unexpected error in get_pubmed_count: {e}")
            return 0
    logging.error(f"Failed to fetch PubMed count for query: {query[:200]} after {max_retries} retries.")
    return 0


def get_list(query, full_text, api_key=None):
    if api_key is None:
        api_key = API_KEY
    num_of_articles = 500
    start_index = 0
    pmids = []
    max_retries = 5

    while True:
        for attempt in range(max_retries):
            try:
                pmid_batch = fetcher.pmids_for_query(
                    query,
                    retstart=start_index,
                    retmax=num_of_articles,
                    pmc_only=full_text, api_key=API_KEY
                )
                pmids.extend(pmid_batch)
                start_index = len(pmids)
                if len(pmid_batch) < num_of_articles:
                    return pmids
                break  # Exit retry loop if successful
            except Exception as e:
                wait_time = 2 ** attempt
                logging.warning(f"Retrying fetch (attempt {attempt + 1}) in {wait_time}s due to error: {e}")
                sleep(wait_time)
        else:
            logging.error(f"Failed to fetch PMIDs after {max_retries} retries for query starting at index {start_index}")
            break
    return pmids

def read_query_from_file(file_path):
    return "cancer[Title]"

def get_fixed_month_interval(start_date):
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


def generate_date_batches(query, start_date, stop_date, target_papers_per_batch=10000, window_sizes=[365, 180, 90, 60, 40, 20, 10, 1], max_workers=2, verbose=True, full_text=False):
    date_ranges = []
    current = start_date
    total_papers = 0

    if verbose:
        print(f"Generating date batches from {start_date} to {stop_date}...\n")

    def get_count_for_window(start, days):
        end = min(start + timedelta(days=days), stop_date)
        q = f'("{start}"[PDat] : "{end}"[PDat]) {query}'
        count = get_pubmed_count(q, full_text=full_text)
        return (start, end, count)

    while current < stop_date:
        fixed_months = get_fixed_month_interval(current)
        if fixed_months:
            batch_end = min(current + relativedelta(months=fixed_months), stop_date)
            count = get_pubmed_count(f'("{current}"[PDat] : "{batch_end}"[PDat]) {query}', full_text=full_text)
            total_papers += count
            date_ranges.append((current, batch_end))
            current = batch_end + timedelta(days=1)
            continue

        candidates = [min(current + timedelta(days=w), stop_date) for w in window_sizes]
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(get_count_for_window, current, (end - current).days) for end in candidates]
            results = []
            for f in as_completed(futures):
                try:
                    results.append(f.result())
                except Exception as e:
                    logging.error(f"Error while evaluating date window: {e}")

        best = None
        for start, end, count in sorted(results, key=lambda x: (x[1] - x[0]).days, reverse=True):
            if verbose:
                print(f"{start} to {end} → {count} papers")
            if count <= target_papers_per_batch:
                best = (start, end)
                total_papers += count
                break

        if best:
            # Check if the batch is a single day and exceeds 10,000 papers
            if (best[1] - best[0]).days == 0 and count > 10000:
                logging.warning(f"1-day batch from {best[0]} contains {count} papers, exceeding the 10,000 threshold.")
            
            date_ranges.append(best)
            current = best[1] + timedelta(days=1)
        else:
            fallback_end = min(current + timedelta(days=1), stop_date)
            logging.warning(f"No acceptable window found at {current}, forcing fallback window.")
            count = get_pubmed_count(f'("{current}"[PDat] : "{fallback_end}"[PDat]) {query}', full_text=full_text)
            total_papers += count
            date_ranges.append((current, fallback_end))
            current = fallback_end + timedelta(days=1)

    if verbose:
        print(f"\nTotal batches generated: {len(date_ranges)}")
        print(f"~Total papers across all batches: {total_papers}")

    return date_ranges

def fetch_pmids_parallel(date_batches, query, full_text, max_workers):
    results = []
    metadata = []

    def fetch_single_batch(batch):
        start, end = batch
        date_query = f'("{start}"[PDat] : "{end}"[PDat])'
        full_query = f"{date_query} {query}"
        pmids = get_list(full_query, full_text)
        return pmids, {"start": str(start), "end": str(end), "count": len(pmids)}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_batch = {executor.submit(fetch_single_batch, batch): batch for batch in date_batches}
        for future in as_completed(future_to_batch):
            try:
                pmids, meta = future.result()
                results.extend(pmids)
                metadata.append(meta)
            except Exception as e:
                logging.error(f"Failed to fetch PMIDs for batch: {e}")

    return results, metadata

def fetch_pmids_over_period(query_file, start="2000-01-01", stop=None, full_text=False, max_workers=2, plot=True):
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read query.")
        return np.array([]), []

    if stop is None:
        stop = datetime.now().strftime("%Y-%m-%d")

    start_date = date.fromisoformat(start)
    stop_date = date.fromisoformat(stop)

    date_batches = generate_date_batches(query, start_date, stop_date, full_text=full_text, max_workers=max_workers)
    if plot:
        batch_metadata_preview = [
            {
                "start": str(start),
                "end": str(end),
                "count": get_pubmed_count(f'("{start}"[PDat] : "{end}"[PDat]) {query}', full_text=full_text)
            }
            for start, end in date_batches
        ]
        plot_density_over_time(batch_metadata_preview)
    pmid_list, batch_metadata = fetch_pmids_parallel(date_batches, query, full_text=full_text, max_workers=max_workers)
    pmid_list = list(set(pmid_list))  # Remove duplicates
    return np.array(pmid_list), batch_metadata




def plot_density_over_time(batch_metadata):
    if not batch_metadata:
        print("No batch metadata to plot.")
        return

    mid_dates = []
    densities = []

    for meta in batch_metadata:
        start = datetime.fromisoformat(meta["start"])
        end = datetime.fromisoformat(meta["end"])
        duration = (end - start).days or 1
        mid = start + (end - start) / 2
        mid_dates.append(mid)
        densities.append(meta["count"] / duration)

    plt.figure(figsize=(12, 6))
    plt.plot(mid_dates, densities, marker='o', linestyle='-', color='blue')
    plt.title("PMID Density Over Time")
    plt.xlabel("Date")
    plt.ylabel("Papers per Day")
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()