
#######################################################################
import re
from metapub import PubMedFetcher, PubMedArticle, pubmedcentral
import logging
import time
import numpy as np
import pandas as pd
import requests
import json
import csv
import os 
import urllib.parse
from time import sleep
from tqdm import tqdm  
from datetime import datetime, timedelta, date
from functools import partial
from dateutil.relativedelta import relativedelta
from urllib.parse import urlparse, parse_qs
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type, stop_after_attempt, wait_exponential
from multiprocessing.pool import ThreadPool
from typing import List, Iterator, Optional, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from eutils import EutilsNCBIError, EutilsRequestError
import xml.etree.ElementTree as ET
from urllib.error import HTTPError
import random
#######################################################################


# Decorator 1 
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(2),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)


fetcher = PubMedFetcher()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############Papers DF######################
def fetch_article(pmcid: str) -> PubMedArticle:
    """Fetch a single article from PubMed by PMCID PMCs should be in the format 'PMCXXXXXX'."""
    article = fetcher.article_by_pmcid(pmcid)
    return article

@retry_on_communication_error()
def fetch_articles(pmcids: List[str], *, processes: Optional[int] = 1) -> Iterator[PubMedArticle]:
    """Fetch multiple articles from PubMed in parallel using a thread pool + tqdm."""
    with ThreadPool(processes=processes) as pool:
        # Wrap with tqdm for progress tracking
        for article in tqdm(
            pool.imap_unordered(fetch_article, pmcids),
            total=len(pmcids),
            desc="Fetching PubMed articles",
            unit="article"
        ):
            if article is not None:
                yield article

def fetch_articles_meta(pmids: List[str]) -> pd.DataFrame:
    """Fetch articles and return them as a pandas DataFrame."""
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

#########################################Fetch PMID from PMCID##########################################################################################################
def fetch_pmid(pmcid):
    """Fetch PMID from PMCID safely with fallback."""
    # Primary: try Metapub
    try:
        from metapub import PubMedFetcher
        fetch = PubMedFetcher()
        article = fetch.article_by_pmcid(pmcid)
        if article and article.pmid:
            return article.pmid
    except Exception as e:
        logging.warning(f"Metapub failed for {pmcid}: {e}")

    # Fallback: use NCBI ID conversion API directly
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={pmcid}&format=json"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        records = data.get("records", [])
        if records and "pmid" in records[0]:
            return records[0]["pmid"]
        else:
            logging.warning(f"No PMID found for {pmcid} in fallback API.")
            return None
    except Exception as e:
        logging.error(f"Fallback API failed for {pmcid}: {e}")
        return None


def get_pmid_for_otherid(pmcid_clean_list):
    """Converts a list of PMCIDs to PMIDs in parallel"""
    PMIDs = {}
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_pmcid = {executor.submit(fetch_pmid, pmcid): pmcid for pmcid in pmcid_clean_list}
        for future in as_completed(future_to_pmcid):
            pmcid = future_to_pmcid[future]
            try:
                pmid = future.result()
                PMIDs[pmcid] = pmid
            except Exception as e:
                logging.error(f"Error processing PMCID {pmcid}: {e}")
                PMIDs[pmcid] = None
    return PMIDs

#################################################################################################################################################################################################
########################################Full text JSON download (BIOC)##################################################################################################################################################
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

###############################################################################################################################################################################################################################
########################  download Supplementary_xml (BIOC) ##########################################################################################################################################################################

def download_supplementary_materials(oa_pmcids, output_dir="supplementary_materials", format="bioc_xml", delay=0.5):
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




################################################Extract text from BioC XML and JSON to DataFrame##########################################################################################################
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

###############################################Extract text from BioC JSON to DataFrame##########################################################################################################

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

                    pmc = base_filename
                    pmid = fetch_pmid(pmc) or "NOPMID"

                    entry_serial = 1

                    def get_entry_id():
                        nonlocal entry_serial
                        entry_id = f"{pmid}_{entry_serial:03}"
                        entry_serial += 1
                        return entry_id

                    for passage in passages:
                        infons = passage.get('infons', {})
                        section_type = infons.get('section_type')
                        subtitle = infons.get('type')
                        if section_type in section_types:
                            norm_text = passage.get('text', '').strip()
                            if norm_text:
                                data_rows.append({
                                    "EntryID": get_entry_id(),
                                    "filename": pmc,
                                    "section_type": section_type,
                                    "subtitle": subtitle,
                                    "text": norm_text
                                })

                    if xml_directory:
                        xml_path = os.path.join(xml_directory, base_filename + ".xml")
                        if os.path.isfile(xml_path):
                            xml_df = extract_text_from_bioc_xml(xml_path)
                            for _, row in xml_df.iterrows():
                                text = row['text']
                                if text:
                                    data_rows.append({
                                        "EntryID": get_entry_id(),
                                        "filename": pmc,
                                        "section_type": "SUPPLEMENT",
                                        "subtitle": row.get('type', None),
                                        "text": text.strip()
                                    })

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




#############################################full pipeline#############################################


def pmc_full_pipeline(
    oa_pmcids: List[str],
    section_types: List[str] = None,
    base_output_dir: str = "./pmc_data",
    download_supplementary: bool = True,
    delay: float = 0.5
) -> pd.DataFrame:
    """Pipeline: Download PMC JSON + supplementary XML, extract text, and save outputs."""

    section_types = section_types or ["ABSTRACT", "INTRO", "METHODS", "DISCUSSION", "TITLE", "CONCL", "FIG", "DISCUSS", "RESULTS"]

    json_dir = os.path.join(base_output_dir, "Full_text_jsons")
    xml_dir = os.path.join(base_output_dir, "supplementary_materials")

    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(xml_dir, exist_ok=True)

    print("\n=== STEP 1: Fetching article metadata ===")
    papers_df = fetch_articles_meta(oa_pmcids)
    papers_df.to_csv(os.path.join(base_output_dir, "Papers.csv"), index=False)
    print(f"Saved article metadata to {os.path.join(base_output_dir, 'Papers.csv')}")

    print("\n=== STEP 2: Downloading full-text JSONs ===")
    success_count, failed_json = download_pmc_articles(oa_pmcids, output_dir=json_dir)

    unsaved_supps = []
    if download_supplementary:
        print("\n=== STEP 3: Downloading supplementary materials ===")
        unsaved_supps = download_supplementary_materials(oa_pmcids, output_dir=xml_dir, format="bioc_xml", delay=delay)

    print("\n=== STEP 4: Extracting content into DataFrame ===")
    df = extract_text_from_json_to_dataframe(json_dir, section_types, xml_directory=(xml_dir if download_supplementary else None))

    df.to_csv(os.path.join(base_output_dir, "Full_text.csv"), index=False)
    print(f"Saved full-text data to {os.path.join(base_output_dir, 'Full_text.csv')}")

    df.attrs['summary'] = {
        'total_pmcids': len(oa_pmcids),
        'json_downloaded': success_count,
        'json_failed': len(failed_json),
        'supp_downloaded': len(oa_pmcids) - len(unsaved_supps) if download_supplementary else 'skipped',
        'supp_failed': len(unsaved_supps) if download_supplementary else 'skipped',
        'entries_extracted': len(df)
    }

    print("\n=== SUMMARY ===")
    for k, v in df.attrs['summary'].items():
        print(f"{k}: {v}")

    return df