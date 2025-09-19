#import built-in libraries
import os
from typing import List, Dict, Tuple, Optional
import requests
import logging
import urllib
import json
import time
from datetime import datetime
import numpy as np 
import csv
import re
import random
#import third-party libraries
from concurrent.futures import ThreadPoolExecutor, as_completed
import matplotlib.pyplot as plt
import pandas as pd 
from tqdm.notebook import tqdm
from ratelimit import limits, sleep_and_retry
#Import DR module from Functions folder
from Functions import DataRetrieval as DR
#import biopython/Eutils libraries
from Bio import Entrez, Medline
from eutils import EutilsNCBIError, EutilsRequestError

def main():


    print("""  _____        _               _____      _        _                 _ 
 |  __ \      | |             |  __ \    | |      (_)               | |
 | |  | | __ _| |_ __ _ ______| |__) |___| |_ _ __ _  _____   ____ _| |
 | |  | |/ _` | __/ _` |______|  _  // _ \ __| '__| |/ _ \ \ / / _` | |
 | |__| | (_| | || (_| |      | | \ \  __/ |_| |  | |  __/\ V / (_| | |
 |_____/ \__,_|\__\__,_|      |_|  \_\___|\__|_|  |_|\___| \_/ \__,_|_|
                                                                       
                                                                       """)
    # === Step 1: Collect all inputs at the start ===
    query_file = input("Enter the path to your query file (.txt): ").strip()
    
    start_date = input("Enter the start date (YYYY-MM-DD) or press Enter to skip: ").strip() or None
    stop_date = input("Enter the stop date (YYYY-MM-DD) or press Enter to skip: ").strip() or None
    email = input("Enter your email for NCBI Entrez: ").strip()
    api_key = input("Enter the NCBI API_KEY (recommended) or press Enter to skip: ").strip() or None
    Entrez.email = email
    Entrez.api_key = api_key


    # Configure max retries on failed requests (optional, default is 3)
    Entrez.max_tries = 5

    plot_choice = input("Do you want to plot the PMID trend? Press Enter to continue or type 'skip' to skip this step: ").strip().lower()
    
    # === Step 2: Run functions conditionally ===

    # Fetch PMIDs
    plot = plot_choice != 'skip'
    print("Fetching PMIDs...")
    pmid_list = DR.fetch_pmids_over_period(
        query_file=query_file,
        start=start_date,
        stop=stop_date or None,
        plot=plot,
        full_text=False,
    )

    # Fetch and parse metadata
    print("Fetching PubMed metadata...")
    data = DR.fetch_parse_pubmed_metadata(pmid_list)


    
    # Process publishers
    print("Processing publishers...")
    data = DR.process_publishers_concurrent(data, email)
    

    # Save to CSV
    filename = f"NCBI_data_{datetime.now().strftime('%Y%m%d')}.csv"
    data.to_csv(filename, index=False)
    print(f"Data saved to {filename}")

if __name__ == "__main__":
    main()