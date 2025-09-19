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
    Entrez.email = email

    plot_choice = input("Do you want to plot the PMID trend? Press Enter to continue or type 'skip' to skip this step: ").strip().lower()
    oa_choice = input("Do you want to flag Open-Access articles? Press Enter to continue or type 'skip' to skip this step: ").strip().lower()

    predatory_choice = input("Do you want to flag suspected predatory papers? skin self hating are commin ").strip().lower() or None

    intact_choice = input("Do you want to add IntAct/Corum interaction flags? Press Enter to continue or type 'skip' to skip this step: ").strip().lower()
    
    complex_choice = input("Do you want to add Complex Portal interaction flags? Press Enter to continue or type 'skip' to skip this step: ").strip().lower()
    complex_dir = None
    if complex_choice != 'skip':
        complex_dir = input("Enter the path to your Complex Portal input directory (or press Enter to skip): ").strip() or None

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
    

    # flag papers of sus predatory journals or publishers 

    if predatory_choice != 'skip':
        print("Adding predatory listed flags")
        data = DR.flag_predatory_publications(data)



    # Add Open Access flag
    if oa_choice != 'skip':
        print("Adding Open Access flags...")
        data = DR.filter_df_oa_database(data)

    # Add IntAct/Corum flags
    if intact_choice != 'skip':
        print("Adding IntAct/Corum interaction flags...")
        data = DR.add_intact_corrum(data, pos_file="intact.txt", neg_file="intact_negative.txt")

    # Add Complex Portal flags
    if complex_choice != 'skip' and complex_dir:
        print("Adding Complex Portal interaction flags...")
        data = DR.add_complex_portal(input_dir=complex_dir, data=data)


    # Save to CSV
    filename = f"NCBI_data_{datetime.now().strftime('%Y%m%d')}.csv"
    data.to_csv(filename, index=False)
    print(f"Data saved to {filename}")

if __name__ == "__main__":
    main()



