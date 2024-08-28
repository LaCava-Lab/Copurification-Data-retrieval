#%%Fetch PMCS list avlbl for full text  

import chat_ip as ip

def fetch_PMC_fulltext(): 
    # Step 1: Read the query from a file
    query_file = "query" 
    query = ip.read_query_from_file(query_file)
    if not query:
        ip.logging.error("Query reading failed. Exiting.")                                                               
        return

    # Step 2: Fetch PMIDs over a specified period
    start_date = "2000-01-01"
    stop_date = None  # Will default to the current date if None
    pmid_array = ip.fetch_pmids_over_period(query_file, start=start_date, stop=stop_date)
    if pmid_array.size == 0:
        ip.logging.error("No PMIDs fetched. Exiting.")
        return

    # Step 3: Save the PMIDs to files
    ip.save_pmids(pmid_array)

    # Step 4: Retrieve PMCIDs for the fetched PMIDs
    pmc_id_list = ip.get_pmcid_for_otherid(pmid_array)
    pmc_ids_filename = "PMCIDS.csv"
    ip.pd.DataFrame(pmc_id_list, columns=["PMCID"]).to_csv(pmc_ids_filename, index=False)
    ip.logging.info(f"PMCIDs saved to {pmc_ids_filename}")

    # Step 5: Filter the OA database using the retrieved PMCIDs
    oa_file_list = "oa_file_list.csv"  # OA file list CSV filename
    pmc_ids_filename = "PMCIDS.csv"
    filtered_df = ip.filter_oa_database(oa_file_list, pmc_ids_filename)

    ip.logging.info("OA database filtering completed.")
    return filtered_df

    # Define the query file and dates and if U want PMC only or not(boolean)
API_KEY = "70faf5cc42501a814dcc4bdb1862acaf3909"


fetch_PMC_fulltext()
