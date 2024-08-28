import chat_ip as ip

def main():
    ip.setup_logging() #setup logging and retry with already set parameters in PY
    ip.setup_retry()

    # Define the query file and dates and if U want PMC only or not(boolean)
    query_file = "query"  # Filename containing the PubMed query
    start_date = "2000-01-01"
    stop_date = "2024-08-01"

    # Step 1: Fetch PMIDs over the specified period using a query
    pmid_array = ip.fetch_pmids_over_period(query_file, start=start_date, stop=stop_date, pmc_only=False)

    # Save PMIDs if any are fetched
    if len(pmid_array) > 0:
        ip.save_pmids(pmid_array)

        # Step 2: Fetch article data from PubMed and save to CSV
        csv_file = "pmid_pmc_title_journal_doi_issn.csv"
        ip.save_articles_to_csv(pmid_array.tolist(), csv_file)

        # Step 3: Add publishers using DOIs and ISSNs
        output_csv = 'FINAL.csv'
        ip.add_publishers_to_csv(csv_file, output_csv, email="m.n.khanji@umcg.nl")
        ip.logging.info("Process completed.")
    else:
        ip.logging.error("No PMIDs were fetched. Check query or date range.")

API_KEY = "70faf5cc42501a814dcc4bdb1862acaf3909"

if __name__ == "__main__":
    main()



