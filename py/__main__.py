def main():
    # Define the query file and dates
    query_file = "query"  # Filename containing the PubMed query
    start_date = "2000-01-01"
    stop_date = "2024-08-01"

    # Step 1: Fetch PMIDs over the specified period using a query
    pmid_array = fetch_pmids_over_period(query_file, start=start_date, stop=stop_date)

    # Save PMIDs if any are fetched
    if len(pmid_array) > 0:
        save_pmids(pmid_array)

        # Step 2: Fetch article data from PubMed and save to CSV
        csv_file = "pmid_pmc_title_journal_doi_issn.csv"
        save_articles_to_csv(pmid_array.tolist(), csv_file)

        # Step 3: Add publishers using DOIs and ISSNs
        output_csv = 'FINAL.csv'
        add_publishers_to_csv(csv_file, output_csv)
        logging.info("Process completed.")
    else:
        logging.error("No PMIDs were fetched. Check query or date range.")

email = "m.n.khanji@umcg.nl"
API_KEY = "70faf5cc42501a814dcc4bdb1862acaf3909"

if __name__ == "__main__":
    main()



