## Execution Steps
1. Run main().py; `python __main__.py`
   
2.	Input Data:
  - Provide email address and an NCBI API key.
  - The script reads a PubMed "query" from file in the same directory.
  - Start and stop date.
3.	Fetching PMIDs:
    - The script fetches PMIDs over a specified date range using the provided query.
4.	Fetching Articles:
    - For each PMID, the script fetches the corresponding article's metadata and saves it into a CSV file.
5.	Adding Publisher Information:
    - The script enhances the CSV file by adding publisher information, leveraging DOIs and ISSNs.
6.	Saving Results:
  - The final CSV file, enriched with publisher data, is saved, process logged in pipe.log.
## Usage Scenario
This script is for researchers and data scientists who need to gather scientific article data from PubMed for analysis, literature reviews, or metadata enrichment.
