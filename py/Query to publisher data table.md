Execution Steps
1.	Input Data:
o	Provide email address and an NCBI API key.
o	The script reads a PubMed query from file in the same directory.
o	Start and stop date.
2.	Fetching PMIDs:
o	The script fetches PMIDs over a specified date range using the provided query.
3.	Fetching Articles:
o	For each PMID, the script fetches the corresponding article's metadata and saves it into a CSV file.
4.	Adding Publisher Information:
o	The script enhances the CSV file by adding publisher information, leveraging DOIs and ISSNs.
5.	Saving Results:
o	The final CSV file, enriched with publisher data, is saved, and the process is logged.
Usage Scenario
This script is for researchers and data scientists who need to gather scientific article data from PubMed for analysis, literature reviews, or metadata enrichment.
