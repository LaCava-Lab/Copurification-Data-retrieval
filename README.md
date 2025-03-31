# Copurification-Data-retrieval

### A. Get started 
*Perform the following steps only the first time, in order to be able to run the code in this repository at your computer (locally).*

**Mac OS**
1. Clone the repository locally accoring to github instructions, go to Terminal (Mac OS) and and have Python installed (v.3.10.10 or earlier).
2. Make sure you have installed virtualenv(or use any other environment or container you are comfortable with):

   `pip3 install --upgrade pip`
   
	`pip3 install virtualenv`
	
3. Create the virtualenv and activate it (skip this step and substitute with appropriate steps, if you used another solution)

	`python3 -m venv venvChatIP`
	
	`source venvChatIP/bin/activate`
4. Build the environment.
 
	`pip3 install --upgrade pip`
	
	`pip3 install ipykernel`
	
	`python3 -m ipykernel install --name=venvChatIP --user`
	
	`brew install wget` # install  wget with brew
	
	`sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"` #  install e-direct with wget

	**Windows VS Code** 
	e-direct is not supported on windows, download wsl2 instead
	`wsl --install` # install command in PowerShell 


5. Install the NCBI oa_text mining DB and move it to the right folder:

	`wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv`
	`mv oa_file_list.csv RefDocs/oa_file_list.csv`

6. ##### After cloning the repo:

Create a `LOGS/` + `Reference_files/` folders and add your `keys.py` and `query`:
   `mkdir Reference_files`
   `mkdir LOGS/`
   `touch Reference_files/keys.py`# ie. API_KEY = "70fath27ifk51a814dcc4bdb1862acaf3909"
   `touch Reference_files/query`

### C. Details
#### Folders



#### Perform specific tasks

Run examples of how to use the code in this repo for different tasks can be seen in the notebook `Notebooks/Examples.ipynb`. Wrappers in bash for specific pipelines coming up.

| Function name | Task | Location | Dependencies | Output(s) |
| :------------- | :---- | :-------- | :------------ | :------ |
| read_query_from_file | read a custom query from a text a certain file and return a query string variable| Functions/DataRetrieval.py | Dependancies | query |
| get_list | | Get a list of hit PMIDs from the query you have | Functions/DataRetrieval.py | Dependancies | PMIds List |
| fetch_pmids_over_period | Fetch PMIDs over a specified period using a query |  Functions/DataRetrieval.py | Dependancies | clean PMIDs np.array | 
| fetch_article | Fetch a single article and return its data as a dict |  Functions/DataRetrieval.py | Dependancies | single article data | 
| fetch_articles_to_dataframe | Fetch a group of articles by leveraging func(fetch_article) | Functions/DataRetrieval.py | Dependancies | Group of article data | 
| fetch_pmcid | converts PMID to it corresponding PMCID | Functions/DataRetrieval.py | Dependancies | PMCID |
| get_pmcid_for_otherid | converts a list of PMIDs to a list of corresponding PMCIDs in parallel using func(fetch_pmcid) multithreading | Functions/DataRetrieval.py | Dependancies | PMCIDs list |
| filter_oa_database | filters a CSV database of PMC articles to only include IDs that are open access | Functions/DataRetrieval.py | pandas (as pd), oa_file_list.csv | Pandas Series of matching "Accession ID" values