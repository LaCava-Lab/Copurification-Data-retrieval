from functools import partial
from tenacity import retry, stop_after_delay, wait_fixed, retry_if_exception_type

# Decorator 1 
retry_on_communication_error = partial(
    retry,
    stop=stop_after_delay(10),  # Maximum 10 seconds wait.
    wait=wait_fixed(0.4),  # Wait 400ms between retries
    retry=retry_if_exception_type((Exception,))  # Use a tuple for exceptions
)

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

@retry_on_communication_error()
def get_list(query):
    """Retrieve all PMIDs for a given query using the PubMedFetcher."""
    num_of_articles = 500
    start_index = 0
    pmids = []
    while True:
        pmid_batch = fetcher.pmids_for_query(query,
                                            retstart=start_index,
                                            retmax=num_of_articles,
                                            pmc_only=True)
        pmids.extend(pmid_batch)
        start_index = len(pmids)
        if len(pmid_batch) < num_of_articles:
            break
    return pmids

@retry_on_communication_error()
def fetch_pmids_over_period(query_file, start="2000-01-01", stop=None):
    """Fetch PMIDs over a specified period using a query read from a file."""
    query = read_query_from_file(query_file)
    if not query:
        logging.error("Failed to read query.")
        return np.array([])

    if stop is None:
        stop = datetime.now().strftime("%Y-%m-%d")

    start_date_str = start
    pmid_list = []

    while True:
        if date.fromisoformat(start_date_str) <= date.fromisoformat("2002-07-01"):
            month_interval = 6
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2005-11-01"):
            month_interval = 5
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2009-11-01"):
            month_interval = 4
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2011-10-01"):
            month_interval = 3
        elif date.fromisoformat(start_date_str) <= date.fromisoformat("2023-01-01"):
            month_interval = 2
        else:
            month_interval = 4

        next_start = date.fromisoformat(start_date_str) + relativedelta(months=month_interval)
        end_date = (next_start - relativedelta(days=1))
        end_date_str = end_date.strftime('%Y-%m-%d')

        date_str = f'''(("{start_date_str}"[Date - Publication] : "{end_date_str}"[Date - Publication]) '''
        pmids = get_list(date_str + query)
        pmid_list.extend(pmids)
        start_date_str = next_start.strftime('%Y-%m-%d')
        if next_start >= date.fromisoformat(stop):
            break

    # Remove duplicates by converting to a set, then back to a list
    pmid_clean_list = list(set(pmid_list))
    logging.info(f"Total PMIDs fetched: {len(pmid_clean_list)}")

    return np.array(pmid_clean_list)
