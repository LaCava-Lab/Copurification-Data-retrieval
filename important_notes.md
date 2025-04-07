### Metapub is a python wrapper for NCBI that we leverage to fetch unique Identifiers of papers from a certain query.

### NCBI has a 10k hit cap per query, we had to make any query fetch in a X month interval to make sure the papers fetched are less than 10K per request.

### NCBI does not include publisher name in it's metadata, Using DOI and ISSN, we retrieved publisher data from Xreff API with *Habanero* a python wrapper for Xreff.

### To retrieve full text available for datamining, BioC, a restfull API provides full text of papers in a simple and annotated format in Json/XML and Unicode/Ascii under 3 categories:
1. Commercial use.
2. Non-commercial use.
3. Auther Manuscript Collection.

## Supplementary data situation:
[This](https://www.researchgate.net/profile/Senay-Kafkas/publication/271194007_Database_citation_in_supplementary_data_linked_to_Europe_PubMed_Central_full_text_biomedical_articles/links/54c0ef260cf28a6324a45661/Database-citation-in-supplementary-data-linked-to-Europe-PubMed-Central-full-text-biomedical-articles.pdf?__cf_chl_tk=Jny2wF0BmKRwZ4ykyCCOkoFC3xeYpbcJNiNromsWGa8-1744050286-1.0.1.1-Z2uL644EuYcY05jviE2AMJuS1xubMKuoc06gzwEtguo) paper uses europePMC to download supplements 