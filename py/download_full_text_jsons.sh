mkdir -p ./Full_text_jsons

while IFS= read -r PMCID; do
    url="https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/${PMCID}/unicode"
    curl -s "${url}" > "./Full_text_jsons/${PMCID}.json"
done < full_text_pmc.txt
