#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

OUTPUT_DIR="./pubmed_articles"
QUERY='immunoprecipitation AND ("2025/07/28"[PDAT] : "2025/07/30"[PDAT])'
LOG_FILE="./pubmed_download.log"
CSV_OUTPUT_FILE="./pubmed_metadata_detailed.csv"

create_output_dir() {
    mkdir -p "$OUTPUT_DIR"
}

initialize_files() {
    : > "$LOG_FILE"
    # Write CSV header
    echo "PMID,Title,Date Received,Date Revised,Date Accepted,Date Completed,Abstract,Journal,PublicationYear,Authors,DOI,ISSN,ChemicalList_Names,MeshHeadingList_Descriptors" > "$CSV_OUTPUT_FILE"
}

log_warning() {
    local message="$1"
    local timestamp;
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    printf "[%s] %s\n" "$timestamp" "$message" >> "$LOG_FILE"
}

fetch_pubmed_ids() {
    local ids;
    if ! ids=$(esearch -db pubmed -query "$QUERY" | efetch -format uid); then
        log_warning "Failed to retrieve PubMed IDs"
        printf "Error: Failed to retrieve PubMed IDs\n" >&2
        return 1
    fi

    if [[ -z "${ids// /}" ]]; then
        log_warning "No PubMed articles found for query: $QUERY"
        printf "No PubMed articles found matching the query\n" >&2
        return 1
    fi

    printf "%s\n" "$ids"
}

count_ids() {
    local id_list="$1"
    local count;
    count=$(printf "%s\n" "$id_list" | grep -c '^[0-9]\+$')
    if [[ -z "$count" || "$count" -eq 0 ]]; then
        log_warning "Counted zero valid PubMed IDs"
        printf "No valid PubMed IDs to process\n" >&2
        return 1
    fi
    printf "Total articles to download: %d\n" "$count"
}

sanitize_csv_field() {
    local input="$1"
    # Escape quotes and remove newlines
    printf '%s' "$input" | sed ':a;N;$!ba;s/\n/ /g' | sed 's/"/""/g'
}

# Helper function to clean element content (remove element names and quotes)
clean_element_content() {
    local content="$1"
    if [[ -n "$content" ]]; then
        # Remove element name prefixes like "ArticleTitle ", "AbstractText ", "ArticleId "
        # Remove surrounding quotes
        echo "$content" | sed 's/^[^"]*"//; s/"[^"]*$//; s/^[[:space:]]*//; s/[[:space:]]*$//'
    else
        echo ""
    fi
}

# Helper function to clean date format
clean_date_format() {
    local date_input="$1"
    if [[ -n "$date_input" ]]; then
        # Remove spaces around hyphens and trim whitespace
        echo "$date_input" | sed 's/[[:space:]]*-[[:space:]]*/-/g; s/^[[:space:]]*//; s/[[:space:]]*$//'
    else
        echo ""
    fi
}

process_article_metadata() {
    local pmid="$1"
    local temp_xml="${OUTPUT_DIR}/PMID${pmid}.tmp.xml"

    # Fetch XML for one article
    if ! efetch -db pubmed -id "$pmid" -format xml > "$temp_xml" 2>tmp_err.log; then
        local error_msg;
        error_msg=$(<tmp_err.log)
        log_warning "EFETCH failed for PMID$pmid: $error_msg"
        rm -f tmp_err.log "$temp_xml"
        return 1
    fi

    if [[ ! -s "$temp_xml" ]]; then
        log_warning "EFETCH returned empty result for PMID$pmid"
        rm -f "$temp_xml"
        return 1
    fi



    # PMID (already known, but extract for consistency)
    local extracted_pmid
    extracted_pmid=$(xtract -input "$temp_xml" -pattern PubmedArticle -element MedlineCitation/PMID 2>/dev/null || echo "$pmid")

    # Title - extract content and clean it
    local title
    title=$(xtract -input "$temp_xml" -pattern ArticleTitle -element . 2>/dev/null || echo "")
    title=$(clean_element_content "$title")

    # Journal info
    local journal pub_year issn
    journal=$(xtract -input "$temp_xml" -pattern Journal -element Title 2>/dev/null || echo "")
    pub_year=$(xtract -input "$temp_xml" -pattern PubDate -element Year 2>/dev/null || echo "")
    issn=$(xtract -input "$temp_xml" -pattern Journal -element ISSN 2>/dev/null || echo "")

    # Abstract - handle multiple sections properly and clean content
    local abstract
    abstract=$(xtract -input "$temp_xml" -pattern AbstractText -element . 2>/dev/null | paste -sd " " - || echo "")
    abstract=$(clean_element_content "$abstract")

    # DOI - specifically target main article DOI only and clean it
    local doi
    doi=$(xtract -input "$temp_xml" -pattern ArticleId -if '@IdType' -equals "doi" -element . 2>/dev/null | head -1 || echo "")
    doi=$(clean_element_content "$doi")


    #  DATE EXTRACTION NEEDS WORK!  

    # Dates from History (using correct pattern)
    local received_date revised_date accepted_date completed_date

    # Extract dates with proper element targeting
    received_date=$(xtract -input "$temp_xml" -pattern PubMedPubDate -if '@PubStatus' -equals "received" -sep "-" -element Year Month Day 2>/dev/null | clean_date_format || echo "")
    revised_date=$(xtract -input "$temp_xml" -pattern PubMedPubDate -if '@PubStatus' -equals "revised" -sep "-" -element Year Month Day 2>/dev/null | clean_date_format || echo "")
    accepted_date=$(xtract -input "$temp_xml" -pattern PubMedPubDate -if '@PubStatus' -equals "accepted" -sep "-" -element Year Month Day 2>/dev/null | clean_date_format || echo "")

    # DateRevised standalone
    local date_revised_standalone
    date_revised_standalone=$(xtract -input "$temp_xml" -pattern DateRevised -sep "-" -element Year Month Day 2>/dev/null | clean_date_format || echo "")
    if [[ -n "$date_revised_standalone" && -z "$revised_date" ]]; then
        revised_date="$date_revised_standalone"
    fi

    # DateCompleted
    completed_date=$(xtract -input "$temp_xml" -pattern DateCompleted -sep "-" -element Year Month Day 2>/dev/null | clean_date_format || echo "")

    # Authors - extract name components properly
    local authors
    authors=$(xtract -input "$temp_xml" -pattern Author -sep " " -element ForeName LastName 2>/dev/null | paste -sd ";" - || echo "")

    # Chemicals
    local chemicals
    chemicals=$(xtract -input "$temp_xml" -pattern Chemical -element NameOfSubstance 2>/dev/null | paste -sd "," - || echo "")

    # Mesh Headings
    local mesh_headings
    mesh_headings=$(xtract -input "$temp_xml" -pattern MeshHeading -sep ": " -element DescriptorName QualifierName 2>/dev/null | paste -sd ";" - || echo "")

    # ================================================
    # Sanitize all fields for CSV
    # ================================================
    extracted_pmid=$(sanitize_csv_field "$extracted_pmid")
    title=$(sanitize_csv_field "$title")
    received_date=$(sanitize_csv_field "$received_date")
    revised_date=$(sanitize_csv_field "$revised_date")
    accepted_date=$(sanitize_csv_field "$accepted_date")
    completed_date=$(sanitize_csv_field "$completed_date")
    abstract=$(sanitize_csv_field "$abstract")
    journal=$(sanitize_csv_field "$journal")
    pub_year=$(sanitize_csv_field "$pub_year")
    authors=$(sanitize_csv_field "$authors")
    doi=$(sanitize_csv_field "$doi")
    issn=$(sanitize_csv_field "$issn")
    chemicals=$(sanitize_csv_field "$chemicals")
    mesh_headings=$(sanitize_csv_field "$mesh_headings")

    # ================================================
    # Write to CSV
    # ================================================
    printf '"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s"\n' \
        "$extracted_pmid" "$title" "$received_date" "$revised_date" "$accepted_date" \
        "$completed_date" "$abstract" "$journal" "$pub_year" "$authors" "$doi" "$issn" \
        "$chemicals" "$mesh_headings" >> "$CSV_OUTPUT_FILE"

    # Cleanup temp file
    rm -f "$temp_xml"
}

main() {
    create_output_dir
    initialize_files

    local pubmed_ids;
    if ! pubmed_ids=$(fetch_pubmed_ids); then
        return 1
    fi

    if ! count_ids "$pubmed_ids"; then
        return 1
    fi

    local id;
    while IFS= read -r id; do
        [[ -z "${id// /}" ]] && continue
        if ! process_article_metadata "$id"; then
            printf "Warning: Failed to retrieve metadata for PMID%s (logged)\n" "$id" >&2
            continue
        fi
    done <<< "$pubmed_ids"
}

main "$@"