import app_data as DR
import os 



def main():

    print("""_____        _               _____      _        _                 _ 
            |  __ \      | |             |  __ \    | |      (_)               | |
            | |  | | __ _| |_ __ _ ______| |__) |___| |_ _ __ _  _____   ____ _| |
            | |  | |/ _` | __/ _` |______|  _  // _ \ __| '__| |/ _ \ \ / / _` | |
            | |__| | (_| | || (_| |      | | \ \  __/ |_| |  | |  __/\ V / (_| | |
            |_____/ \__,_|\__\__,_|      |_|  \_\___|\__|_|  |_|\___| \_/ \__,_|_|
""")
                                                                                  
    base_output_dir = './paper_files'
    os.makedirs(base_output_dir, exist_ok=True)


    # Step 1: Load PMCIDs from text file
    with open('pmcids.txt', 'r') as f:
        pmcids = [line.strip() for line in f if line.strip()]


    # Step 2: Fetch metadata
    papers_df = DR.fetch_articles_meta(pmcids)
    papers_csv_path = os.path.join(base_output_dir, 'Papers.csv')
    papers_df.to_csv(papers_csv_path, index=False)
    print(f"Saved metadata to {papers_csv_path}")


    # Step 3: Download JSON full texts
    json_dir = os.path.join(base_output_dir, 'Full_text_jsons')
    DR.download_pmc_articles(pmcids, output_dir=json_dir)


    # Step 4: Download supplementary materials
    supp_dir = os.path.join(base_output_dir, 'supplementary_materials')
    os.makedirs(supp_dir, exist_ok=True)


    # download_supplementary_materials
    try:
        DR.download_supplementary_materials(pmcids, output_dir=supp_dir)
    except Exception as e:
        print(f"Skipping supplementary downloads due to error: {e}")


    # Step 5: Extract text into DataFrame
    section_types = ["ABSTRACT", "INTRO", "METHODS", "DISCUSSION", "TITLE", "CONCL", "FIG", "DISCUSS", "RESULTS"]
    df = DR.extract_text_from_json_to_dataframe(json_dir, section_types, xml_directory=supp_dir)


    # Step 6: Save extracted text
    fulltext_csv_path = os.path.join(base_output_dir, 'Full_text.csv')
    df.to_csv(fulltext_csv_path, index=False)
    print(f"Saved extracted full text to {fulltext_csv_path}")


if __name__ == '__main__':main()
main()