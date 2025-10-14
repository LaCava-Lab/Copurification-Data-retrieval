


# trypipe.py

trypipe.py imports functions from app_data.py. and imports a list of pmcids from pmcids.txt in the same directory.

starts with pmcids.txt as input, downloads paper.json and then supplementary.xml for the pmcids then produces a df with  the agreed upon format with full text of papers as output.