#!/bin/bash


if [ -d "venvChatIP" ];

# Activate environment
then 
source venvChatIP/bin/activate
pip3 install --upgrade pip
echo | python3 --version

# Create and activate environment
else 
pip3 install --upgrade pip
pip3 install virtualenv 
python3 -m venv venvChatIP

source venvChatIP/bin/activate
pip3 install --upgrade pip

pip3 install ipykernel
python3 -m ipykernel install --name=venvChatIP --user

#sh -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" # install e-direct with curl
brew install wget # install  wget with brew
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"  install e-direct with wget
fi

export PATH=${HOME}/edirect:${PATH}
export NCBI_API_KEY=e401e127b0eefb86674c5dba3a009dd26508

echo | python3 --version
	
pip3 install -r requirements.txt
#There is not yet an option to include argument --no-deps with a requirements file see: https://github.com/pypa/pip/pull/10837
pip3 install 
pip3 install xformers==0.0.21 --no-deps

# Activate Jupyter

jupyter notebook