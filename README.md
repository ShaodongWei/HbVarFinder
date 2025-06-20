# HbVarFinder

# when to use 
This script is to help quickly find hemoglobin structural variants based on tryptic proteomics. So far the tool only support one point mutation (n=890).

# how to use 
Before using, it is recommended to re-run the preparation.py script to re-generate the variants table, so that the table is up to date. The varaints information 
is obtained from the [IthaNet](https://www.ithanet.eu/db/ithagenes?action=list&hcat=0b-)

# Install packages 
´´´
python -m venv hbvarfinder 
source ./hbvarfinder/bin/activate 
source hbvarfinder/Scripts/activate #windows
pip install --upgrade pip 
pip install -r ./script/packages.txt

´´´
