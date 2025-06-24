# HbVarFinder

# When to use 
This script is to help quickly find hemoglobin structural variants based on tryptic proteomics. The var_finder.py script reports the existence of variants if larger than zero abundance of variant-specific peptides are detected in samples. 


# Install packages 
```
python -m venv hbvarfinder  
source ./hbvarfinder/bin/activate  # linux
source hbvarfinder/Scripts/activate #windows 
pip install --upgrade pip  
pip install -r ./scripts/packages.txt
```

# How to use

## Prepare mapping file 
Before using, it is recommended to re-run the generating_mapping_file.py script to re-generate the variants table, so that the table is up to date. The varaints information 
is obtained from the [IthaNet](https://www.ithanet.eu/db/ithagenes?action=list&hcat=0b-)

```
python ./scripts/generating_mapping_file.py -h
python ./scripts/generating_mapping_file.py --missed_cleavage 0 --output_dir output --threads 10

```

## Find variant peptides in samples 
Now we have a ready-to-use mapping file, then we identify variants in our samples
```
python ./scripts/variants_finder.py --peptide_table data_in/peptide_table.csv --index_column peptide --mapping_file output/mapping.csv --output_dir output
```
