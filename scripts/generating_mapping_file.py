import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from Bio import SeqIO
from tqdm import tqdm 
from pyteomics import parser as peptide_parser
import requests
import re
from mygene import MyGeneInfo
import pickle
from io import StringIO
import os
import argparse

# ------------------------------------------------------------
# Functions 
# ------------------------------------------------------------

# this function convert gene name to RefSeq ids, so that the following mutalyzer can recognize
def gene_to_refseq(gene_symbol):
    if pd.isna(gene_symbol):
        return None
    elif re.match(r'[A-Z]{1,3}_[0-9]', gene_symbol):
        match = re.search(r'[A-Z]{1,3}_[0-9]+\.[0-9]', gene_symbol)
        if match:
            return {gene_symbol: match.group()}
    else:
        mg = MyGeneInfo()
        result = mg.query(gene_symbol, species="human", fields="refseq.rna")
        try:
            rna_id = result["hits"][0]["refseq"]["rna"]
            return {gene_symbol: rna_id}
        except:
            return {gene_symbol: None}


# this function is to retriev information from mutalyzer based on hgvs names
def query_mutalyzer(hgvs):
    url = f"https://mutalyzer.nl/api/normalize/{hgvs}"
    try:
        response = requests.get(url)
        data = response.json()
        return {
            "DNA_hgvs": data['input_description'],
            "protein_hgvs": data['protein']['description'],
            "mutation_type": len(data['input_model']['variants']),
            "wild_type_protein": re.sub(r'\*$','',data['protein']['reference']),
            "mutated_protein": re.sub(r'\*$','',data['protein']['predicted'])
        }
    except:
        return None


# this function is to make unique peptides 
def unique_trypsin(wildtype, variant, missed_cleavages):
    if not wildtype or not variant:
        return None
    wt_peps = peptide_parser.cleave(wildtype, peptide_parser.expasy_rules['trypsin'], missed_cleavages=missed_cleavages)
    mt_peps = peptide_parser.cleave(variant, peptide_parser.expasy_rules['trypsin'], missed_cleavages=missed_cleavages)
    return [pep for pep in mt_peps if pep not in wt_peps]

# # this functino is to make unique peptides 
# def compute_uniq_peptide(df, source_col, new_col):
#     df = df.copy()
#     df[new_col] = None
#     df[new_col] = df[new_col].astype(object)
#     for i in df.index:
#         query = df.loc[i,source_col]
#         if not isinstance(query, list):
#             df.at[i, new_col] = None
#             continue
#         other_peptides = list(set(df.drop(index=i)[source_col].dropna().sum()))
#         uniq = [p for p in query if (p not in other_peptides) and (len(p) >= 6)]
#         df.at[i, new_col] = uniq
#     return df


def compute_uniq_peptide(df, source_col, new_col):
    df = df.copy()
    df[new_col] = None
    df[new_col] = df[new_col].astype(object)
    for i in df.index:

        query = df.loc[i,source_col]
        mutation_type = int(df.loc[i, 'mutation_type'])

        if not isinstance(query, list):
                df.at[i, new_col] = None
                continue
    
        # for single point mutation, we aim to use single unique peptide
        if mutation_type == 1:    
            other_peptides = list(set(df.drop(index=i)[source_col].dropna().sum()))
            uniq = [p for p in query if (p not in other_peptides) and (len(p) >= 6)]
            df.at[i, new_col] = uniq

        # for double point mutation, we aim to use multiple unique peptides (each peptide not necessary have mutations or unique)
        elif mutation_type == 2:
            same_peptide_exist = df.drop(index=i)[source_col].apply(lambda x: set(query) <= set(x) if isinstance(x, list) and isinstance(query, list) else False).any()
            # means the peptide combination is not unique, we will not use it
            if same_peptide_exist:
                df.at[i, new_col] = None
                continue
            # means the peptide combination is unique, we will use it 
            df.at[i, new_col] = query
    return df




# this function remove the leading M in the protein sequence 
def remove_leading_M(protein):
    return re.sub('^M','', protein)

missed_cleavages, output_dir, n_jobs = 0, "output", 14

def generating(missed_cleavages, output_dir, n_jobs):
    # ------------------------------------------------------------
    # Step 1: Read IthaGenes Table
    # ------------------------------------------------------------
    url = 'https://www.ithanet.eu/db/ithagenes?action=list&hcat=0b-'
    tables = pd.read_html(url)
    df = tables[0]

    # ------------------------------------------------------------
    # Step 2: Convert Gene Symbols to RefSeq RNA IDs
    # ------------------------------------------------------------


    # Extract gene names and get RefSeq IDs
    df['gene'] = df['HGVS Name'].str.split(':').str[0]
    genes = df['gene'].dropna().unique()
    refid_list = Parallel(n_jobs=n_jobs)(delayed(gene_to_refseq)(gene) for gene in genes)

    # Merge dictionaries
    refid_dict = {}
    for entry in refid_list:
        if entry:
            refid_dict.update(entry)

    # Create new columns for refid and full HGVS
    df['refid'] = df['gene'].apply(lambda g: refid_dict.get(g) if pd.notna(g) else None)
    df['mutation'] = df['HGVS Name'].str.split(':').str[1]
    df['hgvs'] = df['refid'] + ':' + df['mutation']


    # ------------------------------------------------------------
    # Step 3: Query Mutalyzer for Protein Consequences (slow to run, ~ 20 min)
    # ------------------------------------------------------------

    mutalyzer_results = Parallel(n_jobs=n_jobs)(
        delayed(query_mutalyzer)(x) for x in tqdm(df['hgvs'])
    )
    df['mutalyzer'] = mutalyzer_results

    # # save a temporary file
    # import pickle
    # file = 'step3.pkl'
    # with open(file, 'wb') as f:
    #     pickle.dump(df, f)

    # file = 'step3.pkl'
    # with open(file, 'rb') as f:
    #     df = pickle.load(f)

    
    # filter out rows with no mutalyzer information
    df = df[df['mutalyzer'].notna()]
    
    # ------------------------------------------------------------
    # Step 4: Identify Unique Tryptic Peptides
    # ------------------------------------------------------------

    # get unique peitde between mutated proteins and wildtype proteins 
    df['uniq_ref'] = df['mutalyzer'].apply(
        lambda x: unique_trypsin(remove_leading_M(x['wild_type_protein']), remove_leading_M(x['mutated_protein']), missed_cleavages=missed_cleavages) if x else None
    )

    df['uniq_ref'] = df['uniq_ref'].apply(lambda x: [p for p in x if len(p) >= 6] if isinstance(x, list) else None)

    # ------------------------------------------------------------
    # Step 5: Filter Unique Peptides Compared to Other Variants
    # ------------------------------------------------------------

    df['mutation_type'] = df['mutalyzer'].apply(lambda x: x['mutation_type'] if pd.notna(x) else None)
    df = compute_uniq_peptide(df, source_col='uniq_ref', new_col='uniq_var')

    # ------------------------------------------------------------
    # Step 6: Compare to Human Proteome
    # ------------------------------------------------------------

    url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(proteome:UP000005640)+AND+(reviewed:true)"
    response = requests.get(url)
    fasta_io = StringIO(response.text) 
    records = list(SeqIO.parse(fasta_io, "fasta"))

    human_prot = []
    for x in records:
        human_prot.append({
            'header': x.description, 
            'id': x.id,
            'sequence': str(x.seq)
        })

    human_prot = pd.DataFrame(human_prot)
    human_prot['pep'] = human_prot['sequence'].apply(lambda x: list(peptide_parser.cleave(x, peptide_parser.expasy_rules['trypsin'], missed_cleavages=missed_cleavages)))
    human_pep_pool = set(human_prot['pep'].sum())
    
    df['uniq_human'] = None
    df['uniq_human'] = df['uniq_human'].astype(object)
    for i in df.index:
        query = df.loc[i, 'uniq_var']
    
        if not query:
            df.at[i, 'uniq_human'] = []
            continue

        mutation_type = df.loc[i, 'mutation_type']

        if mutation_type == 1:
            # for single point mutation, we only need to check if the peptide is in human proteome
            df.at[i, 'uniq_human'] = [p for p in query if p not in human_pep_pool and len(p) >= 6]
        elif mutation_type == 2:
            mask = all([p for p in query if p not in human_pep_pool ])
            if mask:
                # for double point mutation, we need to check if all peptides are not in human proteome
                df.at[i, 'uniq_human'] = query
            else:
                # if any peptide is in human proteome, we will not use this variant
                df.at[i, 'uniq_human'] = []

    # ------------------------------------------------------------
    # Step 7: Make Mapping Table for Peptide → Variant
    # ------------------------------------------------------------
    mapping = df[['Hb Name', 'uniq_human','mutation_type']].explode('uniq_human')
    mapping = mapping[mapping['uniq_human'].notna()]
    mapping = mapping[mapping['Hb Name'].notna()]
    mapping.to_csv(os.path.join(output_dir, 'mapping.csv'), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generating mapping file", allow_abbrev=False)
    # required arguments
    parser.add_argument('--missed_cleavage', required=True, type=int, choices=[0,1], help="Aollowed missed cleavage with trypsin, higher value leads to more unique peptides.")
    parser.add_argument('--output_dir', required=True, help="Output directory, where mapping file will be generated.")
    parser.add_argument("--threads", required=True, type=int, help="Number of threads.")
    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)    
    generating(missed_cleavages=args.missed_cleavage, output_dir=args.output_dir, n_jobs=args.threads)
