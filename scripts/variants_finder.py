import pandas as pd 
import argparse
import os 

def mapping():
    parser = argparse.ArgumentParser(description="Mapping sample peptides to varaint peptides", allow_abbrev=False)
    # required arguments
    parser.add_argument('--peptide_table', required=True, help="Path to the trypsin peptide table (csv format). Columns are samples.")
    parser.add_argument('--index_column', required=True, help='The index column in peptide_table. It should be peptide amino acids.')
    parser.add_argument('--mapping_file', required=True, help="Path to mapping file (csv).")
    parser.add_argument('--output_dir', required=True, help="Output directory where final mapped variants will be located.")
    parser.add_argument('--variant_min_level', required=True, default=0, type=float, help="The minimal abundance level for variants considered to be present")
    args = parser.parse_args()
    
    # read in mapping file 
    mapping = pd.read_csv(args.mapping_file)

    os.listdir(os.getcwd())
    mapping = pd.read_csv('./output/mapping.csv')

    # read in peptide table 
    dat = pd.read_csv(args.peptide_table, index_col=args.index_column)
    
    # find varaints peptides in sampes 
    def var_finder(df):
        # samples in columns, index is the peptides sequence in each sample
        result = []
        for sample_name in df.columns:
            sample = df[sample_name].dropna()
            if args.variant_min_level == 0:
                sample_pep_pool = list(set(sample.index[sample > 0]))
            else:
                sample_pep_pool = list(set(sample.index[sample >= args.variant_min_level]))
            sample_pep_pool = list(set(sample.index[sample > args.variant_min_level]))
            tmp = mapping[mapping['uniq_human'].isin(sample_pep_pool)].copy()
            tmp['Sample'] = sample_name
            tmp.columns = ['Hb Name', 'Peptides', 'Sample']
            result.append(tmp)
        return pd.concat(result, axis=0)

    # save data 
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    df = var_finder(dat)
    df.to_csv(os.path.join(args.output_dir, 'mapped_variants.csv'), index=False)
    print("\nVariants mapping completed !\n")

if __name__ == '__main__':
    mapping()
