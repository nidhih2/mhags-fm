import pandas as pd
import argparse

def files_for_HLAthena(discordances, id_col, autosomal):
    discordances = pd.read_csv(discordances, sep="\t")
    peptides = discordances[['pep',id_col]]
    peptides = peptides.assign(source='host')
    peptides = peptides.drop_duplicates().reset_index().drop(columns='index')

    if autosomal=="True":
        ref_peptides = discordances[['ref_pep',id_col]]
        ref_peptides.columns=['pep',id_col]
        ref_peptides = ref_peptides.assign(source='donor')

        peptides = pd.concat([peptides,ref_peptides]).drop_duplicates().reset_index().drop(columns='index')

    peps_filename = 'HLAthena_peptides.txt'
    #peps_bucket_path = 'gs://'+BUCKET+'/'+case_id+'_peptides.txt'

    peptides.to_csv(peps_filename, sep='\t', index=False)

def main():
    files_for_HLAthena(discordances=discordances_after_blast, id_col="transid", autosomal=autosomal)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-discordances_after_blast', help = "discordances after blast search (allgenes_puma.txt)")
    parser.add_argument('-autosomal', help=" autosomal = True/False, True for GvL and GvHD, False for Y chromosome")

    args = parser.parse_args()
    discordances_after_blast = args.discordances_after_blast
    autosomal = args.autosomal
    main()
