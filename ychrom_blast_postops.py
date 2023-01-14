import pandas as pd 
import argparse

# we do a blast search of all the possible Y chr peptides against donor proteome and 
# see for peptides that are 100% matching to any peptide in the donor and remove them and 
# and the resulting peptides from above against host proteome to remove any 100% matching peptides 

def ychrom_blast_PostOps(blast_discordances, sample):
    blast_results = pd.read_csv(blast_discordances, sep=",", header=None) # result of blastp 
    blast_results.columns=['protein_id','peptide_id','protein_len','peptide_len','protein_start','protein_end','peptide_start','peptide_end','protein_seq','peptide_seq','evalue','score','alignment_length','p_ident','n_ident']
    blast_results.to_csv("just2.txt", sep="\t", index=False)
    blast_results = blast_results.assign(protein_gene=blast_results.protein_id.str.split('|', expand=True)[0])
    blast_results = blast_results.assign(peptide_gene=blast_results.peptide_id.str.split(':', expand=True)[0])

    blast = blast_results[~((blast_results.peptide_len==blast_results.alignment_length) & (blast_results.p_ident==100))]
    blast = blast[["protein_id", "peptide_id", "peptide_len", "peptide_seq", "peptide_gene"]]
    blast.columns = ["protein_id", "transid", "len", "pep", "gene"]
    blast = blast.drop_duplicates()
    print(blast.shape)
    blast.to_csv(sample+"_ychrom_blastPostOps_peptides.txt", sep="\t", index=False) # req outputs in WDL : donor_ychrom_blast_peptides.txt and host_ychrom_blast_peptides.txt
    return blast


def donor_peptides_fasta_gen(host_peps, switch):
    with open("ychrom_hostPostBlast_peptides.fa",'w') as f:
        for index, row in host_peps.iterrows():
            if len(row.pep) >= 8 & len(row.pep) <= 11:
                f.write('>'+row.transid)
                f.write("\n")
                f.write(row.pep)
                f.write("\n")
    f.close()


def main():
    ychrom_blast = ychrom_blast_PostOps(blast_discordances=blast_discordances, sample=sample)
    donor_peptides_fasta_gen(host_peps=ychrom_blast, switch=switch)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    parser.add_argument("-blast_discordances")
    parser.add_argument("-sample")
    parser.add_argument("-switch", required=False)

    args = parser.parse_args()
    blast_discordances = args.blast_discordances
    sample = args.sample
    switch = args.switch
    main()