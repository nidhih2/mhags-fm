import pandas as pd 

# Kari's file - has 52 gene/possible peptides which we don't need 
ychrom = pd.read_csv("gs://mhags-data/ychrom_8_11mer_peptides.txt", sep="\t") 
print("Peptides",len(set(ychrom["pep"])))
# ychrom["gene"].drop_duplicates().to_csv("just_the_ychrom_genes.txt", sep="\t", index=False)

# these are the 9 Y chr genes we are interested in 
ychrom_needed_genes = ["KDM5D", "UTY", "USP9Y", "ZFY", "DDX3Y",	"EIF1AY", "RPS4Y1", "NLGN4Y", "TMSB4Y"]
print(len(ychrom_needed_genes))

# filtering out the 9 from 52 gene list
ychrom_needed_peptides = ychrom[ychrom.gene.isin(ychrom_needed_genes)]
print(ychrom_needed_peptides)
print(ychrom_needed_peptides.gene.drop_duplicates())
ychrom_needed_peptides.to_csv("gs://mhags-data/ychrom_9genes_8_11mer_peptides.txt", sep="\t", index=False)
ychrom_needed_peptides.to_csv("ychrom_9genes_8_11mer_peptides.txt", sep="\t", index=False)
# [31298 rows x 5 columns] 9 genes with possible 31298 peptides 

# next file : ychrom_peptides_fasta_gen.py -> for generation of fasta file of all the peptides 



