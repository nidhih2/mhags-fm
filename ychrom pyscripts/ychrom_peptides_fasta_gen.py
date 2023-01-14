import pandas as pd

ychrom_peptides = pd.read_csv("gs://mhags-data/ychrom_9genes_8_11mer_peptides.txt", sep="\t")

with open('ychrom_9genes_peptides_for_blast.fa','a') as f:
    for index, row in ychrom_peptides.iterrows():
        if len(row.pep) >= 8 & len(row.pep) <= 11:
            f.write('>'+row.pep_id)
            f.write("\n")
            f.write(row.pep)
            f.write("\n")
f.close()

# file present on gcloud bucket : gs://mhags-data/ychrom_9genes_peptides_for_blast.fa