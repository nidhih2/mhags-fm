import pandas as pd
import argparse

def process_ychrom_predictions(hlathena_predictions, sample_name):
	hlathena_predictions = pd.read_csv(hlathena_predictions, sep="\t")
	hlathena_predictions = hlathena_predictions[['pep', 'transid', 'len', 'assign.MSi_ranks', 'assign.MSi_allele']]

	hlathena_weak = hlathena_predictions[hlathena_predictions.loc[:, 'assign.MSi_allele'] != "unknown"]
	hlathena_weak = hlathena_predictions[hlathena_predictions.loc[:, 'assign.MSi_ranks'] < 0.5]

	hlathena_weak_transid = hlathena_weak["transid"].str.split(";", expand=True)[0]
	gene_name = hlathena_weak_transid.str.split(":", expand = True)[0]
	hlathena_weak["gene"] = gene_name

	peptide_pos = hlathena_weak_transid.str.split(":", expand = True)[2]
	peptide_start = peptide_pos.str.split("_", expand=True)[0]
	peptide_end = peptide_pos.str.split("_", expand = True)[1]
	hlathena_weak["Peptide-Start-Pos"] = peptide_start
	hlathena_weak["Peptide-End-Pos"] = peptide_end

	hlathena_weak = hlathena_weak.drop_duplicates()

	hlathena_weak = hlathena_weak[["gene", "Peptide-Start-Pos", "Peptide-End-Pos", "pep", "len", 'assign.MSi_allele', 'assign.MSi_ranks']]
	hlathena_weak.columns = ["Gene", "Peptide-Start-Pos", "Peptide-End-Pos", "Y-Chr-Peptide", "Peptide-Length", "assign.MSi_allele", "assign.MSi_ranks"]
	hlathena_weak = hlathena_weak.sort_values(by = ["Peptide-Length", "assign.MSi_allele"])

	hlathena_weak.to_csv(sample_name+"_host_ychrom_weak_peptides.txt", sep="\t", index=False)

	#hlathena_weak_merge = yhost_annotations.merge(hlathena_weak, how = "left")
	#hlathena_weak_merge.to_csv("yhost_annotations_with_predictions.txt", sep="\t", index = False)

	gene_count = hlathena_weak["Gene"].value_counts().to_dict()
	allele_count = hlathena_weak["assign.MSi_allele"].value_counts().to_dict()
	final_dict = {"Tot_peptides_Y_chr":len(set(hlathena_weak["Y-Chr-Peptide"])), "gene_count":[gene_count], "allele":[allele_count]}
	final_dict_df = pd.DataFrame(final_dict, index=[0])
	final_dict_df.to_csv(sample_name+"_ychrom_stats.csv", index=False)

	return(hlathena_weak)

def main():
	process_ychrom_predictions(hlathena_predictions = hlathena_predictions, sample_name = sample_name)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-hlathena_predictions")
	parser.add_argument("-sample_name")
	#parser.add_argument("-yhost_annotations")
	args = parser.parse_args()
	hlathena_predictions = args.hlathena_predictions
	sample_name = args.sample_name
	#yhost_annotations = args.yhost_annotations
	main()












