Folder: mhags-pipeline-fm (Female to Male Transplant Pipeline)  

####### Content:
1. 7 python scripts + 1 test.py 
2. 1 DockerFile
3. 1 req.txt 

####### Description: 
The pipeline takes in Donor and Patient WES samples, Sex, Patient HLA typing (tested on HLA-Matched Donor and Patients) - performs Variant Calling using DeepVariant, Variant Annotation using Funcotator, Merges Donor and Patient annotated data simultaneous *construction of Donor and Patient Custom Proteome, Gene/Variant Filtration steps: 1) Retains Variants specific to Patients, 2) Looks for Non-Synonymous Variants (Protein AA seq altering variants - Missense, Non-stop, Nonsense, Frame Shifts), 3) Looks for Variants present in "Tissues of Interest" (Lung, Liver, Skin, Colon, Lacrimal and Oral Mucosa - **scRNA Seq exp) ***AND ADDS ALL THE PEPTIDES ARISING FROM 9 MALE Y-SEX CHROMOSOME (gcloud path: gs://mhags-data/ychrom_9genes_peptides_for_blast.fa) ****Peptides - checks if remnant peptides are present in Donor or Patient Custom Proteomes, if not present then predicts binding stability using HLAthena with the Patient HLAs (Peptide-HLA binding prediction), HLAthena post processing: picks out peptides that have a binding stability of <0.5

#######
* construction of Donor and Patient Custom Proteome = Explanation in presentation

** scRNA Seq exp = For every tissue (organ) = 
1. Removed immune cells 
2. Reclustered and annotated tissue specific cell types 
3. Strategically eliminated drop-out effect 
(Path to scRNA Exp: mhags-pipeline-work/sc_analysis) 

*** MSY Peptides: Explanation on slide 27 to 29

**** Peptides = Peptides are taken from Funcotator Transcript subfolder (after Variant annotation)

####### Storage details and other info:
1. Docker name: "nidhihookeri/minors_pipeline_fm" (docker.com)
2. TERRA Workspace name: broad-fireclous-wuclonal/mHAgs_pipeline_nidhi_FM
3. TERRA Workflow name/Firecloud: mhags_pipeline_fm
4. GCP Bucket info: wu-lab-archives/wld-instance-2/mhags-data (gcloud compute ssh wld-instance-2 --project wu-lab-archives --zone us-central1-a) 
5. Presentation: mhags-pipeline-work/mhags_pipeline.pptx
6. Final WDL: mhags-pipeline-work/WDL/mhags_pipeline_fm.11.wdl
7. Associated inputs JSON:  mhags-pipeline-work/WDL/mhags_pipeline_fm.11.json

####### ORDER OF RUN (IF LOCALLY): 
1. bmt-simulation.py (includes as child files: functions.py, helpers.py, translation_helpers.py)
2. blast_preprocessing.py
3. blast_postprocessing.py
4. ychrom_blast_postops.py
5. HLAthena_preprocessing.py
6. HLAthena_postprocessing.py
7. HLAthena_postProcessing_ychrom.py

1. 
python3 mhags-pipeline/bmt-simulation.py -donor_vcf /path/to/donorVCF/from/Funcotator -host_vcf /path/to/hostVCF/from/Funcotator -host_sex M/F -donor_sex M/F -tissue AML/CML/CLL -gvhd yes/no 

OUTPUT CASES TAKEN INTO CONSIDERATION = (F->M): 
A. AML_specific_discordances.txt
B. GvHD_specific_discordances.txt
C. Ychrom_minor_antigens.txt (doesn't output in the final outputs section in the pipeline: static file)

2. 
2A. python3 mhags-pipeline/blast_preprocessing.py -specific_discordances /for/GvL/AML_specific_discordances.txt 
OUTPUT: 
A. peptides_for_blast.fa (FASTA file for GvL output)

2B. python3 mhags-pipeline/blast_preprocessing.py -specific_discordances /for/GvL/GvHD_specific_discordances.txt 
OUTPUT: 
A. peptides_for_blast.fa (FASTA file for GvHD output)

3. (Maintain sequence of Donor to Host - important)
# DONOR GvL
3Aa.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /for/GvL/AML_specific_discordances.txt -blastp_result /after/makeblastdb/and/blastp/blastp_result.csv -sample donor_basename -sample_type "donor" -condition GvL

3Ab MAKEBLASTDB AND BLASTP OF GvL PEPTIDES

# HOST GvL
3Ac.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /discordances/from/3Aa -blastp_result /after/makeblastdb/and/blastp/after/donor/blastp_result.csv -sample host_basename -sample_type "host" -condition GvL

# DONOR GvHD
3Ba.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /for/GvHD/GvHD_specific_discordances.txt -blastp_result /after/makeblastdb/and/blastp/blastp_result.csv -sample donor_basename -sample_type "donor" -condition GvHD

3Bb MAKEBLASTDB AND BLASTP FOR GvHD PEPTIDES

# HOST GvHD
3Bc  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /discordances/from/3Ba -blastp_result /after/makeblastdb/and/blastp/after/donor/blastp_result.csv -sample host_basename -sample_type "host" -condition GvHD

# Y CHROM ANTIGENS
3Ca  MAKEBLASTDB AND BLASTP OF YCHROM ANTIGENS - HOST
3Cb  python3 ychrom_blast_postops.py -blast_discordances /output/from/3Ca -sample host -switch True
3Cc  MAKEBLASTDB AND BLASTP OF YCHROM ANTIGENS - DONOR
3Cd  python3 ychrom_blast_postops.py -blast_discordances /output/from/3CC -sample donor -switch False

4. 
4A  python3 mhags-pipeline/HLAthena_preprocessing.py -discordances_after_blast /discordances/from/3Ac -autosomal "True"
4B  python3 mhags-pipeline/HLAthena_preprocessing.py -discordances_after_blast /discordances/from/3Bc -autosomal "True"
4C  python3 mhags-pipeline/HLAthena_preprocessing.py -discordances_after_blast /discordances/from/3Bc -autosomal "False"

5. HLAthena 

6. 
6A  python3 mhags-pipeline/HLAthena_postprocessing.py -sample_predictions /HLAthena/output/sample_predictions/for/GvL -discordances_after_blast /discordances/output/from/3Ac/for/GvL

6B  python3 mhags-pipeline/HLAthena_postprocessing.py -sample_predictions /HLAthena/output/sample_predictions/for/GvHD -discordances_after_blast /discordances/output/from/3Bc/for/GvHD

6C  python3 mhags-pipeline/HLAthena_postProcessing_ychrom.py -hlathena_predictions /HLAthena/output/sample_predictions/for/ychrom -sample_name file_output_name (variable)

FINAL OUTPUTS: 
A. GvL HLAthena postprocessing binding putative Minor Antigens (text file)
B. GvHD HLAthena postprocessing binding putative Minor Antigens (text file)
C. ychrom HLAthena postPorcessing binding putative Minor Antigens (text file)
