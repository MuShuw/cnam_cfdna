# README

## Calculating the genome-wide, fragmentomic profiles

Script `src/00_BamToChrScore.py`.

## Calculating the gene-centric, fragmentomic profiles

Script `src/src/01_ScorePerSite.py`.

<dl>
<dt>Input</dt>  <dd>Path to a folder containing the chromosome scoring</dd>
<dt>Output</dt> <dd>2 files containing the FD and WPS scores around TSSs of a given GeneClass (-gc)</dd>
</dl>

Coding genes

```bash
python3 src/01_ScorePerSite.py -i data/01_chr_score/ -o data/02_tss_score/ -gc coding_mRNA -ft5 data/90_fantom/intermediary/
```

Long intergenic, non-coding RNAL

```bash
python3 src/01_ScorePerSite.py -i data/01_chr_score/ -o data/02_tss_score/ -gc lncRNA_intergenic -ft5 data/90_fantom/intermediary/
```

## 02 Features generation
The script 02_FeatsExtraction.py calculate the features for each genes

input : TSS's WPS and FD table for a GeneClass 
output : 1 file containing the features for each gene of a specific GeneClass

to be done for both coding_mRNA and lncRNA_intergenic

usage :
python3 src/02_FeatsExtraction.py -iw data/02_tss_score/WPS_coding_mRNAscoreTable.tsv -if data/02_tss_score/FD_coding_mRNAscoreTable.tsv -o data/03_features/mRNA_features.tsv
python3 src/02_FeatsExtraction.py -iw data/02_tss_score/WPS_lncRNA_intergenicscoreTable.tsv -if data/02_tss_score/FD_lncRNA_intergenicscoreTable.tsv -o data/03_features/lincRNA_features.tsv

# 03a SVC scoring 


group="PFb"
feats="Ulz"
sort="mann_whitney_pval"
rank=100
NAME="$group"_"$sort"_"$feats"
python3 src/03b_ClassificationScoring.py -r $rank -i data/90_fantom/ontho/ -g $group -f $feats -if data/03_features/lincRNA_features.tsv -lbf data/03b_control_features/Leuko_feats.tsv -u data/03b_control_features/Ubi_feats.tsv -sb mann_whitney_pval -o data/04_SVC_performances/"$NAME"_svc_Intermediary -v

# 03b parallel SVC scoring 


## Running linc scoring 
sudo bash src/032_parallel_ClassificationScoring.sh data/03_features/lincRNA_features.tsv data/04_SVC_performances/lincrna/ 100 data/90_fantom/ontho/ data/03b_control_features/Leuko_feats.tsv data/03b_control_features/Ubi_feats.tsv

## Running mRNA scoring
sudo bash src/032_parallel_ClassificationScoring.sh data/03_features/mRNA_features.tsv data/04_SVC_performances/mrna/ 100 data/90_fantom/ontho/ data/03b_control_features/Leuko_feats.tsv data/03b_control_features/Ubi_feats.tsv

# 041

## Joining SVM output in one folder
bash src/041_join_SVM_output.sh

# 042

## Calculating the means et medians of lymphomyeloid and non-lymphomyeloid ontologies

python3 src/042_MeanMedianCacl.py -i data/04_SVC_performances/ -o data/042_SVC_performances_means/


