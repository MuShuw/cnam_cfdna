# README

## Calculating the genome-wide, fragmentomic profiles

Script `src/00_BamToChrScore.py`.

<dl>
<dt>Input</dt>  <dd>Path to the BAM file</dd>
<dt>Output</dt> <dd>Files with genome-wide fragmentomic profiles, one per chromosome</dd>
</dl>

Command line

```bash
python 00_BamToChrScore.py -i path/to/BH01.bam
```

The resulting `*.npz` files are then moved to the directory `01_chr_score`.

## Calculating the gene-centric, fragmentomic profiles

Script `src/01_ScorePerSite.py`.

<dl>
<dt>Input</dt>  <dd>Path to a folder containing the chromosome profiles</dd>
<dt>Output</dt> <dd>2 files containing the FD and WPS scores around TSSs of a given GeneClass (-gc)</dd>
</dl>

Coding genes

```bash
python src/01_ScorePerSite.py -i data/01_chr_score/ -o data/02_tss_score/ -gc coding_mRNA -ft5 data/90_fantom/intermediary/
```

Long intergenic, non-coding RNAL

```bash
python src/01_ScorePerSite.py -i data/01_chr_score/ -o data/02_tss_score/ -gc lncRNA_intergenic -ft5 data/90_fantom/intermediary/
```

## 02 Features extraction

The script `02_FeatsExtraction.py` calculate the features for each genes

<dl>
<dt>Input</dt> <dd>Gene-centric WPS and FD table for a either class (coding or intergenic lncRNAs)</dd>
<dt>Output</dt> <dd>A file containing the features for each gene of a specific class</dd>
</dl>

```bash
python 02_FeatsExtraction.py \
    -iw data/02_tss_score/WPS_coding_mRNAscoreTable.tsv \
    -if data/02_tss_score/FD_coding_mRNAscoreTable.tsv  \
    -o  data/03_features/mRNA_features.tsv

python 02_FeatsExtraction.py \
    -iw data/02_tss_score/WPS_lncRNA_intergenicscoreTable.tsv \
    -if data/02_tss_score/FD_lncRNA_intergenicscoreTable.tsv  \
    -o  data/03_features/lincRNA_features.tsv
```

## Model training

A *model configuration* is given by choosing

* between coding genes or intergenic lncRNA genes
* an operational definition of gene class A (eg, top 100 genes based on positive fold-changes)
* an operational definition of gene class B (eg, ubiquitous genes)
* a set of fragmentomic features to include

We test many model configurations.

The script `032_parallel_ClassificationScoring.sh` is used to train the gene classification model for each FANTOM tissue and each model configurations.

### LincRNA configurations

```bash
bash 032_parallel_ClassificationScoring.sh    \
    data/03_features/lincRNA_features.tsv     \ # linc features (from previous step)
    data/04_SVC_performances/lincrna/         \ # output directory
    100                                       \ # number of genes in each class
    data/90_fantom/onto/                      \ # ontologies tables
    data/03b_control_features/Leuko_feats.tsv \ # control features 1
    data/03b_control_features/Ubi_feats.tsv   \ # features of ubiquitous genes
```

### Running mRNA scoring

sudo bash src/032_parallel_ClassificationScoring.sh data/03_features/mRNA_features.tsv data/04_SVC_performances/mrna/ 100 data/90_fantom/ontho/ data/03b_control_features/Leuko_feats.tsv data/03b_control_features/Ubi_feats.tsv

## 03a Supervised learning 

group="PFb"
feats="Ulz"
sort="mann_whitney_pval"
rank=100
NAME="$group"_"$sort"_"$feats"
python3 src/03b_ClassificationScoring.py -r $rank -i data/90_fantom/ontho/ -g $group -f $feats -if data/03_features/lincRNA_features.tsv -lbf data/03b_control_features/Leuko_feats.tsv -u data/03b_control_features/Ubi_feats.tsv -sb mann_whitney_pval -o data/04_SVC_performances/"$NAME"_svc_Intermediary -v


## 041

## Joining SVM output in one folder
bash src/041_join_SVM_output.sh

## 042

### Calculating the means et medians of lymphomyeloid and non-lymphomyeloid ontologies

python3 src/042_MeanMedianCacl.py -i data/04_SVC_performances/ -o data/042_SVC_performances_means/
