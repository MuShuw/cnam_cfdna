# README

For this work, we relied on
Python (3.6.9),
NumPy (1.19.5),
Pandas (1.1.5),
Scikit-Learn (0.24.2) and
Pysam (0.19.1).

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
* an operational definition of gene class A (eg, top 100 genes with positive fold-changes)
* an operational definition of gene class B (eg, ubiquitous genes)
* a set of fragmentomic features to include

We test many model configurations.

The script `032_parallel_ClassificationScoring.sh` is used to train the gene classification model for each FANTOM ontology and each model configurations. It relies, internally, on the `03b_ClassificationScoring.py` script, which trains an SVM model for each model configuration (see below).

### LincRNA configurations

```bash
bash 032_parallel_ClassificationScoring.sh     \
     data/03_features/lincRNA_features.tsv     \ # linc features (from previous step)
     data/04_SVC_performances/lincrna_alt/     \ # output directory
     100                                       \ # number of genes in each class
     data/90_fantom/onto/                      \ # ontologies
     data/03b_control_features/Leuko_feats.tsv \ # features with leukocyte p-values/FC
     data/03b_control_features/Ubi_feats.tsv     # features of ubiquitous genes
```

### Running mRNA scoring

```bash
bash 032_parallel_ClassificationScoring.sh     \
     data/03_features/mRNA_features.tsv        \ # coding features (from previous step)
     data/04_SVC_performances/mrna/            \ # output directory
     100                                       \
     data/90_fantom/onto/                      \
     data/03b_control_features/Leuko_feats.tsv \
     data/03b_control_features/Ubi_feats.tsv
```

### Classification model

The classification model used in this work is a kernel support vector machine (implemented by Scikit-Learn).

The configuration is parametered by the variables `cat`, `group`, `sort` and `feats`:

```bash
cat="data/03_features/lincRNA_features.tsv"
group="PFb"
sort="mann_whitney_pval"
feats="Ulz"
rank=100
NAME="$group"_"$sort"_"$feats"

python3 03b_ClassificationScoring.py -v            \
    -r $rank                                       \ # rank=100
    -i data/90_fantom/onto/                        \
    -g $group                                      \ # group="PFb"
    -sb $sort                                      \ # sort="mann_whitney_pval" 
    -f $feats                                      \ # feats="Ulz"
    -if $cat      				   \
    -lbf data/03b_control_features/Leuko_feats.tsv \
    -u data/03b_control_features/Ubi_feats.tsv     \
    -o data/04_SVC_performances/"$NAME"_svc_Intermediary # output file
```

In this example:

* Genes are sorted by $p$-value
* Class A: 100 positive fold-change genes at the top of the list (most significant $p$-values)
* Class B: 100 genes at the bottom of the list (least significant)
* Training uses the feature set know as ‘Ulz’

## Postprocessing

## Joining SVM output in one folder

bash src/041_join_SVM_output.sh

## 042

### Calculating the means et medians of lymphomyeloid and non-lymphomyeloid ontologies

python3 src/042_MeanMedianCacl.py -i data/04_SVC_performances/ -o data/042_SVC_performances_means/

## Creating the groups


```bash
mrna="data/03_features/mRNA_features.tsv"
linc="data/03_features/lincRNA_features.tsv"


rank=100
onto="data/90_fantom/onto/"
leuko="data/03b_control_features/Leuko_feats.tsv"
ubi="data/03b_control_features/Ubi_feats.tsv"

sort=(mann_whitney_pval fold mann_whitney_pval fold)
table=($mrna $mrna $linc $linc)
rna=(coding_mRNA coding_mRNA lncRNA_intergenic lncRNA_intergenic)

for i in 0 1 2 3; do
    NAME=ARN"${rna[$i]}"_sortedBy_"${sort[$i]}"
    python3 src/90_GroupsForProfil.py -v               \
        -r $rank                                       \ # rank=100
        -i $onto                                       \ # ontologies folder
        -o $NAME                                       \ # output name
        -if ${table[$i]}                               \ # target RNA features
        -lbf $leuko				       \ # leukocyte bot feats 
        -u $ubi					       \ # fantom ubi feats
        -sb ${sort[$i]}                                \ # sort
```

## Supplementary Information

### Bokeh interactive plots

```bash
cd data/042_SVC_performances_means
bokeh serve --show ../../src/051_generate_interactivePlot.py
```
