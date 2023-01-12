library(data.table)
library(dplyr)
library(plotly)
library(wesanderson)
library(Rtsne)
library(htmltools)
library(R.utils)
set.seed(1609)
## connection to plotly
Sys.setenv("plotly_username"='muxxu')
Sys.setenv("plotly_api_key"="fO9e7r3LiZIkjpEMEkaH")
flag_t_dt_expression_plot = 0


# EXPRESSION DATAFRAME IMPORT AND PRE-PROCESSING

## Here we import the file containing the Tag Per Million
dt_expression = data.table::fread("../0002_tsne_onOptimal_results/FANTOM_CAT.expression_atlas.gene.lv3_robust.MEDIAN_tpm_PER_ONTHO_ELYAS.tsv.gz",
                                  stringsAsFactors = F)
## For futur filtering we import the IDs of the lincRNA and mRNA
### Import the tables containing the IDs
dt_linc_IDs = data.table::fread("../0002_tsne_onOptimal_results/supp_table_11.cell_type_gene_association_lincrna.tsv",
                                stringsAsFactors = F)
dt_mrna_IDs = data.table::fread("../0002_tsne_onOptimal_results/supp_table_11.cell_type_gene_association_mrna.tsv",
                                stringsAsFactors = F)
### Extract the IDs
linc_IDs = dt_linc_IDs$CAT_geneID
mrna_IDs = dt_mrna_IDs$CAT_geneID
geneIds = dt_expression$geneID

## We tranpose the expression matrix and name the cols and rows 
cellLines = colnames(dt_expression)[-1]
dt_expression$geneID = NULL
dt_expression = as.data.frame(t(dt_expression))
colnames(dt_expression)=geneIds
rownames(dt_expression)
dim(dt_expression)

## Eliminate null variation gene
### Ca  lculating the genes with var==0
var_fold_gene = apply(dt_expression, MARGIN = 2, var)
names(var_fold_gene) = colnames(dt_expression)
cst_gene = names(var_fold_gene[var_fold_gene==0])

### Filtering the constant genes from the expression table
### And updating the lincRNA and mRNA IDs
dt_expression<-setDT(dt_expression)
dt_expression = dt_expression[,-..cst_gene]
linc_IDs <- intersect(linc_IDs,colnames(dt_expression))
mrna_IDs <- intersect(mrna_IDs,colnames(dt_expression))

## We set two alternative table for lincRNA and mRNA expression tables
dt_expression_linc <- dt_expression[,..linc_IDs]
dt_expression_mrna <- dt_expression[,..mrna_IDs]

## We calculate the log of the TPM for the t-SNE
## You can choose the use one of the 3 intermidiary table
## all genes / mRNA genes / lincRNA genes
## But is see no reason not to use the full table for all
dt_expression_bis = log(dt_expression+0.0001)
# dt_expression_bis = log(dt_expression_linc+0.0001)
# dt_expression_bis = log(dt_expression_mrna+0.0001)
tsne <- Rtsne(as.matrix(dt_expression_bis),
              check_duplicates=FALSE,
              pca=TRUE,
              perplexity=30,
              theta=0.5,
              dims=2)
tsne

## Import a performance file
perf_filename='inputs_output_features_lncRNA:sort_mann_whitney_pval:feats_UlzW_mixed1_sny:group_tPFvsb:rank_100:_svc_Intermediary_simple.csv'

dt_perf = read.csv(perf_filename,
				   sep=',', row.names=1, header=T)
dt_perf = as.data.frame(t(dt_perf))

## We name the row 
rownames(dt_expression_bis)=sort(rownames(dt_perf))

## we use the dt_perf rownames because of mismatch but the order is the same
## can be checked with :
cellLines[cellLines!=sort(rownames(dt_perf))]
sort(rownames(dt_perf))[cellLines!=sort(rownames(dt_perf))]


# Plot

## Set element for plot
accu = dt_perf$accuracy_mean
names(accu)=rownames(dt_perf)
col_from_accu = accu
size_from_accu = accu
labels_from_accu = accu
names(col_from_accu) = names(size_from_accu) = names(labels_from_accu) = names(accu)

## Set colors 
colors_from_plot = as.vector(wes_palette("Zissou1", 9, type = "continuous"))
col_from_accu[accu>0.9]=">.9"
col_from_accu[accu<=0.9]="=<.9"
col_from_accu[accu<=0.85]="=<.85"
col_from_accu[accu<=0.8]="=<.8"
col_from_accu[accu<=0.75]="=<.75"
col_from_accu[accu<=0.7]="=<.7"
col_from_accu[accu<=0.65]="=<.65"
col_from_accu[accu<=0.6]="=<.6"
col_from_accu[accu<=0.55]="=<.55"
n_colors = length(unique(col_from_accu))

## Titre and axis labels
f <- list(
  family = "Courier New, monospace",
  size = 18,
  color = "#7f7f7f ")
titre <- paste("tSNE CV mean accuracy\n",perf_filename)
x <- list(
  title = "tSNE dimension 1",
  titlefont = f)
y <- list(
  title = "tSNE dimension 2",
  titlefont = f)


# Produce the plot
j <- plot_ly(x=tsne$Y[,1], y=tsne$Y[,2],
             text=paste(rownames(dt_expression_bis),
                        round(accu[rownames(dt_expression_bis)],
                              5),
                        'accuracy'),
             color = col_from_accu[rownames(dt_expression_bis)],
             colors = colors_from_plot[1:n_colors],
             size = round((accu[rownames(dt_expression_bis)]+0.5)**20) ) %>% 
  layout(title = titre, xaxis = x, yaxis= y)
j


ggsave(filename = "tsnePerf.png", plot = j, dpi = 300,
  width = 10,
  height = 6)
out_file = "tmp"
save_html(j, paste(out_file,".html",collapse = ""), libdir = out_file )

