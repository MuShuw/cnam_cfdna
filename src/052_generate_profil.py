import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

TARGET_ONTO="../Ontho_full/CL_0000771_eosinophil.tsv.gz"

TARGET_ONTO="../Ontho_full/CL_0000775_neutrophil.tsv.gz"

TARGET_ONTO="../Ontho_full/CL_0000767_basophil.tsv.gz"

TARGET_ONTO="../Ontho_full/CL_0000094_granulocyte.tsv.gz"

TARGET_ONTO="../Ontho_full/CL_0000763_myeloid_cell.tsv.gz"

TARGET_ONTO="../Ontho_full/UBERON_0002193_hemolymphoid_system.tsv.gz"


UBI="generic_ubi_sample"
LB="generic_leuko"
RENAME = {
          'first_100_values':'top',
          'first_100_values_with_neg_FC':'topNF',
          'first_100_values_with_pos_FC':'topPF',
          'last_100_values':'bot',
          'middle_100_values':'mid'
          }

dt_group_fold = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_lncRNA_intergenic_sortedBy_fold_groupings", sep="\t", index_col=0)
dt_group_pval = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_lncRNA_intergenic_sortedBy_mann_whitney_pval_groupings", sep="\t", index_col=0)


dt_group_fold = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_coding_mRNA_sortedBy_fold_groupings", sep="\t", index_col=0)
dt_group_pval = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_coding_mRNA_sortedBy_mann_whitney_pval_groupings", sep="\t", index_col=0)


lb = dt_group_pval[(dt_group_pval["onthology"]==LB)]
ubi = dt_group_pval[(dt_group_pval["onthology"]==UBI)]

dt_group_fold=dt_group_fold[dt_group_fold["onthology"]==TARGET_ONTO].reset_index(drop=True)
dt_group_pval=dt_group_pval[dt_group_pval["onthology"]==TARGET_ONTO].reset_index(drop=True)

dt_profil_linc = pd.read_csv("data/02_tss_score/FD_lncRNA_intergenicscoreTable.tsv", sep="\t")
dt_profil_mrna = pd.read_csv("data/02_tss_score/FD_coding_mRNAscoreTable.tsv", sep="\t")
dt_profil=pd.concat([dt_profil_linc,dt_profil_mrna]).reset_index(drop=True)

dt_pval=pd.merge(dt_profil, dt_group_pval, how="right")
dt_fold=pd.merge(dt_profil, dt_group_fold, how="right")
lb=pd.merge(dt_profil, lb, how="right")
ubi=pd.merge(dt_profil, ubi, how="right")

lb["geneGrouping"]="lbot"
ubi["geneGrouping"]="ubi"

dt_pval=pd.concat([dt_pval, lb, ubi])
dt_fold=pd.concat([dt_fold, lb, ubi])
dt_pval=dt_pval.drop(['chrom', 'strand', 'start', 'CAT_geneClass', 'fold', 'mann_whitney_pval', 'onthology'], axis=1)
dt_fold=dt_fold.drop(['chrom', 'strand', 'start', 'CAT_geneClass', 'fold', 'mann_whitney_pval', 'onthology'], axis=1)

dt_fold=dt_fold[dt_fold["geneGrouping"]!="first_100_values"].reset_index(drop=True)
dt_pval=dt_pval[dt_pval["geneGrouping"]!="first_100_values"].reset_index(drop=True)

dt_pval=dt_pval.melt(id_vars=["geneID","geneGrouping"], var_name='position', value_name='FD')
dt_fold=dt_fold.melt(id_vars=["geneID","geneGrouping"], var_name='position', value_name='FD')
dt_pval["sorting"]="pval"
dt_fold["sorting"]="fold"

dt=pd.concat([dt_pval, dt_fold])

for key in RENAME:
    dt.loc[dt["geneGrouping"]==key,"geneGrouping"]=RENAME[key]

# plot 1
ax=sns.lineplot(data=dt_pval, x="position", y="FD", hue="geneGrouping", ci=None)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
plt.show()

# plot 2
g = sns.FacetGrid(dt_pval, row="geneGrouping", hue="geneGrouping")
g.map_dataframe(sns.lineplot,x="position",y="FD",ci=None)
for ax in g.axes.flat:
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
    ax.set_title("")
g.add_legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.show()

# plot 3
g = sns.FacetGrid(
                  dt,
                  col="sorting",
                  row="geneGrouping",
                  hue="geneGrouping",
                  margin_titles=True
                )

g.map_dataframe(
                sns.lineplot,
                x="position",
                y="FD",
                ci=None,
                estimator=np.mean
                )
for ax in g.axes.flat:
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2500))
    plt.setp(ax.texts, text="")

g.set_titles(row_template = '', col_template = '{col_name}')
g.add_legend()
plt.show()

