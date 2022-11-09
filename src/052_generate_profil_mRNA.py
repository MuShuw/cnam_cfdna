import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

TARGET_ONTO="../Ontho_full/CL_0000771_eosinophil.tsv.gz"
LB="../Ontho_full/CL_0000738_leukocyte.tsv.gz"
RENAME = {
          'first_100_values':'top',
          'first_100_values_with_neg_FC':'topNF',
          'first_100_values_with_pos_FC':'topPF',
          'last_100_values':'bot',
          'middle_100_values':'mid'
          }

dt_group_fold = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_coding_mRNA_sortedBy_fold_groupings", sep="\t", index_col=0)
dt_group_pval = pd.read_csv("data/91_groupsForPubli/all_groups_of_size_100_ARN_coding_mRNA_sortedBy_mann_whitney_pval_groupings", sep="\t", index_col=0)

lb = dt_group_pval[(dt_group_pval["onthology"]==LB) & (dt_group_pval["geneGrouping"]=="last_100_values")]


dt_group_fold=dt_group_fold[dt_group_fold["onthology"]==TARGET_ONTO].reset_index(drop=True)
dt_group_pval=dt_group_pval[dt_group_pval["onthology"]==TARGET_ONTO].reset_index(drop=True)

dt_profil = pd.read_csv("data/02_tss_score/FD_coding_mRNAscoreTable.tsv", sep="\t")


dt_pval=pd.merge(dt_profil, dt_group_pval, how="right")
dt_fold=pd.merge(dt_profil, dt_group_fold, how="right")
lb=pd.merge(dt_profil, lb, how="right")

lb["geneGrouping"]="lbot"
dt_pval=pd.concat([dt_pval, lb])
dt_fold=pd.concat([dt_fold, lb])
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


