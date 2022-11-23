import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

TARGET_ONTO="../Ontho_full/CL_0000775_neutrophil.tsv.gz"


UBI="generic_ubi_sample"

LB="generic_leuko"

RENAME = {
          'first_100_values':'top',
          'first_100_values_with_neg_FC':'NF',
          'first_100_values_with_pos_FC':'PF',
          'last_100_values':'b',
          'middle_100_values':'mid'
          }

rna_type=["mrna","linc"]
srt_type=["fold","pval"]
rna_name=["coding_mRNA","lncRNA_intergenic"]
srt_name=["fold","mann_whitney_pval"]
rna_titl=["Coding","LincRNAs"]
srt_titl=["FC","P"]
rna_dic={rna_type[i]:rna_name[i] for i in range(len(rna_type))}
srt_dic={srt_type[i]:srt_name[i] for i in range(len(srt_type))}


dt_group={}
for i in rna_type:
    dt_group[i]={}
    for j in srt_type:
        dt_group[i][j]=pd.read_csv(f"data/91_groupsForPubli/all_groups_of_size_100_ARN_{rna_dic[i]}_sortedBy_{srt_dic[j]}_groupings", sep="\t", index_col=0)

lb = dt_group[i][j][(dt_group[i][j]["onthology"]==LB)]
ubi = dt_group[i][j][(dt_group[i][j]["onthology"]==UBI)]


for i in rna_type:
    for j in srt_type:
        dt_group[i][j]=dt_group[i][j][dt_group[i][j]["onthology"]==TARGET_ONTO].reset_index(drop=True)

dt_prof={}
dt_prof["linc"] = pd.read_csv("data/02_tss_score/FD_lncRNA_intergenicscoreTable.tsv", sep="\t")
dt_prof["mrna"] = pd.read_csv("data/02_tss_score/FD_coding_mRNAscoreTable.tsv", sep="\t")

dt_profil=pd.concat([dt_prof["linc"], dt_prof["mrna"]]).reset_index(drop=True)

lb=pd.merge(dt_profil, lb, how="right")
ubi=pd.merge(dt_profil, ubi, how="right")

lb["geneGrouping"]="leuko bottom"
ubi["geneGrouping"]="u"

dt={}
for i in rna_type:
    dt[i]={}
    for j in srt_type:
        dt[i][j]=pd.concat([pd.merge(dt_profil, dt_group[i][j], how="right"), lb, ubi])

##
for i in rna_type:
    for j in srt_type:
        dt[i][j]=dt[i][j].drop(['chrom', 'strand', 'start', 'CAT_geneClass', 'fold', 'mann_whitney_pval', 'onthology'], axis=1)
        dt[i][j]=dt[i][j][dt[i][j]["geneGrouping"]!="first_100_values"].reset_index(drop=True)

for i in rna_type:
    for j in srt_type:
        dt[i][j]=dt[i][j].melt(id_vars=["geneID","geneGrouping"], var_name='position', value_name='FD')
        dt[i][j]["sorting"]=str(srt_titl[srt_type.index(j)])+\
        " & "+\
        str(rna_titl[rna_type.index(i)])
        
DT=pd.concat([dt["mrna"]["pval"], dt["mrna"]["fold"], dt["linc"]["pval"], dt["linc"]["fold"]])

for key in RENAME:
    DT.loc[DT["geneGrouping"]==key,"geneGrouping"]=RENAME[key]

DT["position"]=DT["position"].astype("int") 
DT["geneGrouping"] = pd.Categorical(DT['geneGrouping'], ["PF", "NF", "mid", "u", "b"])
DT["sorting"] = pd.Categorical(DT['sorting'], ['P & Coding', 'FC & Coding', 'P & LincRNAs', 'FC & LincRNAs'])
DT=DT.sort_values(["position","geneGrouping","sorting"], ascending=[True,True,True])

#df=DT[(DT["position"]>-50) &
#      (DT["position"]<50) &
#      ~(DT["geneGrouping"].isna())
#     ].reset_index(drop=True)

df=DT[~(DT["geneGrouping"].isna())
     ].reset_index(drop=True)

    
sns.set_theme(font_scale=1.2, style="whitegrid")
with sns.color_palette("husl") :#and sns.plotting_context("paper",font_scale=2):
    g = sns.FacetGrid(
                  df,
                  col="sorting",
                  row="geneGrouping",
                  hue="geneGrouping",
                  margin_titles=True,
                  gridspec_kws={"wspace":0.05, "hspace":0.05}
                )
    g.map_dataframe(
                    sns.lineplot,
                    x="position",
                    y="FD",
                    n_boot=100,
                    #err_style="bars",
                    #err_kws={"errorevery":500},
                    #ci="sd",
                    #estimator=np.mean
                    )
             
for ax in g.axes.flat:
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2500)) # normal value 2500
    plt.setp(ax.texts, text="")
    ax.tick_params(axis='x', rotation=45)

g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
#g.add_legend()
#plt.show()
plt.savefig("neutrophil_profil.png")
