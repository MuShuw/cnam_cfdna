import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from bioinfokit import visuz

dt = pd.read_csv("data/90_fantom/onto/CL_0000775_neutrophil.tsv.gz", sep="\t")
#dt = dt[dt["fold"]!=0]
dt = dt.reset_index(drop=True)
info = pd.read_csv("data/90_fantom/intermediary/FANTOM_CAT.lv3_robust.info_table.gene.tsv.gz", sep="\t")
info = info[["geneID","geneClass"]]
df = pd.merge(dt, info, on="geneID", how="right")
#dt["log2FC"] = np.log2(dt["fold"].values)
visuz.GeneExpression.volcano(df=dt, lfc='log2FC', pv='mann_whitney_pval')

visuz.GeneExpression.volcano(df=df[df["geneClass"]=="lncRNA_intergenic"], lfc='fold', pv='mann_whitney_pval', show=True)
visuz.GeneExpression.volcano(df=df[df["geneClass"]=="coding_mRNA"], lfc='fold', pv='mann_whitney_pval', show=True)
