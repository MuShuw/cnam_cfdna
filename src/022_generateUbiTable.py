import numpy as np
import pandas as pd
UBIS="data/90_fantom/intermediary/HonEtAl_SuppTable4.csv"
INFO="data/90_fantom/intermediary/FANTOM_CAT.lv3_robust.info_table.gene.tsv.gz"

fant=pd.read_csv(UBIS, sep=";")
info=pd.read_csv(INFO, sep="\t")

fant = fant.rename({"DPI peak":"DPI"}, axis=1)
info = info.rename({"strongest_DPIClstrID":"DPI"}, axis=1)

dt = pd.merge(info, fant)

dt = dt[dt["geneClass"].isin(["coding_mRNA","lncRNA_intergenic"])] 

def mySplit(ser):
    return ser.split(':')[0]

X=dt["loc"].map(mySplit)

dt["chrom"]=X.copy()

dt=dt[~dt["chrom"].isin(["chrX","chrY"])]

# save that as you wish
