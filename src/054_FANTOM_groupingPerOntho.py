#! /usr/bin/python3

import statistics
import re
import numpy as np
import pandas as pd

# To be modified to incorporate args for table
def get_mean_for_ontho(pd_serie):
    tmp_table = tpm_sample[pd_serie[0].split(',')]
    tmp_mean = tmp_table.apply(statistics.median, axis=1)
    print(pd_serie[1])
    return pd.Series(tmp_mean, name = pd_serie[1])

MYREGEX = re.compile('[ ,-]')
onthology_sample_association = pd.read_csv('data/expression/FANTOM/supp/supp_table_10.sample_ontology_information.csv', sep='\t', header=0)
real_ontho = onthology_sample_association.apply( lambda x : '_'.join([*x[0].split(':'), *MYREGEX.split(x[1])]), axis=1)
onthology_sample_association=onthology_sample_association.assign(real_ontho=real_ontho)
onthology_sample_association=onthology_sample_association[['CAGE_lib_ID','real_ontho']]   

tpm_sample = pd.read_csv('results/FANTOM_CAT.expression_atlas.gene.lv3_robust.tpm_ELYAS.tsv',sep='\t', index_col=0)

x_tmp = onthology_sample_association.apply( get_mean_for_ontho, axis=1)
x_tmp.index=onthology_sample_association.real_ontho
x_tmp=np.transpose(x_tmp)

x_tmp.to_csv('results/FANTOM_CAT.expression_atlas.gene.lv3_robust.MEDIAN_tpm_PER_ONTHO_ELYAS.tsv', sep='\t')

