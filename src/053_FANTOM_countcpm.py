#! usr/bin/python3

import statistics
import numpy as np
import pandas as pd

# First get the read count for FANTOMv5
df = pd.read_csv('data/expression/FANTOM/FANTOM_CAT.expression_atlas.gene.lv3_robust.count.tsv.gz', sep='\t', index_col=0)

sample_depth = dt.apply(sum, axis=0)
factor_ratio = sample_depth/1_000_000

tpm_df = df/factor_ratio

tpm_df.to_csv('results/FANTOM_CAT.expression_atlas.gene.lv3_robust.tpm_ELYAS.tsv', sep='\t')

