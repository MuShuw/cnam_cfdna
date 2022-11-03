#! /usr/bin/python3
import sys
import numpy as np
import pandas as pd
import statistics
import argparse
import glob
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.preprocessing import scale
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import ShuffleSplit
from sklearn.pipeline import make_pipeline

def reduce_mem_usage(props):
    start_mem_usg = props.memory_usage() / 1024**2 
    print("Memory usage of properties dataframe is :",start_mem_usg," MB")
    if props.dtype != object:  # Exclude strings

        # Print current column type
        print("******************************")
        print("Column: ",props.name)
        print("dtype before: ",props.dtype)

        # make variables for Int, max and min
        IsInt = False
        mx = props.max()
        mn = props.min()

        # test if column can be converted to an integer
        asint = props.fillna(0).astype(np.int64)
        result = (props - asint)
        result = result.sum()
        if result > -0.01 and result < 0.01:
            IsInt = True


        # Make Integer/unsigned Integer datatypes
        if IsInt:
            if mn >= 0:
                if mx < 255:
                    props = props.astype(np.uint8)
                elif mx < 65535:
                    props = props.astype(np.uint16)
                elif mx < 4294967295:
                    props = props.astype(np.uint32)
                else:
                    props = props.astype(np.uint64)
            else:
                if mn > np.iinfo(np.int8).min and mx < np.iinfo(np.int8).max:
                    props = props.astype(np.int8)
                elif mn > np.iinfo(np.int16).min and mx < np.iinfo(np.int16).max:
                    props = props.astype(np.int16)
                elif mn > np.iinfo(np.int32).min and mx < np.iinfo(np.int32).max:
                    props = props.astype(np.int32)
                elif mn > np.iinfo(np.int64).min and mx < np.iinfo(np.int64).max:
                    props = props.astype(np.int64)    

        # Make float datatypes 32 bit
        else:
            props = props.astype(np.float32)

        # Print new column type
        print("dtype after: ",props.dtype)
        print("******************************")
    
    # Print final result
    print("___MEMORY USAGE AFTER COMPLETION:___")
    mem_usg = props.memory_usage() / 1024**2 
    print("Memory usage is: ",mem_usg," MB")
    print("This is ",100*mem_usg/start_mem_usg,"% of the initial size")
    return props

TSS = 5000
PERFORMANCES=['accuracy','f1','precision','recall']
CV = ShuffleSplit(n_splits=100, test_size=0.3, random_state=42) 

parser = argparse.ArgumentParser(description='Feature table creation')
parser.add_argument('-o','--outfile', type=str, help='Output file for the scores')
parser.add_argument('-iw','--inwps', type=str, help='Input table of FD scores')
parser.add_argument('-if','--infd', type=str, help='Input table of FD scores')
parser.add_argument('-v','--verbose', action='store_true', help='Boolean (default=True) Whether to print information about processing')

args = parser.parse_args()

# Get the FD score table #################################################################################################################
features=[]
try:
    print('Ouverture de '+args.infd)
    score_table = pd.read_csv(args.infd, sep='\t')
except:
    print('Unproper table for FD.')
    sys.exit()
# Extract FD features 
if args.verbose:
    print('Extractcting FD features\n')
FDmean = score_table.iloc[:,:10001].apply(statistics.mean, axis=1)
FDmean.name='FDmean'
features.append(FDmean)
TSS2kFD = score_table.iloc[:,(TSS-1000):(TSS+1000+1)].apply(statistics.mean, axis=1) 
TSS2kFD.name='TSS2kFD'
features.append(TSS2kFD)
TSSFD = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(statistics.mean, axis=1)
TSSFD.name='TSSFD'
features.append(TSSFD)
TSSvar = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(statistics.variance, axis=1)
TSSvar.name='TSSvar'
features.append(TSSvar)
PseudoVarFD = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(lambda x: max(x)-min(x), axis=1)
PseudoVarFD.name='PseudoVarFD'
features.append(PseudoVarFD)
if args.verbose:
    print('Done\n')

# Get the geneID, CAT_CAT_geneClass and chrom
ids = score_table.iloc[:,10001:]

# Merging feats with score_table
for feat in features: 
    score_table = score_table.join(feat)

# Create feat_table
feat_table = score_table.iloc[:,10001:]

# Prepare for scaling
score_table = score_table.iloc[:,0:10001].apply(lambda x : x.astype(float))

# Scaling FD #############################################################################################################################
if args.verbose:
    print('Scaling FD')
score_table = score_table.iloc[:,:10001].apply(scale, axis = 1 ).apply(pd.Series).apply(reduce_mem_usage, axis=0).join(ids)
if args.verbose:
    print('Done\n')
# Extract scaled FD features
features=[]
if args.verbose:
    print('Extractcting scaled FD features\n')
TSS2kFD_scaled = score_table.iloc[:,(TSS-1000):(TSS+1000+1)].apply(statistics.mean, axis=1) 
TSS2kFD_scaled.name='TSS2kFD_scaled'
features.append(TSS2kFD_scaled)
TSSFD_scaled = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(statistics.mean, axis=1)
TSSFD_scaled.name='TSSFD_scaled'
features.append(TSSFD_scaled)
TSSvar_scaled = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(statistics.variance, axis=1)
TSSvar_scaled.name='TSSvar_scaled'
features.append(TSSvar_scaled)
PseudoVarFD_scaled = score_table.iloc[:,(TSS-150):(TSS+50+1)].apply(lambda x: max(x)-min(x), axis=1)
PseudoVarFD_scaled.name='PseudoVarFD_scaled'
features.append(PseudoVarFD_scaled)
if args.verbose:
    print('Done\n')

# Get the geneID, CAT_CAT_geneClass and chrom
ids = score_table.iloc[:,10001:]

# Merging feats with score_table
for feat in features: 
    score_table = score_table.join(feat)

# append feat_table
feat_table = feat_table.merge(score_table.iloc[:,10001:], on=['geneID', 'chrom', 'strand', 'start', 'CAT_geneClass'], how='left')

# Get the WPS score table #################################################################################################################
features=[]
try:
    print('Ouverture de '+args.inwps)
    score_table = pd.read_csv(args.inwps, sep='\t')
except:
    print('Unproper table for WPS.')
    sys.exit()
# Extract WPS features 
if args.verbose:
    print('Extractcting WPS features\n')
pseudoVarPlus1 = score_table.iloc[:,TSS:(TSS+300+1)].apply(lambda x: max(x[0:(200+1)])-min(x[(150):(300+1)]), axis=1)
pseudoVarPlus1.name='pseudoVarPlus1'
if args.verbose:
    print('Done\n')

# Get the geneID, CAT_CAT_geneClass and chrom
ids = score_table.iloc[:,10001:]

# Merging feats with score_table
score_table=score_table.join(pseudoVarPlus1)

# append feat_table
feat_table = feat_table.merge(score_table.iloc[:,10001:], on=['geneID', 'chrom', 'strand', 'start', 'CAT_geneClass'], how='left')

# Prepare for scaling
score_table = score_table.iloc[:,0:10001].apply(lambda x : x.astype(float))

# Scaling WPS #############################################################################################################################
features=[]
if args.verbose:
    print('Scaling WPS')
score_table = score_table.iloc[:,:10001].apply(scale, axis = 1 ).apply(pd.Series).apply(reduce_mem_usage, axis=0).join(ids)
if args.verbose:
    print('Done\n')
# Extract scaled FD features
if args.verbose:
    print('Extractcting scaled FD features\n')
pseudoVarPlus1_scaled = score_table.iloc[:,TSS:(TSS+300+1)].apply(lambda x: max(x[0:(200+1)])-min(x[(150):(300+1)]), axis=1)
pseudoVarPlus1_scaled.name='pseudoVarPlus1_scaled'
if args.verbose:
    print('Done\n')

# Get the geneID, CAT_CAT_geneClass and chrom
ids = score_table.iloc[:,10001:]

# Merging feats with score_table
score_table=score_table.join(pseudoVarPlus1_scaled)

# append feat_table
feat_table = feat_table.merge(score_table.iloc[:,10001:], on=['geneID', 'chrom', 'strand', 'start', 'CAT_geneClass'], how='left')

# freeing memory ??
del(score_table)

# typing feature table
feat_table['TSS2kRATIO']=(feat_table.TSS2kFD/feat_table.FDmean)
feat_table = feat_table.drop(feat_table[feat_table.TSS2kRATIO.isnull()].index)
print(feat_table) 
print(feat_table.columns)
feat_table.iloc[:,5:] = feat_table.iloc[:,5:].apply(lambda x : x.astype(float))
feat_table.to_csv(args.outfile, sep='\t')

