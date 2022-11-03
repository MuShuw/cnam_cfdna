#! /usr/bin/python3
# This one will use the preset feature tables
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
from sklearn import svm


PERFORMANCES=['accuracy','f1','precision','recall']
CV = ShuffleSplit(n_splits=100, test_size=0.3, random_state=42) 
GENE_GROUPS=['PFb','NFb','PFNF','PFu','NFu','bu','PFlb','NFlb','blb','PFmid','NFmid','bmid','bNF']
GENE_TYPE_SEL=['sortcb','cbsort']
FEAT_SELECTIONS=['Ulz',
                 'Ulzs',
                 'Ulzm-var',
                 'Ulzm-vars',
                 'Ulzm-peak',
                 'Ulzm-peaks',
                 'Ulzm-var-peak',
                 'Ulzm-vars-peak',
                 'Ulzm-var-peaks',
                 'Ulzm-vars-peaks',
                 'Full',
                 'final',
                 'all']
FEATS_NAMES = ['FDmean', 
               'TSS2kFD', 'TSS2kFD_scaled',
               'TSSFD', 'TSSFD_scaled',
               'TSSvar', 'TSSvar_scaled',
               'PseudoVarFD', 'PseudoVarFD_scaled',
               'pseudoVarPlus1', 'pseudoVarPlus1_scaled',
               'TSS2kRATIO'] # To edit accordingly to the feat table creation
FEAT_DIC={
    FEAT_SELECTIONS[0]:[FEATS_NAMES[i] for i in [1,3]],
    FEAT_SELECTIONS[1]:[FEATS_NAMES[i] for i in [2,4]],
    FEAT_SELECTIONS[2]:[FEATS_NAMES[i] for i in [1,4,5]],
    FEAT_SELECTIONS[3]:[FEATS_NAMES[i] for i in [1,4,6]],
    FEAT_SELECTIONS[4]:[FEATS_NAMES[i] for i in [1,4,9]],
    FEAT_SELECTIONS[5]:[FEATS_NAMES[i] for i in [1,4,10]],
    FEAT_SELECTIONS[6]:[FEATS_NAMES[i] for i in [1,4,5,9]],
    FEAT_SELECTIONS[7]:[FEATS_NAMES[i] for i in [1,4,6,9]],
    FEAT_SELECTIONS[8]:[FEATS_NAMES[i] for i in [1,4,5,10]],
    FEAT_SELECTIONS[9]:[FEATS_NAMES[i] for i in [1,4,6,10]],
    FEAT_SELECTIONS[10]:FEATS_NAMES,# to edit once fd ratio has been added
    FEAT_SELECTIONS[11]:[FEATS_NAMES[i] for i in [0,11,4,6,8,10]]
}

def set_generic_group(ranks, leuko_table, ubi_table):
    leuko_feats = leuko_table[-ranks:]
    leuko_feats["CAT_geneClass"]=leuko_table["geneClass"]
    ubi_sample = ubi_table.sample(ranks, random_state=42)
    ubi_sample["mann_whitney_pval"]=None
    ubi_sample['fold']=None
    return leuko_feats, ubi_sample

def set_onthologic_group(ranks, feat_table):
    top_val = feat_table[:ranks]
    bot_val = feat_table[-ranks:]
    mid_val = feat_table[round((feat_table.shape[0]-ranks)/2):round((feat_table.shape[0]+ranks)/2)]
    pos = feat_table[feat_table.fold>=0][:ranks]
    neg = feat_table[feat_table.fold<0][:ranks]
#    bot_pos = feat_table[feat_table.fold>=0][-ranks:]
#    bot_neg = feat_table[feat_table.fold<0][-ranks:]
    return top_val, bot_val, pos, neg, mid_val

def X_creation(group, pos, neg, top_val, bot_val, mid_val, leuko_bot, ubi_sample):
    if group == 'PFb':
        X = pd.concat([pos,bot_val])
    if group == 'NFb':
        X = pd.concat([neg,bot_val])
    if group == 'PFNF':
        X = pd.concat([pos,neg])
    if group == 'PFu':
        X = pd.concat([pos,ubi_sample])
    if group == 'NFu':
        X = pd.concat([neg,ubi_sample])
    if group == 'bu':
        X = pd.concat([bot_val,ubi_sample])
    if group == 'PFlb':
        X = pd.concat([pos,leuko_bot])
    if group == 'NFlb':
        X = pd.concat([neg,leuko_bot])
    if group == 'blb':
        X = pd.concat([bot_val,leuko_bot])
    if group == 'PFmid':
        X = pd.concat([pos,mid_val])
    if group == 'NFmid':
        X = pd.concat([neg,mid_val])
    if group == 'bmid':
        X = pd.concat([bot_val,mid_val])
    if group == 'bNF':
        X = pd.concat([bot_val,neg])
        
    return X

def X_feat_selection(X, feat):
    print(FEAT_DIC[feat])
    return X[FEAT_DIC[feat]]

parser = argparse.ArgumentParser(description='Classification score')
parser.add_argument('-r','--rank', type=int, default=100, help='The rank of top/bot gene to use', nargs='+')
parser.add_argument('-i','--indir', type=str, help='Input dir for onthologies tables')
parser.add_argument('-g','--grouping', type=str, default='PFb', help='The groups of genes to compare', nargs='+')
parser.add_argument('-f','--feats', type=str, default='Ulz', help='Choose which features to use', nargs='+')
parser.add_argument('-o','--outfile', type=str, help='Output file for the scores')
parser.add_argument('-if','--infeats', type=str, help='Input table of features')
parser.add_argument('-lbf','--lbinfeats', type=str, help='Input table of features sorted for leuko')
parser.add_argument('-u','--ubigenes', type=str, help='Input table of ubiquitious gene')
parser.add_argument('-sb','--sortby', type=str, default='mann_whitney_pval', help="Sort by 'fold' or 'mann_whitney_pval'(==default)")
parser.add_argument('-v','--verbose', action='store_true', help='Boolean (default=True) Whether to print information about processing')
args = parser.parse_args()

# We set a constant concerning the sorting
if args.sortby not in ['fold','mann_whitney_pval','both']:
    args.sortby='mann_whitney_pval'
    if args.verbose:
        print('Incorrect --sortby argument, sort set by default(mann_whitney_pval)')
if args.sortby=='fold':
    SORT='F'
elif args.sortby=='mann_whitney_pval':
    SORT="MW"
else:
    SORT='B'

# Load feat_table
try :
    feat_table = pd.read_csv(args.infeats, sep='\t', index_col=0)
except :
    print("Unproper table for features")
    sys.exit()

# Load feat_table sorted for leuko 
try :
    leuko_table = pd.read_csv(args.lbinfeats, sep='\t', index_col=0)
except :
    print("Unproper table for features ubis")
    sys.exit()

# Load feat_table for ubi
try :
    ubi_table = pd.read_csv(args.ubigenes, sep='\t', index_col=0)
except :
    print("Unproper table for features")
    sys.exit()

# Check if grouping are correct
classif_grouping=[]
for group in args.grouping:
    if group in GENE_GROUPS:
        classif_grouping.append(group)

# Check if feats selections are correct
feats_sel=[]
if 'all' in args.feats:
    feats_sel = FEAT_SELECTIONS[:-1]
else :
    for feat in args.feats:
        if feat in FEAT_SELECTIONS:
            feats_sel.append(feat)

# Building dictionnaries
perf_per_ontho={}
for group in classif_grouping:
    perf_per_ontho[group]={}
    for rank in args.rank:
        perf_per_ontho[group][rank]={}
        for feat in feats_sel:
            perf_per_ontho[group][rank][feat]={}


for onthology in glob.glob(args.indir+'*'):
    try :
        if args.verbose:
            print('Trying to open '+onthology+' table')
        rank_ontho = pd.read_csv(onthology, sep='\t')
        if args.verbose:
            print('Opened '+onthology+' table')
    except :
        if args.verbose:
            print("Couldn't open "+onthology)
        pass
    
    filename = onthology.split('/')[-1].split('.')[0]

    print('Processing '+filename)
    # Preparation du fd
    ## First drop previous ranks
    try:
        feat_table = feat_table.drop(['fold','mann_whitney_pval'], axis=1)
    except:
        if args.verbose:
            print('No fold')
    ## Add new ranks
    feat_table = feat_table.merge(rank_ontho[['fold','geneID','mann_whitney_pval']], on='geneID', how='left')
    ## Sort with new rank
    if SORT!='B':
        if SORT=='MW':
            feat_table = feat_table.sort_values(args.sortby, ascending=True)
        elif SORT=='F':
            feat_table = feat_table.sort_values(args.sortby, ascending=False)
    for group in classif_grouping:
        for ranks in args.rank:
            # Label creation
            y = pd.concat([ pd.Series(np.repeat(1,ranks)), pd.Series(np.repeat(0,ranks))])
            # Selection des features a implementer
            top_val, bot_val, pos, neg, mid_val = set_onthologic_group(ranks, feat_table)
            leuko_feats, ubi_sample = set_generic_group(ranks, leuko_table, ubi_table)
            # Concat de X
            X = X_creation(group, pos, neg, top_val, bot_val, leuko_feats, ubi_sample, mid_val)
            for feat in feats_sel:
                perf_per_ontho[group][ranks][feat][filename]={}
                Xb = X_feat_selection(X, feat)
                # Creation du SVC
                clf = make_pipeline(preprocessing.StandardScaler(), svm.SVC(kernel='rbf',C=1))

                # Enregistrement des scores
                for perf in PERFORMANCES:
                    scores = cross_val_score(clf, Xb, y, cv=CV, scoring=perf)
                    perf_per_ontho[group][ranks][feat][filename][perf+'_mean'] = scores.mean()
                    perf_per_ontho[group][ranks][feat][filename][perf+'_std'] = scores.std()

for group in perf_per_ontho.keys():
    for rank in perf_per_ontho[group].keys():
        for feat in perf_per_ontho[group][rank].keys():
            Scores=pd.DataFrame.from_dict(perf_per_ontho[group][rank][feat])
            print('lincRNA_rank'+str(rank)+'_'+args.outfile)
            Scores.to_csv(args.outfile+str(rank)+'.tsv', sep='\t')
