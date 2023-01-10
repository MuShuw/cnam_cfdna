import sys
import numpy as np
import pandas as pd
import argparse
import glob


def X_feat_selection(X, feat):
    print(FEAT_DIC[feat])
    return X[FEAT_DIC[feat]]

def set_generic_group(ranks, leuko_table, ubi_table):
    leuko_feats = leuko_table[-ranks:][['geneID','fold','mann_whitney_pval','chrom']]
    leuko_feats["CAT_geneClass"]=leuko_table["geneClass"]
    leuko_feats["geneGrouping"]="last_{}_leuko".format(ranks)
    leuko_feats["ontology"]="generic_leuko"
    ubi_sample = ubi_table.sample(ranks, random_state=42)[['geneID','chrom','CAT_geneClass']]
    ubi_sample["mann_whitney_pval"]=None
    ubi_sample['fold']=None
    ubi_sample["geneGrouping"]="{}_random_ubi".format(ranks)
    ubi_sample["ontology"]="generic_ubi_sample"
    return pd.concat([leuko_feats, ubi_sample])

def set_onthologic_group(ranks, feat_table, ontho):
    top_val = feat_table[:ranks][['geneID','fold','mann_whitney_pval','chrom','CAT_geneClass']]
    top_val["geneGrouping"]="first_{}_values".format(ranks)
    bot_val = feat_table[-ranks][['geneID','fold','mann_whitney_pval','chrom','CAT_geneClass']]
    bot_val["geneGrouping"]="last_{}_values".format(ranks)
    mid_val = feat_table[round((feat_table.shape[0]-ranks)/2):round((feat_table.shape[0]+ranks)/2)][['geneID','fold','mann_whitney_pval','chrom','CAT_geneClass']]
    mid_val["geneGrouping"]="middle_{}_values".format(ranks)
    pos = feat_table[feat_table.fold>=0][:ranks][['geneID','fold','mann_whitney_pval','chrom','CAT_geneClass']]
    pos["geneGrouping"]="first_{}_values_with_pos_FC".format(ranks)
    neg = feat_table[feat_table.fold<0][:ranks][['geneID','fold','mann_whitney_pval','chrom','CAT_geneClass']]
    neg["geneGrouping"]="first_{}_values_with_neg_FC".format(ranks)
    top_val["ontology"]=bot_val["ontology"]=mid_val["ontology"]=pos_val["ontology"]=neg_val["ontology"]=ontho
    return pd.concat([top_val, bot_val, pos, neg, mid_val])



parser = argparse.ArgumentParser(description='Classification score')
parser.add_argument('-r','--rank', type=int, default=100, help='The rank of top/bot gene to use', nargs='+')
parser.add_argument('-i','--indir', type=str, help='Input dir for onthologies tables')
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
    print("Unproper table for features leuko")
    sys.exit()

# Load feat_table for ubi
try :
    ubi_table = pd.read_csv(args.ubigenes, sep='\t', index_col=0)
except :
    print("Unproper table for features ubis")
    sys.exit()

Y=set_generic_group(args.rank[0], leuko_table, ubi_table)
for ontology in glob.glob(args.indir+'*'):
    try :
        if args.verbose:
            print('Trying to open '+ontology+' table')
        rank_ontho = pd.read_csv(ontology, sep='\t')
        if args.verbose:
            print('Opened '+ontology+' table')
    except :
        if args.verbose:
            print("Couldn't open "+ontology)
        pass
    
    filename = ontology.split('/')[-1].split('.')[0]

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
    for ranks in args.rank:
        X = set_onthologic_group(ranks, feat_table, ontology)
        Y = pd.concat([Y,X])
            # Creation du SVC

Y.to_csv('all_groups_of_size_'+str(args.rank[0])+'_'+args.outfile, sep='\t')
