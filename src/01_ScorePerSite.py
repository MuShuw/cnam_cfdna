#!/usr/bin/env python3

"""
This script is calculating the score for all TSS provided.

The scores per site are then saved in pickle format.
"""
import argparse
import copy
import glob
import pickle
import re
import statistics
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import pandas as pd


# define the name of the directory for temporary processing
PATH_TMP = "/tmp/table_creation"  # Put this into args

# define the position of the TSS
RANGE = 5000

# define the temporary files prefixes
PICKLE_PREFIX = "tmp_chr_"

# define possible gene class
GENE_CLASSES = ['lncRNA_intergenic', 'coding_mRNA']


def set_parser():
    """
    Use this function to build the script parser.

    Returns
    -------
    The parsed args.

    """
    parser = argparse.ArgumentParser(
        description='Sequence Score Table Creation')
    parser.add_argument('-p', '--pos',
                        type=int,
                        default=0,
                        help='The assumed position of TSSs \n 0 for top cap count\
                        nt in cluster \n 1 for start of cluster \n 2 for end of\
                        cluster \n else for mid of cluster')
    parser.add_argument('-i', '--indir',
                        type=str,
                        help='Input folder for sequences results')
    parser.add_argument('-o', '--outdir',
                        type=str,
                        help='Output folder for the scores tables')
    parser.add_argument('-ft5', '--fantom5',
                        type=str,
                        help='Input folder for the Fantom tables ')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Boolean (default=True) Whether to print information\
                        about processing')
    parser.add_argument('-gc', '--geneclass',
                        type=str,
                        default='lncRNA_intergenic',
                        help="Select the FANTOM gene class selected (default == \
                        lncRNA_intergenic)")
    args = parser.parse_args()
    return(args)


def extract_strand_start(pd_df_line, pos=0):
    """
    Cette fonction récupère le point le plus amont du transcript
    selon qu'il soit forward ou reverse.
    """
    splt = re.compile(":|\.\.|,")  # expression used to split the DPIcluster
    chrom, start, end, strand = splt.split(pd_df_line['strongest_DPIClstrID'])
    if strand == '-':
        tmp = int(start)
        start = int(end)
        end = tmp
    else:
        start = int(start)
        end = int(end)
    if pos == 0:
        TSS_start = pd_df_line.TSS
    elif pos == 1:
        TSS_start = start
    elif pos == 2:
        TSS_start = end
    else:
        TSS_start = round(statistics.mean([start, end]))
    to_return = pd.Series([pd_df_line.geneID, chrom[3:], strand, TSS_start],
                          index=['geneID', 'chrom', 'strand', 'start'])
    return(to_return)


def check_strand(pd_df_line):
    """
    Retourne un range autour du tss.
    """
    to_return = pd_df_line['geneID'], pd_df_line['chrom'], int(
        pd_df_line['start'])-RANGE, int(pd_df_line['start'])+RANGE+1, pd_df_line['strand']
    return(pd.Series(to_return, index=['geneID', 'chrom', 'first_pos', 'last_pos', 'strand']))


def sort_ranges(range_list):
    sorted_ranges = []
    tmp_range_series = pd.Series(range_list)
    tmp_range_list = pd.DataFrame(tmp_range_series, columns=['range']).assign(
        min_val=tmp_range_series.apply(min)).sort_values('min_val')['range']
    return(list(tmp_range_list))


def merge_range(list_range):
    merged = []
    for higher in list_range:
        if merged:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if min(higher) <= max(lower)+1:
                upper_bound = max(max(lower), max(higher))
                # replace by merged interval
                merged[-1] = range(min(min(higher), min(lower)), upper_bound+1)
            else:
                merged.append(copy.deepcopy(higher))
        else:
            merged.append(copy.deepcopy(higher))
    return merged


def clean_score_range(score_range):
    ref_chrom = [str(x) for x in list(range(1, 23))]
    return(score_range.loc[score_range.chrom.isin(ref_chrom)])


def get_usable_range_withoutX(score_range):
    # This is valable only if we sorted out the X, Y and Mt chroms
    # /!\ BUT DOES NOT WORK SINCE chrom columns is str type
    # chroms = list(sorted(score_range.chrom.unique().astype(int)))

    # Else we use
    chroms = list(sorted(score_range.chrom.unique()))
    # First we list the chromosoms concerned
    plop = dict()
    for chrom in chroms:
        plop[chrom] = list((score_range[score_range.chrom == chrom]['geneID'].values,
                            score_range[score_range.chrom ==
                                        chrom]['first_pos'].values,
                            score_range[score_range.chrom ==
                                        chrom]['last_pos'].values,
                            score_range[score_range.chrom == chrom]['strand'].values))
    # We create the final iterator
    iterators = dict()
    for key in plop.keys():
        iterators[key] = list(
            zip(plop[key][0], plop[key][1], plop[key][2], plop[key][3]))

    # Then we return the iterator
    return(iterators)


def extract_for_chr(chrom, iterator):
    try:
        chrom_files = glob.glob(PATH_IN+'*_'+chrom+'.npz')
        if len(chrom_files) != 1:
            print("Wrong number of files ending with chr"+chrom +
                  ".npz \n Skipping chromosome number"+chrom)
            return()
        file_name = chrom_files[0]
        file_chr = pd.DataFrame(np.load(file_name)["arr_0"])
        print(chrom, ' ', 'fichier : ', file_name)
        final_range_chr = iterator
        WPS = []
        FD = []
        for eTSS in final_range_chr:
            if eTSS[3] == '+':
                wps = list(file_chr.wps[eTSS[1]:eTSS[2]])
                fd = list(file_chr.fd[eTSS[1]:eTSS[2]])
            else:
                wps = list(reversed(list(file_chr.wps[eTSS[1]:eTSS[2]])))
                fd = list(reversed(list(file_chr.fd[eTSS[1]:eTSS[2]])))
            wps.append(eTSS[0])
            fd.append(eTSS[0])
            WPS.append(wps)
            FD.append(fd)
        WPS_table = pd.DataFrame(WPS)
        FD_table = pd.DataFrame(FD)
        cols = list(range(-5000, 5001))
        cols.append("geneID")
        WPS_table.columns = FD_table.columns = cols
        print('Chromosome '+chrom+' done.\n')
        return WPS_table, FD_table
    except IOError:
        print('fail')
        return()
    return()


def extract_WPS_FD_pos(score_range):
    iterator = get_usable_range_withoutX(score_range)
    p = Pool(8)
    p.starmap(extract_process, zip(iterator.keys(), iterator.values()))
    p.close()
    p.join()
    return


def extract_process(chrom, iterator):
    tmp_score = extract_for_chr(chrom, iterator)
    print('chr'+chrom+' done\n')
    pickle.dump(tmp_score, open(PICKLE_PREFIX+chrom+'.pickle', "wb"))
    return


def save_table(table, chroms):
    fd_saved_score = dict()
    wps_saved_score = dict()
    for chrom in chroms:
        wps_saved_score[chrom], fd_saved_score[chrom] = \
            pickle.load(open(PICKLE_PREFIX+chrom+'.pickle', "rb"))

    wps_saved_score = pd.concat(wps_saved_score.values())
    fd_saved_score = pd.concat(fd_saved_score.values())

    wps_saved_score = wps_saved_score.merge(
        table, how='left', left_index=True, on='geneID'
    )
    fd_saved_score = fd_saved_score.merge(
        table, how='left', left_index=True, on='geneID'
    )

    wps_saved_score.to_csv(PATH_OUT+'/WPS_'+GENE_CLASS +
                           'scoreTable.tsv', sep='\t', index=False)
    fd_saved_score.to_csv(PATH_OUT+'/FD_'+GENE_CLASS +
                          'scoreTable.tsv', sep='\t', index=False)

    return


if __name__ == '__main__':
    # Check if in dir in exists:
    args = set_parser()
    p = Path(args.indir)
    if not p.exists() or not p.is_dir():
        print("The input folder can't be found. Correct that.")
        sys.exit()
    else:
        PATH_IN = args.indir

    # Check if given gene class is correct
    if not (args.geneclass in GENE_CLASSES):
        print("The given gene class is incorrect. Coorect that.")
        sys.exit()
    else:
        GENE_CLASS = args.geneclass

    # Check if output file does exist
    p = Path(args.outdir)
    if p.exists():
        print("The output file name is already taken. Correct that.")
        # sys.exit() check a proper way to handle that
        PATH_OUT = args.outdir
    else:
        PATH_OUT = args.outdir

    # Check if fantom folder does exist
    p = Path(args.indir)
    if not p.exists() or not p.is_dir():
        print("The fantom folder can't be found. Correct that.")
        sys.exit()
    else:
        PATH_FT = args.fantom5
    # Creating variable for fantom files
    AD_INFO_TABLE = PATH_FT+'FANTOM_CAT.lv3_robust.info_table.gene.tsv.gz'
    AD_ASSO_TABLE = PATH_FT+'supp_table_11.cell_type_gene_association.tsv'
    AD_CAGE_BED = PATH_FT+'FANTOM_CAT.lv3_robust.CAGE_cluster.bed.gz'

    # Loading fantom5 tables
    # Add try blocks
    robust_major_tss = pd.read_csv(AD_INFO_TABLE, sep='\t')
    cell_type_gene_association = pd.read_csv(AD_ASSO_TABLE, sep='\t')
    cell_type_gene_association.columns = ['geneID',
                                          'CAT_geneName',
                                          'CAT_geneCategory',
                                          'CAT_geneClass',
                                          'top_associated_sample_ontology',
                                          'num_associated_sample_ontology',
                                          'CAT_browser_link',
                                          'associated_sample_ontology']
    bed_point = pd.read_csv(AD_CAGE_BED, sep='\t', header=None)

    bed_point = bed_point[[3, 7]]
    bed_point.columns = ['strongest_DPIClstrID', 'TSS']
    robust_major_tss = robust_major_tss.merge(
        bed_point, on='strongest_DPIClstrID', how='left')
    # We filter the biotype
    geneTypeSelection = cell_type_gene_association[
        cell_type_gene_association.CAT_geneClass == GENE_CLASS]
    # We join the filtered table with the TSSs table (still not sure why
    # I didn't filtered directly the TSS table)
    detailed_table = geneTypeSelection.merge(
        robust_major_tss, how='left', left_on='geneID', right_on='geneID')

    # Let's get thoose TSS
    simple_table = detailed_table.apply(extract_strand_start, axis=1)
    # The next step is to filter some chromosomes
    simple_table = simple_table.loc[~simple_table.chrom.isin(['Y', 'X', 'MT'])]

    # Let's get thoose ranges
    score_range = simple_table.apply(check_strand, axis=1)
    # This is a repeat of the chromosome filtering
    # @TODO check which one to keep
    score_range = clean_score_range(score_range)

    # This is the list of tuples we get back :
    # range_iterator=get_usable_range_withoutX(score_range)

    # We extract the range around all gene's TSS for each chrom
    # And dump them (WPS, FD) in pickle files
    extract_WPS_FD_pos(score_range)

    # We modify simple table to use it later on on the score tables
    # So to add the CAT_geneClass and chrom to the table
    simple_table = simple_table.merge(
        detailed_table[['geneID', 'CAT_geneClass']],
        how='left', left_index=True, on='geneID')

    # We combine all chromosome in one table for WPS and FD
    # And save them as tsv files
    save_table(simple_table, score_range.chrom.unique())
