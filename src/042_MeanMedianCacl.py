#! /usr/bin/python3
import pandas as pd
import argparse
import glob


parser = argparse.ArgumentParser(description='Classification score')
parser.add_argument('-i','--indir', type=str, help='Input dir for SVC-scores files')
parser.add_argument('-o','--outdir', type=str, help='Output dir for grouped SVC-scores files')
args = parser.parse_args()


col_names = [
        "CL_0000037_hematopoietic_stem_cell",
        "CL_0000049_common_myeloid_progenitor",
        "CL_0000071_blood_vessel_endothelial_cell",
        "CL_0000084_T_cell","CL_0000094_granulocyte",
        "CL_0000097_mast_cell",
        "CL_0000127_astrocyte",
        "CL_0000234_phagocyte",
        "CL_0000235_macrophage",
        "CL_0000236_B_cell",
        "CL_0000451_dendritic_cell",
        "CL_0000542_lymphocyte",
        "CL_0000558_reticulocyte",
        "CL_0000576_monocyte",
        "CL_0000623_natural_killer_cell",
        "CL_0000624_CD4_positive_alpha_beta_T_cell",
        "CL_0000625_CD8_positive_alpha_beta_T_cell",
        "CL_0000630_supportive_cell",
        "CL_0000738_leukocyte",
        "CL_0000763_myeloid_cell",
        "CL_0000766_myeloid_leukocyte",
        "CL_0000767_basophil",
        "CL_0000771_eosinophil",
        "CL_0000775_neutrophil",
        "CL_0000784_plasmacytoid_dendritic_cell",
        "CL_0000789_alpha_beta_T_cell",
        "CL_0000798_gamma_delta_T_cell",
        "CL_0000837_hematopoietic_multipotent_progenitor_cell",
        "CL_0000840_immature_conventional_dendritic_cell",
        "CL_0000842_mononuclear_cell",
        "CL_0000860_classical_monocyte",
        "CL_0000893_thymocyte",
        "CL_0000988_hematopoietic_cell",
        "CL_0000990_conventional_dendritic_cell",
        "CL_0002138_endothelial_cell_of_lymphatic_vessel",
        "CL_0002396_CD14_low_CD16_positive_monocyte",
        "CL_0002397_CD14_positive_CD16_positive_monocyte",
        "CL_0002554_fibroblast_of_lymphatic_vessel",
        "UBERON_0000178_blood",
        "UBERON_0001473_lymphatic_vessel",
        "UBERON_0001638_vein",
        "UBERON_0001981_blood_vessel",
        "UBERON_0002193_hemolymphoid_system",
        "UBERON_0002390_hematopoietic_system",
        "UBERON_0002405_immune_system",
        "UBERON_0002465_lymphoid_system",
        "UBERON_0004177_hemopoietic_organ"]

GROUPS=['PFb','NFb','PFNF','PFu','NFu','bu','PFlb','NFlb','blb','PFmid','NFmid','bmid','bNF']
RNA_TYPE=["mrna","lincrna"]
SORTING=["fold","mann_whitney_pval"]
FEAT_SELECTIONS=['Ulzm-vars-peaks',
                 'Ulzm-vars-peak',
                 'Ulzm-var-peaks',
                 'Ulzm-var-peak',
                 'Ulzm-peaks',
                 'Ulzm-peak',
                 'Ulzm-vars',
                 'Ulzm-var',
                 'Ulzs',
                 'Ulz',
                 'Full',
                 'final',
                 'all']
DT=[]
for file in glob.glob(args.indir+'*.tsv'):
    data = pd.read_csv(file, sep='\t', index_col=0)
    filename = file.split('/')[-1].split('.')[0]
    #print(filename)
    All={}
    All["mean"]=data.mean(axis=1)
    All["median"]=data.median(axis=1)
    All["min"]=data.min(axis=1)
    All["max"]=data.max(axis=1)


    lymph_myeloid = data[data.columns.intersection(col_names)]
    lym={}
    lym["mean"]=lymph_myeloid.mean(axis=1)
    lym["median"]=lymph_myeloid.median(axis=1)
    lym["min"]=lymph_myeloid.min(axis=1)
    lym["max"]=lymph_myeloid.max(axis=1)

    no_lymph_myeloid = data.drop(data[data.columns.intersection(col_names)],axis=1)
    nly={}
    nly["mean"]=no_lymph_myeloid.mean(axis=1)
    nly["median"]=no_lymph_myeloid.median(axis=1)
    nly["min"]=no_lymph_myeloid.min(axis=1)
    nly["max"]=no_lymph_myeloid.max(axis=1)

    dict={"all":All,"lymphoMyeloid":lym,"nonLymphoMyeloid":nly}

    dt_list=[]
    for key in dict:
        dt=pd.DataFrame(dict[key])
        dt["ontologyGroup"]=key
        dt_list.append(dt.copy())

    dt=pd.concat(dt_list) 
    dt=dt.reset_index().rename({"index":"performance"}, axis=1) 
    #print(concat_Table.round(2))
    for group in GROUPS:
        if group in filename:
            dt["group"]=group
            break
    for sorting in SORTING:
        if sorting in filename:
            dt["sortingMethod"]=sorting
            break
    for feat in FEAT_SELECTIONS:
        if feat in filename:
            dt["features"]=feat
            break
    for rna in RNA_TYPE:
        if rna in filename:
            dt["rna_type"]=rna
    
    DT.append(dt.copy())

DF=pd.concat(DT)
DF.to_csv(args.outdir+'grouped_performance.tsv', sep='\t')

print('Output files exported')
