#!/usr/bin/env bash

#Les paramétres à donner en arguments 

infeats=$1              # Features table
outdir=$2               # Output directory
rank=$3                 # The rank of top/bot gene to use
indir_ontology=$4       # Input directory for ontologies tables
lbinfeats=$5            # Input features table 
ubigenes=$6             # Sort by fold, mann_withney_pval 


SVC_script=src/031_ClassificationScoring.py


#Les paramétres "grouping","gene_type", et "feats" sont inclus dans le script et non en arguments.

func() {

	grouping=$1
	feats=$2
	sort_by=$3
	NAME="$grouping"_"$sort_by"_"$feats"

	python3 $SVC_script -r $rank -i $indir_ontology -g $grouping -f $feats -if $infeats -lbf $lbinfeats -u $ubigenes -sb $sort_by -o "$outdir"/"$NAME"_svc_Intermediary -v &> "$outdir"/"$NAME"_svc_Intermediary.log 
	}

####################
 
max_jobs=22

for grouping in PFb NFb PFNF PFu NFu bu PFlb NFlb blb PFmid NFmid bmid bNF ;do
    echo $grouping
    for sort_by in mann_whitney_pval fold ;do
	echo $sort_by
        for feats in Ulz Ulzs Ulzm-var Ulzm-vars Ulzm-peak Ulzm-peaks Ulzm-var-peak Ulzm-vars-peak Ulzm-var-peaks Ulzm-vars-peaks Full final ;do
	    echo $feats
            if (( $(jobs | wc -l) >= max_jobs )); then
                wait -n 
            fi
	    func "$grouping" "$feats" "$sort_by" & 
        done
    done
done


