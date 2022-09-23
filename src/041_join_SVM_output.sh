#!/bin/bash

for i in $(find data/04_SVC_performances/* | grep tsv); do
	newName=$(echo $i | sed 's/\(.*\)\//\1_/')
	mv -f $i $newName
done
