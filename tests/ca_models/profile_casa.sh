#!/bin/bash

# script to get times

pdbs=`ls *.pdb`
for pdb in $pdbs
do		
	(time casa -predict ${pdb} > /dev/null) &> out
	timing=`grep real out | sed 's/m/ /g'`
	nres=`grep CA ${pdb}  | uniq | wc -l`
	echo $pdb $nres $timing
done