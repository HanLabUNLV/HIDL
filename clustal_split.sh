#!/bin/bash
#this script replaces fasta files with aligned fasta files in the SplitTrees directory
#the first command line argument is the directory over which you want to iterate
#the second command line argument is the directory where you want to store aligned
for direct in "/home/dbarth/Work/GeneTrees/SplitTrees/*/"
	do
	for filename in $direct/*.fasta
		do
		clustalo -i $filename -o $filename.fas >$filename_stdout.txt 2>$filename_stderr.txt 
	done
done
exit


