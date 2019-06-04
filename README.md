# Homologous Inter-domain Linkers
Repository for the work I'm doing in the Han Lab

The file json_2_fasta uses the python interpreter of json to parse through files downloaded from a library of uniprot/ensembl seqences. This file is in .json format. To convert this odd file format to a more malleable one, .fasta, for alignment or other purposes run this script on the machine with the .json files on it. IMPORTANT NOTE: ONLY WORKS FOR VERSION 90 JSON FILES, NEW FORMAT AS OF V.92

The file phylo_2_fasta gets around this problem by use of the Bio.Phylo package from BioPython. It imports .xml files obtained from the ensembl database and outputs .fasta formats for multiple sequence alignments of the gene tree. Works as of V.92.
