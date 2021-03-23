file=$1
name=${file%.fasta.fas}
iLearn-protein-basic.py --file $file --method AAC --format tsv --out features/$name.AAC.tsv
iLearn-protein-basic.py --file $file --method GAAC --format tsv --out features/$name.GAAC.tsv
#iLearn-protein-basic.py --file $file --method CTDC --format tsv --out features/$name.CTDC.tsv
#iLearn-protein-basic.py --file $file --method CTDT --format tsv --out features/$name.CTDT.tsv
#iLearn-protein-basic.py --file $file --method CTDD --format tsv --out features/$name.CTDD.tsv
#iLearn-protein-PseKRAAC.py --file $file --method type2 --type 8 --gap_lambda 2  --format tsv --out features/$name.PseKRAAC.tsv
#iLearn-protein-PseKRAAC.py --file $file --method type2 --type 8 --gap_lambda 1  --format tsv --out features/$name.PseKRAAC.g1.tsv
#iLearn-protein-basic.py --file $file --method DisorderC --path examples/predictedProteinProperty
