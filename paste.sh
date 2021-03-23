dir=$1
file=$2
name=${file%.fasta.fas}
paste $dir/features/$name* > /mnt/ramdisk/tmp/$name.features.tsv
