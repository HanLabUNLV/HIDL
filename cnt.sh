file=$1
grep '>' $file | wc  | awk '{print $1}' > cnt/$file.cnt
