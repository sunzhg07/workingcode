x=$(wc -l < "$1")
echo $x > tmpf
cat $1 >> tmpf
mv tmpf $1
