# single-quoted awk program with braces, $-fields, and a literal $(...) that must NOT run
echo "1 2 3" | awk '{ s=$1+$2+$3; print "sum="s }' > awk.out
echo 'literal $(DO NOT RUN) and $1 field' >> awk.out
printf 'ran\n' > ran.txt
