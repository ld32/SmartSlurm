# embedded UNQUOTED here-doc: $(...) runs, \$x stays literal
cat <<INNER > hd.out
sub=$(echo RAN)
literal=\$notexpanded
INNER
printf 'ran\n' > ran.txt
