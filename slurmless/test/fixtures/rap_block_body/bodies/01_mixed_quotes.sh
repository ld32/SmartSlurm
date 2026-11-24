# mixed nested quotes, escaped quotes, literal backslashes
echo "a'b\"c\\d" > q.out
echo 'single with "double" and $notvar and `nobacktick`' >> q.out
printf 'ran\n' > ran.txt
