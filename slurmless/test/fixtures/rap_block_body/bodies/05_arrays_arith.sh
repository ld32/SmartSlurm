# arrays, arithmetic, parameter expansion with special chars
declare -a arr=(one two three)
n=$(( 3 * (2 + 1) ))
echo "count=${#arr[@]} n=$n last=${arr[-1]}" > arr.out
val="abc.def.ghi"; echo "${val%.*} ${val##*.}" >> arr.out
printf 'ran\n' > ran.txt
