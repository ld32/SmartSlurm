# case/esac, globs, &&/||, subshells, grouping, semicolons
x=match
case "$x" in
  match) r=yes ;;
  *)     r=no  ;;
esac
( cd / && ls >/dev/null ) || r=err
{ echo "r=$r"; echo "glob chars: *.txt ? [a-z] {a,b}"; } > ctl.out
printf 'ran\n' > ran.txt
