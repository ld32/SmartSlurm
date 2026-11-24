# ANSI-C quoting, tabs, unicode, comment line, a non-marker '# @' line
# this is a comment that must survive
# @1 this looks marker-ish but is a comment (space after #)
printf '%b\n' 'tab\tend' > ansi.out
echo $'café \x21' >> ansi.out
printf 'ran\n' > ran.txt
