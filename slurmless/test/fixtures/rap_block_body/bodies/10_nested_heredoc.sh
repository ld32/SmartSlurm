# doubly/triply nested here-doc: an OUTER heredoc whose body contains a $(...) that
# itself runs an INNER heredoc containing a further $(...). Verbatim capture must
# copy every line byte-for-byte AND the nesting must resolve correctly on the node.
cat <<OUTER > nest.out
before
$(cat <<INNER
mid=$(echo DEEP)
INNER
)
after
OUTER
printf 'ran\n' > ran.txt
