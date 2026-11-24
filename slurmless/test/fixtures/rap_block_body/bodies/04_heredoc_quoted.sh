# embedded QUOTED here-doc: everything literal, nothing expands
cat <<'INNER' > hdq.out
no expand $(echo X) `echo Y` $HOME
dollar $signs and "quotes" kept
INNER
printf 'ran\n' > ran.txt
