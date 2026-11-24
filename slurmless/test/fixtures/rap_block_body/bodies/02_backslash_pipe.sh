# backslash line-continuation across a multi-stage pipeline
printf '%s\n' delta alpha charlie \
  | sort \
  | head -1 > first.out
printf 'ran\n' > ran.txt
