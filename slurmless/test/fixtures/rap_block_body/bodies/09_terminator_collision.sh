# a body line that collides with the here-doc terminator must NOT truncate cmd.sh
echo one > coll.out
SMARTSLURM_CMD_EOF
echo two >> coll.out
printf 'ran\n' > ran.txt
