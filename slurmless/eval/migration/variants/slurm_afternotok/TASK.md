# Task: port this SLURM pipeline to SmartSlurm / runAsPipeline

`submit.sh` fetches each word (one `fetch` job per position in `words.txt`). The
fetch of position 5 fails; an `afternotok` cleanup (`recover_5`) then recovers that
word, and `assemble` waits for the good fetches and the recovery, producing
`quote.txt`.

Convert this to an annotated `runAsPipeline` script with the migration-assistant
skill, validate it in the slurmless sandbox, and leave `pipeline.sh` here. The
pipeline's observable output (quote.txt) must match the original's. Keep step
payloads in their original form.
