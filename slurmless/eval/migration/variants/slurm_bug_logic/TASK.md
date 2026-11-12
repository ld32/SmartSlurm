# Task: port this SLURM pipeline to SmartSlurm / runAsPipeline

`submit.sh` is a hand-rolled SLURM submission script. It runs one `lower` job per
word position (from `words.txt`), then an `assemble` job that depends on all of
them, producing `quote.txt`.

Using the **migration-assistant** skill, convert this to an annotated
`runAsPipeline` (rAP) script, validate it in the slurmless sandbox, and leave the
rAP `pipeline.sh` (plus any step scripts it needs) in this directory. Keep the step
payloads in their original form.
